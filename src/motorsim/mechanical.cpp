#include "motorsim/mechanical.hpp"

#include <cmath>
#include <iostream>
#include <limits>

namespace motorsim {

namespace {
constexpr double kPi = 3.14159265358979323846;
}

void MechanicalSimulator::initialize(const ScenarioSpec& baseSpec,
                                     const std::vector<ScenarioFrame>& /*frames*/) {
    baseSpec_ = &baseSpec;
    rotors_.clear();
    active_ = false;

    if (!baseSpec.mechanical.has_value()) {
        return;
    }

    const ScenarioSpec::MechanicalSystem& mechanical = *baseSpec.mechanical;
    if (mechanical.rotors.empty()) {
        return;
    }

    rotors_.reserve(mechanical.rotors.size());
    for (const auto& rotorConfig : mechanical.rotors) {
        if (rotorConfig.rotorIndex >= baseSpec.rotors.size()) {
            continue;
        }
        RotorState state{};
        state.name = !rotorConfig.name.empty() ? rotorConfig.name : baseSpec.rotors[rotorConfig.rotorIndex].name;
        state.rotorIndex = rotorConfig.rotorIndex;
        state.inertia = rotorConfig.inertia;
        state.damping = rotorConfig.damping;
        state.loadTorque = rotorConfig.loadTorque;
        state.angleRad = rotorConfig.hasInitialAngle ? rotorConfig.initialAngleDeg * kPi / 180.0 : 0.0;
        state.omega = rotorConfig.hasInitialSpeed ? rotorConfig.initialSpeedRadPerSec : 0.0;
        state.torqueProbeId = rotorConfig.torqueProbeId;
        state.hasTorqueProbe = rotorConfig.hasTorqueProbe;
        state.lastTorque = 0.0;
        rotors_.push_back(state);
    }

    if (!rotors_.empty()) {
        active_ = true;
    }
}

void MechanicalSimulator::apply_state(ScenarioFrame& frame) const {
    if (!active_ || !baseSpec_) {
        return;
    }
    for (const auto& state : rotors_) {
        apply_state_to_frame(frame, state);
    }
}

void MechanicalSimulator::handle_solved_frame(
    const ScenarioFrame& frame,
    const std::unordered_map<std::string, StressTensorResult>& torqueResults,
    ScenarioFrame* nextFrame) {
    if (!active_ || !baseSpec_) {
        return;
    }

    double dt = 0.0;
    if (nextFrame) {
        dt = nextFrame->time - frame.time;
    }

    for (auto& state : rotors_) {
        double torque = state.hasTorqueProbe ? state.lastTorque : 0.0;
        if (state.hasTorqueProbe) {
            const auto it = torqueResults.find(state.torqueProbeId);
            if (it != torqueResults.end()) {
                torque = it->second.torqueZ;
                state.lastTorque = torque;
                state.missingTorqueWarning = false;
            } else if (!state.missingTorqueWarning) {
                std::cerr << "MechanicalSimulator: torque probe '" << state.torqueProbeId
                          << "' not found for rotor '" << state.name << "'\n";
                state.missingTorqueWarning = true;
            }
        }

        if (!(dt > 0.0) || !(state.inertia > 0.0)) {
            continue;
        }

        const auto accel = [&](double omega) {
            const double netTorque = torque - state.loadTorque - state.damping * omega;
            return netTorque / state.inertia;
        };

        const double k1Omega = accel(state.omega);
        const double k1Theta = state.omega;

        const double omegaK2 = state.omega + 0.5 * dt * k1Omega;
        const double k2Omega = accel(omegaK2);
        const double k2Theta = omegaK2;

        const double omegaK3 = state.omega + 0.5 * dt * k2Omega;
        const double k3Omega = accel(omegaK3);
        const double k3Theta = omegaK3;

        const double omegaK4 = state.omega + dt * k3Omega;
        const double k4Omega = accel(omegaK4);
        const double k4Theta = omegaK4;

        state.omega += (dt / 6.0) * (k1Omega + 2.0 * k2Omega + 2.0 * k3Omega + k4Omega);
        state.angleRad += (dt / 6.0) * (k1Theta + 2.0 * k2Theta + 2.0 * k3Theta + k4Theta);

        if (nextFrame) {
            apply_state_to_frame(*nextFrame, state);
        }
    }
}

std::optional<double> MechanicalSimulator::rotor_angle_deg(const std::string& name) const {
    for (const auto& state : rotors_) {
        if (state.name == name) {
            return state.angleRad * 180.0 / kPi;
        }
    }
    return std::nullopt;
}

std::optional<double> MechanicalSimulator::rotor_speed_rad_s(const std::string& name) const {
    for (const auto& state : rotors_) {
        if (state.name == name) {
            return state.omega;
        }
    }
    return std::nullopt;
}

void MechanicalSimulator::apply_state_to_frame(ScenarioFrame& frame, const RotorState& state) const {
    if (!baseSpec_ || state.rotorIndex >= baseSpec_->rotors.size() ||
        state.rotorIndex >= frame.spec.rotors.size()) {
        return;
    }
    rotateRotorGeometry(baseSpec_->rotors[state.rotorIndex], state.angleRad, frame.spec);
}

void MechanicalSimulator::rotateRotorGeometry(const ScenarioSpec::Rotor& rotor,
                                              double angleRad,
                                              ScenarioSpec& workingSpec) const {
    if (!baseSpec_) {
        return;
    }

    const double cosA = std::cos(angleRad);
    const double sinA = std::sin(angleRad);

    for (std::size_t idx : rotor.polygonIndices) {
        if (idx >= baseSpec_->polygons.size() || idx >= workingSpec.polygons.size()) {
            continue;
        }
        const auto& basePoly = baseSpec_->polygons[idx];
        auto& poly = workingSpec.polygons[idx];
        if (basePoly.xs.size() != basePoly.ys.size() || basePoly.xs.empty()) {
            continue;
        }
        poly.xs.resize(basePoly.xs.size());
        poly.ys.resize(basePoly.ys.size());
        double minX = std::numeric_limits<double>::infinity();
        double maxX = -std::numeric_limits<double>::infinity();
        double minY = std::numeric_limits<double>::infinity();
        double maxY = -std::numeric_limits<double>::infinity();
        for (std::size_t v = 0; v < basePoly.xs.size(); ++v) {
            const double dx = basePoly.xs[v] - rotor.pivotX;
            const double dy = basePoly.ys[v] - rotor.pivotY;
            const double rx = rotor.pivotX + dx * cosA - dy * sinA;
            const double ry = rotor.pivotY + dx * sinA + dy * cosA;
            poly.xs[v] = rx;
            poly.ys[v] = ry;
            minX = std::min(minX, rx);
            maxX = std::max(maxX, rx);
            minY = std::min(minY, ry);
            maxY = std::max(maxY, ry);
        }
        poly.min_x = minX;
        poly.max_x = maxX;
        poly.min_y = minY;
        poly.max_y = maxY;
    }

    for (std::size_t idx : rotor.magnetIndices) {
        if (idx >= baseSpec_->magnetRegions.size() || idx >= workingSpec.magnetRegions.size()) {
            continue;
        }
        const auto& baseMagnet = baseSpec_->magnetRegions[idx];
        auto& magnet = workingSpec.magnetRegions[idx];

        std::vector<double> baseXs = baseMagnet.xs;
        std::vector<double> baseYs = baseMagnet.ys;
        if (baseXs.size() != baseYs.size() || baseXs.size() < 3) {
            baseXs = {baseMagnet.min_x, baseMagnet.max_x, baseMagnet.max_x, baseMagnet.min_x};
            baseYs = {baseMagnet.min_y, baseMagnet.min_y, baseMagnet.max_y, baseMagnet.max_y};
        }
        magnet.xs.resize(baseXs.size());
        magnet.ys.resize(baseYs.size());
        double minX = std::numeric_limits<double>::infinity();
        double maxX = -std::numeric_limits<double>::infinity();
        double minY = std::numeric_limits<double>::infinity();
        double maxY = -std::numeric_limits<double>::infinity();
        for (std::size_t v = 0; v < baseXs.size(); ++v) {
            const double dx = baseXs[v] - rotor.pivotX;
            const double dy = baseYs[v] - rotor.pivotY;
            const double rx = rotor.pivotX + dx * cosA - dy * sinA;
            const double ry = rotor.pivotY + dx * sinA + dy * cosA;
            magnet.xs[v] = rx;
            magnet.ys[v] = ry;
            minX = std::min(minX, rx);
            maxX = std::max(maxX, rx);
            minY = std::min(minY, ry);
            maxY = std::max(maxY, ry);
        }
        magnet.min_x = minX;
        magnet.max_x = maxX;
        magnet.min_y = minY;
        magnet.max_y = maxY;

        const double mx = baseMagnet.Mx;
        const double my = baseMagnet.My;
        magnet.Mx = mx * cosA - my * sinA;
        magnet.My = mx * sinA + my * cosA;
    }

    for (std::size_t idx : rotor.wireIndices) {
        if (idx >= baseSpec_->wires.size() || idx >= workingSpec.wires.size()) {
            continue;
        }
        const auto& baseWire = baseSpec_->wires[idx];
        auto& wire = workingSpec.wires[idx];
        const double dx = baseWire.x - rotor.pivotX;
        const double dy = baseWire.y - rotor.pivotY;
        wire.x = rotor.pivotX + dx * cosA - dy * sinA;
        wire.y = rotor.pivotY + dx * sinA + dy * cosA;
    }
}

}  // namespace motorsim
