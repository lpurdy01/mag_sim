#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "motorsim/ingest.hpp"
#include "motorsim/probes.hpp"

namespace motorsim {

class MechanicalSimulator {
public:
    MechanicalSimulator() = default;

    struct RotorSample {
        double time{0.0};
        double angleRad{0.0};
        double omega{0.0};
        double torque{0.0};
    };

    void initialize(const ScenarioSpec& baseSpec, const std::vector<ScenarioFrame>& frames);

    [[nodiscard]] bool is_active() const { return active_; }

    void apply_state(ScenarioFrame& frame) const;

    void handle_solved_frame(const ScenarioFrame& frame,
                             const std::unordered_map<std::string, StressTensorResult>& torqueResults,
                             ScenarioFrame* nextFrame);

    [[nodiscard]] std::optional<double> rotor_angle_deg(const std::string& name) const;
    [[nodiscard]] std::optional<double> rotor_speed_rad_s(const std::string& name) const;

    [[nodiscard]] const std::unordered_map<std::string, std::vector<RotorSample>>& history() const {
        return history_;
    }

private:
    struct RotorState {
        std::string name;
        std::size_t rotorIndex{0};
        double inertia{0.0};
        double damping{0.0};
        double loadTorque{0.0};
        double angleRad{0.0};
        double omega{0.0};
        double lastTorque{0.0};
        std::string torqueProbeId;
        bool hasTorqueProbe{false};
        bool missingTorqueWarning{false};
    };

    void apply_state_to_frame(ScenarioFrame& frame, const RotorState& state) const;
    void rotateRotorGeometry(const ScenarioSpec::Rotor& rotor,
                             double angleRad,
                             ScenarioSpec& workingSpec) const;

    const ScenarioSpec* baseSpec_{nullptr};
    std::vector<RotorState> rotors_;
    std::unordered_map<std::string, std::vector<RotorSample>> history_;
    bool active_{false};
};

}  // namespace motorsim
