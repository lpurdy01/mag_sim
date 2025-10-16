#include "motorsim/circuit.hpp"

#include "motorsim/probes.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace motorsim {
namespace {

constexpr double kTiny = 1e-12;

double safe_orientation(double orientation) {
    if (std::abs(orientation) < kTiny) {
        return 1.0;
    }
    return orientation;
}

}  // namespace

void CircuitSimulator::initialize(const std::vector<ScenarioFrame>& frames) {
    circuits_.clear();
    active_ = false;

    if (frames.empty()) {
        return;
    }

    const ScenarioSpec& baseSpec = frames.front().spec;
    if (baseSpec.circuits.empty()) {
        return;
    }

    circuits_.reserve(baseSpec.circuits.size());

    for (const auto& specCircuit : baseSpec.circuits) {
        CircuitData data{};
        data.id = specCircuit.id;
        data.nodeNames = specCircuit.nodes;
        data.nodeCount = specCircuit.nodes.size();
        data.nodeVarCount = data.nodeCount > 0 ? data.nodeCount - 1 : 0;
        data.voltageVarCount = specCircuit.voltageSources.size();

        data.resistors.reserve(specCircuit.resistors.size());
        for (const auto& res : specCircuit.resistors) {
            if (!(res.resistance > 0.0)) {
                throw std::runtime_error("Circuit '" + specCircuit.id + "' resistor has non-positive resistance");
            }
            ResistorData entry{};
            entry.nodePos = res.nodePos;
            entry.nodeNeg = res.nodeNeg;
            entry.conductance = 1.0 / res.resistance;
            data.resistors.push_back(entry);
        }

        data.voltageSources.reserve(specCircuit.voltageSources.size());
        for (const auto& src : specCircuit.voltageSources) {
            VoltageSourceData entry{};
            entry.nodePos = src.nodePos;
            entry.nodeNeg = src.nodeNeg;
            data.voltageSources.push_back(entry);
        }

        data.inductors.reserve(specCircuit.inductors.size());
        for (const auto& ind : specCircuit.inductors) {
            if (!(ind.inductance > 0.0)) {
                throw std::runtime_error("Circuit '" + specCircuit.id + "' inductor has non-positive inductance");
            }
            InductorState state{};
            state.nodePos = ind.nodePos;
            state.nodeNeg = ind.nodeNeg;
            state.inductance = ind.inductance;
            state.initialCurrent = ind.initialCurrent;
            state.current = ind.initialCurrent;
            state.prevFlux = 0.0;
            state.prevFluxTime = frames.front().time;
            state.lastFlux = 0.0;
            state.lastFluxTime = frames.front().time;
            state.backEmf = 0.0;
            data.inductors.push_back(std::move(state));
        }

        data.coilLinks.reserve(specCircuit.coilLinks.size());
        for (const auto& link : specCircuit.coilLinks) {
            if (link.inductorIndex >= data.inductors.size()) {
                throw std::runtime_error("Circuit '" + specCircuit.id + "' coil_link references unknown inductor");
            }
            if (!(link.turns > 0.0)) {
                throw std::runtime_error("Circuit '" + specCircuit.id + "' coil_link requires positive turns");
            }
            if (link.regionIndex >= baseSpec.currentRegions.size()) {
                throw std::runtime_error("Circuit '" + specCircuit.id + "' coil_link region index out of range");
            }
            CoilLinkData entry{};
            entry.inductorIndex = link.inductorIndex;
            entry.regionIndex = link.regionIndex;
            entry.turns = link.turns;
            data.inductors[entry.inductorIndex].coilLinkIndices.push_back(data.coilLinks.size());
            data.coilLinks.push_back(entry);
        }

        const std::size_t totalUnknowns = data.nodeVarCount + data.voltageVarCount;
        data.lhsTemplate.assign(totalUnknowns * totalUnknowns, 0.0);

        const auto addToMatrix = [&](std::size_t row, std::size_t col, double value) {
            data.lhsTemplate[row * totalUnknowns + col] += value;
        };

        for (const auto& res : data.resistors) {
            const std::size_t pos = res.nodePos;
            const std::size_t neg = res.nodeNeg;
            if (pos != 0 && pos < data.nodeCount) {
                addToMatrix(pos - 1, pos - 1, res.conductance);
            }
            if (neg != 0 && neg < data.nodeCount) {
                addToMatrix(neg - 1, neg - 1, res.conductance);
            }
            if (pos != 0 && neg != 0 && pos < data.nodeCount && neg < data.nodeCount) {
                addToMatrix(pos - 1, neg - 1, -res.conductance);
                addToMatrix(neg - 1, pos - 1, -res.conductance);
            }
        }

        for (std::size_t vsIdx = 0; vsIdx < data.voltageSources.size(); ++vsIdx) {
            const auto& src = data.voltageSources[vsIdx];
            const std::size_t col = data.nodeVarCount + vsIdx;
            if (src.nodePos != 0 && src.nodePos < data.nodeCount) {
                const std::size_t row = src.nodePos - 1;
                addToMatrix(row, col, 1.0);
                addToMatrix(col, row, 1.0);
            }
            if (src.nodeNeg != 0 && src.nodeNeg < data.nodeCount) {
                const std::size_t row = src.nodeNeg - 1;
                addToMatrix(row, col, -1.0);
                addToMatrix(col, row, -1.0);
            }
        }

        circuits_.push_back(std::move(data));
    }

    const auto& regions = frames.front().spec.currentRegions;
    for (auto& circuit : circuits_) {
        for (auto& state : circuit.inductors) {
            if (!state.coilLinkIndices.empty() && !regions.empty()) {
                double accum = 0.0;
                std::size_t count = 0;
                for (std::size_t linkIndex : state.coilLinkIndices) {
                    if (linkIndex >= circuit.coilLinks.size()) {
                        continue;
                    }
                    const auto& link = circuit.coilLinks[linkIndex];
                    if (link.regionIndex >= regions.size()) {
                        continue;
                    }
                    const auto& region = regions[link.regionIndex];
                    const double orientation = safe_orientation(region.orientation);
                    accum += region.current / orientation;
                    ++count;
                }
                if (count > 0) {
                    state.current = accum / static_cast<double>(count);
                }
            }
        }
    }

    active_ = !circuits_.empty();
}

bool CircuitSimulator::solve_dense_linear_system(std::vector<double>& matrix,
                                                 std::vector<double>& rhs,
                                                 std::size_t dim) {
    if (dim == 0) {
        return true;
    }

    for (std::size_t pivot = 0; pivot < dim; ++pivot) {
        std::size_t pivotRow = pivot;
        double maxVal = std::abs(matrix[pivot * dim + pivot]);
        for (std::size_t row = pivot + 1; row < dim; ++row) {
            const double value = std::abs(matrix[row * dim + pivot]);
            if (value > maxVal) {
                maxVal = value;
                pivotRow = row;
            }
        }
        if (maxVal < kTiny) {
            return false;
        }
        if (pivotRow != pivot) {
            for (std::size_t col = 0; col < dim; ++col) {
                std::swap(matrix[pivot * dim + col], matrix[pivotRow * dim + col]);
            }
            std::swap(rhs[pivot], rhs[pivotRow]);
        }
        const double pivotVal = matrix[pivot * dim + pivot];
        for (std::size_t row = pivot + 1; row < dim; ++row) {
            const double factor = matrix[row * dim + pivot] / pivotVal;
            if (std::abs(factor) < kTiny) {
                continue;
            }
            for (std::size_t col = pivot; col < dim; ++col) {
                matrix[row * dim + col] -= factor * matrix[pivot * dim + col];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    for (std::size_t row = dim; row-- > 0;) {
        double sum = rhs[row];
        for (std::size_t col = row + 1; col < dim; ++col) {
            sum -= matrix[row * dim + col] * rhs[col];
        }
        rhs[row] = sum / matrix[row * dim + row];
    }

    return true;
}

void CircuitSimulator::compute_rhs(const CircuitData& circuit,
                                   const std::vector<double>& inductorCurrents,
                                   std::vector<double>& rhs) const {
    std::fill(rhs.begin(), rhs.end(), 0.0);
    for (std::size_t idx = 0; idx < circuit.inductors.size(); ++idx) {
        const auto& state = circuit.inductors[idx];
        const double current = inductorCurrents[idx];
        if (state.nodePos != 0 && state.nodePos < circuit.nodeCount) {
            rhs[state.nodePos - 1] -= current;
        }
        if (state.nodeNeg != 0 && state.nodeNeg < circuit.nodeCount) {
            rhs[state.nodeNeg - 1] += current;
        }
    }
}

void CircuitSimulator::evaluate_derivatives(const CircuitData& circuit,
                                            const std::vector<double>& inductorCurrents,
                                            const std::vector<double>& voltages,
                                            std::vector<double>& derivs) const {
    const std::size_t totalUnknowns = circuit.nodeVarCount + circuit.voltageVarCount;
    if (totalUnknowns == 0) {
        std::fill(derivs.begin(), derivs.end(), 0.0);
        return;
    }

    std::vector<double> matrix = circuit.lhsTemplate;
    std::vector<double> rhs(totalUnknowns, 0.0);
    compute_rhs(circuit, inductorCurrents, rhs);
    for (std::size_t idx = 0; idx < circuit.voltageVarCount; ++idx) {
        rhs[circuit.nodeVarCount + idx] = (idx < voltages.size()) ? voltages[idx] : 0.0;
    }

    if (!solve_dense_linear_system(matrix, rhs, totalUnknowns)) {
        throw std::runtime_error("Circuit solver encountered singular system");
    }

    for (std::size_t idx = 0; idx < circuit.inductors.size(); ++idx) {
        const auto& state = circuit.inductors[idx];
        const double va = (state.nodePos == 0 || state.nodePos >= circuit.nodeCount)
                              ? 0.0
                              : rhs[state.nodePos - 1];
        const double vb = (state.nodeNeg == 0 || state.nodeNeg >= circuit.nodeCount)
                              ? 0.0
                              : rhs[state.nodeNeg - 1];
        const double branchVoltage = va - vb;
        derivs[idx] = (branchVoltage - state.backEmf) / state.inductance;
    }
}

void CircuitSimulator::apply_currents(ScenarioFrame& frame) {
    if (!active_ || frame.spec.currentRegions.empty()) {
        return;
    }

    for (std::size_t circuitIndex = 0; circuitIndex < circuits_.size(); ++circuitIndex) {
        const auto& circuit = circuits_[circuitIndex];
        for (const auto& link : circuit.coilLinks) {
            if (link.regionIndex >= frame.spec.currentRegions.size()) {
                continue;
            }
            auto& region = frame.spec.currentRegions[link.regionIndex];
            const double orientation = safe_orientation(region.orientation);
            const double current = circuit.inductors[link.inductorIndex].current;
            region.current = orientation * current;
        }
    }
}

void CircuitSimulator::record_solved_frame(const ScenarioFrame& frame, const Grid2D& grid) {
    if (!active_) {
        return;
    }
    if (grid.Bx.empty() || grid.By.empty()) {
        return;
    }

    for (auto& circuit : circuits_) {
        for (std::size_t idx = 0; idx < circuit.inductors.size(); ++idx) {
            auto& state = circuit.inductors[idx];
            double linkage = 0.0;
            for (std::size_t linkIndex : state.coilLinkIndices) {
                if (linkIndex >= circuit.coilLinks.size()) {
                    continue;
                }
                const auto& link = circuit.coilLinks[linkIndex];
                if (link.regionIndex >= frame.spec.currentRegions.size()) {
                    continue;
                }
                const auto& region = frame.spec.currentRegions[link.regionIndex];
                const double flux = integrate_polygon_flux_component(
                    grid, frame.spec.originX, frame.spec.originY, frame.spec.dx, frame.spec.dy, region.xs,
                    region.ys, region.min_x, region.max_x, region.min_y, region.max_y, FluxComponent::Bmag);
                linkage += link.turns * safe_orientation(region.orientation) * flux;
            }
            state.prevFlux = state.lastFlux;
            state.prevFluxTime = state.lastFluxTime;
            state.lastFlux = linkage;
            state.lastFluxTime = frame.time;
            const double dt = state.lastFluxTime - state.prevFluxTime;
            if (std::abs(dt) > kTiny) {
                state.backEmf = -(state.lastFlux - state.prevFlux) / dt;
            }
        }
    }
}

void CircuitSimulator::prepare_next_frame(std::size_t /*frameIndex*/,
                                          const ScenarioFrame& frame,
                                          const ScenarioFrame* nextFrame) {
    if (!active_ || !nextFrame) {
        return;
    }

    const double t0 = frame.time;
    const double t1 = nextFrame->time;
    const double dt = t1 - t0;
    if (std::abs(dt) < kTiny) {
        return;
    }

    for (std::size_t circuitIndex = 0; circuitIndex < circuits_.size(); ++circuitIndex) {
        auto& circuit = circuits_[circuitIndex];
        if (circuit.inductors.empty()) {
            continue;
        }

        std::vector<double> currents(circuit.inductors.size(), 0.0);
        for (std::size_t idx = 0; idx < circuit.inductors.size(); ++idx) {
            currents[idx] = circuit.inductors[idx].current;
        }

        const auto& currentCircuitSpec = frame.spec.circuits[circuitIndex];
        const auto& nextCircuitSpec = nextFrame->spec.circuits[circuitIndex];

        std::vector<double> v0(circuit.voltageVarCount, 0.0);
        std::vector<double> v1(circuit.voltageVarCount, 0.0);
        for (std::size_t idx = 0; idx < circuit.voltageVarCount; ++idx) {
            v0[idx] = currentCircuitSpec.voltageSources[idx].value;
            v1[idx] = nextCircuitSpec.voltageSources[idx].value;
        }

        auto sampleVoltages = [&](double time, std::vector<double>& out) {
            const double alpha = std::clamp((time - t0) / dt, 0.0, 1.0);
            for (std::size_t idx = 0; idx < out.size(); ++idx) {
                out[idx] = v0[idx] + alpha * (v1[idx] - v0[idx]);
            }
        };

        std::vector<double> voltages(circuit.voltageVarCount, 0.0);
        std::vector<double> k1(circuit.inductors.size(), 0.0);
        std::vector<double> k2(circuit.inductors.size(), 0.0);
        std::vector<double> k3(circuit.inductors.size(), 0.0);
        std::vector<double> k4(circuit.inductors.size(), 0.0);
        std::vector<double> tempCurrents = currents;

        sampleVoltages(t0, voltages);
        evaluate_derivatives(circuit, currents, voltages, k1);

        for (std::size_t idx = 0; idx < tempCurrents.size(); ++idx) {
            tempCurrents[idx] = currents[idx] + 0.5 * dt * k1[idx];
        }
        sampleVoltages(t0 + 0.5 * dt, voltages);
        evaluate_derivatives(circuit, tempCurrents, voltages, k2);

        for (std::size_t idx = 0; idx < tempCurrents.size(); ++idx) {
            tempCurrents[idx] = currents[idx] + 0.5 * dt * k2[idx];
        }
        sampleVoltages(t0 + 0.5 * dt, voltages);
        evaluate_derivatives(circuit, tempCurrents, voltages, k3);

        for (std::size_t idx = 0; idx < tempCurrents.size(); ++idx) {
            tempCurrents[idx] = currents[idx] + dt * k3[idx];
        }
        sampleVoltages(t1, voltages);
        evaluate_derivatives(circuit, tempCurrents, voltages, k4);

        for (std::size_t idx = 0; idx < currents.size(); ++idx) {
            currents[idx] += (dt / 6.0) * (k1[idx] + 2.0 * k2[idx] + 2.0 * k3[idx] + k4[idx]);
            circuit.inductors[idx].current = currents[idx];
        }
    }
}

std::optional<double> CircuitSimulator::inductor_current(const std::string& circuitId,
                                                          std::size_t inductorIndex) const {
    for (const auto& circuit : circuits_) {
        if (circuit.id == circuitId) {
            if (inductorIndex < circuit.inductors.size()) {
                return circuit.inductors[inductorIndex].current;
            }
            return std::nullopt;
        }
    }
    return std::nullopt;
}

}  // namespace motorsim

