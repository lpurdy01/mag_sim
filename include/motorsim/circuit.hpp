#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>
#include <limits>

#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"

namespace motorsim {

class MechanicalSimulator;

class CircuitSimulator {
public:
    CircuitSimulator() = default;

    explicit CircuitSimulator(const std::vector<ScenarioFrame>& frames) { initialize(frames); }

    void initialize(const std::vector<ScenarioFrame>& frames);

    [[nodiscard]] bool is_active() const { return active_; }

    void update_for_frame(ScenarioFrame& frame, const MechanicalSimulator* mechanical);

    void apply_currents(ScenarioFrame& frame);

    void record_solved_frame(const ScenarioFrame& frame, const Grid2D& grid);

    void prepare_next_frame(std::size_t frameIndex,
                            const ScenarioFrame& frame,
                            const ScenarioFrame* nextFrame);

    [[nodiscard]] std::optional<double> inductor_current(const std::string& circuitId,
                                                         std::size_t inductorIndex) const;

    struct TraceSample {
        double time{0.0};
        double current{0.0};
        double effectiveCurrent{0.0};
        double ampereTurns{0.0};
        double orientation{0.0};
        double backEmf{0.0};
    };

    struct TraceMetadata {
        std::string circuitId;
        std::string coilId;
        std::string regionId;
        double turns{0.0};
    };

    [[nodiscard]] const std::unordered_map<std::string, std::vector<TraceSample>>& history() const { return history_; }
    [[nodiscard]] const std::unordered_map<std::string, TraceMetadata>& trace_metadata() const {
        return traceMetadata_;
    }

private:
    struct ResistorData {
        std::size_t nodePos{0};
        std::size_t nodeNeg{0};
        double conductance{0.0};
    };

    struct VoltageSourceData {
        std::size_t nodePos{0};
        std::size_t nodeNeg{0};
    };

    struct CoilLinkData {
        std::size_t inductorIndex{0};
        std::size_t regionIndex{0};
        double turns{0.0};
        double baseOrientation{1.0};
        double currentOrientation{1.0};
        std::string coilId;
        std::string regionId;
        std::string traceKey;
        struct CommutatorSegment {
            double startDeg{0.0};
            double endDeg{0.0};
            double orientation{1.0};
        };
        struct CommutatorData {
            bool active{false};
            std::size_t rotorIndex{std::numeric_limits<std::size_t>::max()};
            std::string rotorName;
            double defaultOrientation{1.0};
            std::vector<CommutatorSegment> segments;
        } commutator;
    };

    struct InductorState {
        std::size_t nodePos{0};
        std::size_t nodeNeg{0};
        double inductance{0.0};
        double current{0.0};
        double initialCurrent{0.0};
        double prevFlux{0.0};
        double prevFluxTime{0.0};
        double lastFlux{0.0};
        double lastFluxTime{0.0};
        double backEmf{0.0};
        std::vector<std::size_t> coilLinkIndices;
    };

    struct CircuitData {
        std::string id;
        std::vector<std::string> nodeNames;
        std::vector<ResistorData> resistors;
        std::vector<VoltageSourceData> voltageSources;
        std::vector<InductorState> inductors;
        std::vector<CoilLinkData> coilLinks;
        std::vector<double> lhsTemplate;
        std::size_t nodeCount{0};
        std::size_t nodeVarCount{0};
        std::size_t voltageVarCount{0};
    };

    std::vector<CircuitData> circuits_;
    const ScenarioSpec* baseSpec_{nullptr};
    bool active_{false};
    std::unordered_map<std::string, std::vector<TraceSample>> history_;
    std::unordered_map<std::string, TraceMetadata> traceMetadata_;

    static bool solve_dense_linear_system(std::vector<double>& matrix,
                                          std::vector<double>& rhs,
                                          std::size_t dim);

    static double evaluate_commutator(const CoilLinkData& link, double angleDeg);

    void compute_rhs(const CircuitData& circuit,
                     const std::vector<double>& inductorCurrents,
                     std::vector<double>& rhs) const;

    void evaluate_derivatives(const CircuitData& circuit,
                              const std::vector<double>& inductorCurrents,
                              const std::vector<double>& voltages,
                              std::vector<double>& derivs) const;
};

}  // namespace motorsim

