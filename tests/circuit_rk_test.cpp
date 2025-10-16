#include "motorsim/circuit.hpp"
#include "motorsim/ingest.hpp"

#include <cassert>
#include <cmath>

namespace {

motorsim::ScenarioSpec make_single_branch_spec(double resistance, double inductance, double voltage) {
    motorsim::ScenarioSpec spec{};
    spec.version = "0.2";
    spec.Lx = 0.1;
    spec.Ly = 0.1;
    spec.nx = 3;
    spec.ny = 3;
    spec.dx = spec.Lx / static_cast<double>(spec.nx - 1);
    spec.dy = spec.Ly / static_cast<double>(spec.ny - 1);

    motorsim::ScenarioSpec::Circuit circuit{};
    circuit.id = "test";
    circuit.nodes = {"gnd", "p"};

    motorsim::ScenarioSpec::Circuit::Resistor resistor{};
    resistor.id = "R";
    resistor.nodePos = 1;
    resistor.nodeNeg = 0;
    resistor.resistance = resistance;
    circuit.resistors.push_back(resistor);

    motorsim::ScenarioSpec::Circuit::Inductor inductor{};
    inductor.id = "L";
    inductor.nodePos = 1;
    inductor.nodeNeg = 0;
    inductor.inductance = inductance;
    circuit.inductors.push_back(inductor);

    motorsim::ScenarioSpec::Circuit::VoltageSource source{};
    source.id = "V";
    source.nodePos = 1;
    source.nodeNeg = 0;
    source.value = voltage;
    circuit.voltageSources.push_back(source);

    spec.circuits.push_back(circuit);
    return spec;
}

}  // namespace

int main() {
    const double resistance = 2.0;
    const double inductance = 0.5;
    const double voltage = 10.0;
    const double dt = 0.01;

    motorsim::ScenarioFrame frame0{};
    frame0.index = 0;
    frame0.time = 0.0;
    frame0.spec = make_single_branch_spec(resistance, inductance, voltage);

    motorsim::ScenarioFrame frame1{};
    frame1.index = 1;
    frame1.time = dt;
    frame1.spec = frame0.spec;
    frame1.spec.circuits[0].voltageSources[0].value = voltage;

    std::vector<motorsim::ScenarioFrame> frames{frame0, frame1};
    motorsim::CircuitSimulator simulator(frames);
    assert(simulator.is_active());

    simulator.prepare_next_frame(0, frames[0], &frames[1]);

    const double expected = (voltage / resistance) * (1.0 - std::exp(-resistance * dt / inductance));
    const auto current = simulator.inductor_current("test", 0);
    assert(current.has_value());
    const double error = std::abs(current.value() - expected);
    const double tolerance = 5e-3;
    assert(error < tolerance);
    return 0;
}

