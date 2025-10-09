#include "motorsim/ingest.hpp"
#include "motorsim/io_csv.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    const fs::path tempDir = fs::temp_directory_path() / "output_quantity_test";
    std::error_code ec;
    fs::create_directories(tempDir, ec);

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/output_quantity_test.json").lexically_normal();

    ScenarioSpec spec = loadScenarioFromJson(scenarioPath.string());
    if (spec.outputs.fieldMaps.size() != 2 || spec.outputs.lineProbes.size() != 1) {
        std::cerr << "Unexpected output counts parsed from scenario\n";
        return 1;
    }
    if (spec.outputs.fieldMaps[0].quantity != "BH" ||
        spec.outputs.fieldMaps[1].quantity != "energy_density") {
        std::cerr << "Field map quantities not parsed correctly\n";
        return 1;
    }
    if (spec.outputs.lineProbes[0].quantity != "Hmag") {
        std::cerr << "Line probe quantity not parsed correctly\n";
        return 1;
    }

    const std::vector<double> xs{0.0, 0.01, 0.02};
    const std::vector<double> ys{0.0, 0.01, 0.02};
    const std::vector<double> hx{1.0, 2.0, 3.0};
    const std::vector<double> hy{0.5, 0.25, 0.125};
    const std::vector<double> energy{0.4, 0.9, 1.2};
    const fs::path csvPath = tempDir / "field.csv";
    write_csv_field_map(csvPath.string(), xs, ys,
                        {{"Hx", &hx}, {"Hy", &hy}, {"EnergyDensity", &energy}});

    std::ifstream ifs(csvPath);
    if (!ifs.is_open()) {
        std::cerr << "Failed to read generated CSV\n";
        return 1;
    }

    std::string header;
    std::getline(ifs, header);
    if (header != "x,y,Hx,Hy,EnergyDensity") {
        std::cerr << "Unexpected CSV header: " << header << '\n';
        return 1;
    }

    std::string firstRow;
    std::getline(ifs, firstRow);
    std::stringstream ss(firstRow);
    std::string token;
    std::getline(ss, token, ',');
    const double x0 = std::stod(token);
    std::getline(ss, token, ',');
    const double y0 = std::stod(token);
    std::getline(ss, token, ',');
    const double hx0 = std::stod(token);
    std::getline(ss, token, ',');
    const double hy0 = std::stod(token);
    std::getline(ss, token, ',');
    const double e0 = std::stod(token);

    const auto approxEqual = [](double a, double b) {
        return std::abs(a - b) < 1e-12;
    };

    if (!approxEqual(x0, xs[0]) || !approxEqual(y0, ys[0]) || !approxEqual(hx0, hx[0]) ||
        !approxEqual(hy0, hy[0]) || !approxEqual(e0, energy[0])) {
        std::cerr << "CSV first row mismatch: " << firstRow << '\n';
        return 1;
    }

    return 0;
}
