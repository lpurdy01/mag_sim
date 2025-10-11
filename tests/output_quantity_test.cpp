#include "motorsim/ingest.hpp"
#include "motorsim/io_csv.hpp"
#include "motorsim/io_vtk.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
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
    if (spec.outputs.fieldMaps.size() != 3 || spec.outputs.lineProbes.size() != 1) {
        std::cerr << "Unexpected output counts parsed from scenario\n";
        return 1;
    }
    if (spec.outputs.fieldMaps[0].quantity != "BH" ||
        spec.outputs.fieldMaps[1].quantity != "energy_density" ||
        spec.outputs.fieldMaps[2].format != "vti") {
        std::cerr << "Field map quantities not parsed correctly\n";
        return 1;
    }
    if (spec.outputs.fieldMaps[2].path != "outputs/vtk_map.vti") {
        std::cerr << "Default VTK output path not generated as expected\n";
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

    const fs::path vtkPath = tempDir / "field.vti";
    const std::size_t gridNx = 3;
    const std::size_t gridNy = 3;
    const std::size_t nodeCount = gridNx * gridNy;
    std::vector<double> nodeBx(nodeCount, 1.0);
    std::vector<double> nodeBy(nodeCount, 0.0);
    std::vector<double> nodeHx(nodeCount, 0.5);
    std::vector<double> nodeHy(nodeCount, 0.0);
    try {
        write_vti_field_map(vtkPath.string(), gridNx, gridNy, 0.0, 0.0, 0.01, 0.02, nodeBx, nodeBy,
                            &nodeHx, &nodeHy, true, true);
    } catch (const std::exception& ex) {
        std::cerr << "Failed to write VTK field: " << ex.what() << '\n';
        return 1;
    }

    std::ifstream vti(vtkPath, std::ios::binary);
    if (!vti.is_open()) {
        std::cerr << "Failed to read generated VTK file\n";
        return 1;
    }

    const std::string contents((std::istreambuf_iterator<char>(vti)), std::istreambuf_iterator<char>());
    if (contents.find("Name=\"B\"") == std::string::npos ||
        contents.find("Name=\"Bx\"") == std::string::npos ||
        contents.find("Name=\"H\"") == std::string::npos ||
        contents.find("Name=\"|B|\"") == std::string::npos ||
        contents.find("Name=\"energy_density\"") == std::string::npos) {
        std::cerr << "VTK header missing expected arrays\n";
        return 1;
    }
    const std::size_t appendedPos = contents.find("<AppendedData");
    if (appendedPos == std::string::npos) {
        std::cerr << "VTK file missing AppendedData section\n";
        return 1;
    }
    const std::size_t underscorePos = contents.find('_', appendedPos);
    if (underscorePos == std::string::npos) {
        std::cerr << "VTK file missing appended payload marker\n";
        return 1;
    }
    const char* cursor = contents.data() + underscorePos + 1;
    std::size_t remaining = contents.size() - (underscorePos + 1);
    const auto readBlock = [&](std::size_t expectedCount, std::size_t components) {
        if (remaining < sizeof(std::uint64_t)) {
            throw std::runtime_error("Unexpected end of VTK payload header");
        }
        std::uint64_t byteCount = 0;
        std::memcpy(&byteCount, cursor, sizeof(byteCount));
        cursor += sizeof(byteCount);
        remaining -= sizeof(byteCount);
        const std::size_t expectedBytes = expectedCount * components * sizeof(double);
        if (byteCount != expectedBytes || remaining < byteCount) {
            throw std::runtime_error("VTK payload size mismatch");
        }
        std::vector<double> values(expectedCount * components);
        std::memcpy(values.data(), cursor, expectedBytes);
        cursor += expectedBytes;
        remaining -= expectedBytes;
        return values;
    };

    const std::size_t cellCount = (gridNx - 1) * (gridNy - 1);
    try {
        const auto bVecCells = readBlock(cellCount, 3);
        const auto bxCells = readBlock(cellCount, 1);
        const auto byCells = readBlock(cellCount, 1);
        const auto bmagCells = readBlock(cellCount, 1);
        const auto hVecCells = readBlock(cellCount, 3);
        const auto hxCells = readBlock(cellCount, 1);
        const auto hyCells = readBlock(cellCount, 1);
        const auto hmagCells = readBlock(cellCount, 1);
        const auto energyCells = readBlock(cellCount, 1);

        for (std::size_t i = 0; i < cellCount; ++i) {
            const std::size_t vecIdx = i * 3;
            if (!approxEqual(bVecCells[vecIdx + 0], 1.0) || !approxEqual(bVecCells[vecIdx + 1], 0.0) ||
                !approxEqual(bVecCells[vecIdx + 2], 0.0) || !approxEqual(bxCells[i], 1.0) ||
                !approxEqual(byCells[i], 0.0) || !approxEqual(bmagCells[i], 1.0) ||
                !approxEqual(hVecCells[vecIdx + 0], 0.5) || !approxEqual(hVecCells[vecIdx + 1], 0.0) ||
                !approxEqual(hVecCells[vecIdx + 2], 0.0) || !approxEqual(hxCells[i], 0.5) ||
                !approxEqual(hyCells[i], 0.0) || !approxEqual(hmagCells[i], 0.5) ||
                !approxEqual(energyCells[i], 0.25)) {
                std::cerr << "Unexpected cell values in VTK payload\n";
                return 1;
            }
        }
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 1;
    }

    return 0;
}
