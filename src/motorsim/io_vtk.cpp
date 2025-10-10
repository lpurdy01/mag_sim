#include "motorsim/io_vtk.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace motorsim {
namespace {

bool isLittleEndian() {
    const std::uint16_t value = 1;
    return reinterpret_cast<const std::uint8_t*>(&value)[0] == 1;
}

std::size_t cellCount(std::size_t nx, std::size_t ny) {
    if (nx < 2 || ny < 2) {
        throw std::invalid_argument("VTK export requires at least a 2x2 grid");
    }
    return (nx - 1) * (ny - 1);
}

struct DataArrayView {
    std::string name;
    const std::vector<double>* values{nullptr};
};

}  // namespace

void write_vti_field_map(const std::string& path,
                         std::size_t nx,
                         std::size_t ny,
                         double originX,
                         double originY,
                         double dx,
                         double dy,
                         const std::vector<double>& nodeBx,
                         const std::vector<double>& nodeBy,
                         const std::vector<double>* nodeHx,
                         const std::vector<double>* nodeHy,
                         bool includeH,
                         bool includeEnergyDensity) {
    const std::size_t nodeCount = nx * ny;
    if (nodeBx.size() != nodeCount || nodeBy.size() != nodeCount) {
        throw std::invalid_argument("VTK export requires node-aligned B field arrays");
    }
    if (includeH) {
        if (nodeHx == nullptr || nodeHy == nullptr) {
            throw std::invalid_argument("VTK export requested H field but pointers were null");
        }
        if (nodeHx->size() != nodeCount || nodeHy->size() != nodeCount) {
            throw std::invalid_argument("VTK export requires node-aligned H field arrays");
        }
    }
    if (includeEnergyDensity && !includeH) {
        throw std::invalid_argument("Energy density output requires H field data");
    }

    const std::size_t cells = cellCount(nx, ny);
    const std::size_t cellNx = nx - 1;

    std::vector<double> cellBx(cells);
    std::vector<double> cellBy(cells);
    std::vector<double> cellBmag(cells);
    std::vector<double> cellHx;
    std::vector<double> cellHy;
    std::vector<double> cellHmag;
    std::vector<double> cellEnergy;
    if (includeH) {
        cellHx.resize(cells);
        cellHy.resize(cells);
        cellHmag.resize(cells);
    }
    if (includeEnergyDensity) {
        cellEnergy.resize(cells);
    }

    for (std::size_t j = 0; j + 1 < ny; ++j) {
        for (std::size_t i = 0; i + 1 < nx; ++i) {
            const std::size_t idxCell = j * cellNx + i;
            const std::size_t idx00 = j * nx + i;
            const std::size_t idx10 = j * nx + (i + 1);
            const std::size_t idx01 = (j + 1) * nx + i;
            const std::size_t idx11 = (j + 1) * nx + (i + 1);

            const double bxAvg = 0.25 *
                                 (nodeBx[idx00] + nodeBx[idx10] + nodeBx[idx01] + nodeBx[idx11]);
            const double byAvg = 0.25 *
                                 (nodeBy[idx00] + nodeBy[idx10] + nodeBy[idx01] + nodeBy[idx11]);
            cellBx[idxCell] = bxAvg;
            cellBy[idxCell] = byAvg;
            cellBmag[idxCell] = std::hypot(bxAvg, byAvg);

            if (includeH) {
                const double hxAvg = 0.25 * ((*nodeHx)[idx00] + (*nodeHx)[idx10] + (*nodeHx)[idx01] +
                                              (*nodeHx)[idx11]);
                const double hyAvg = 0.25 * ((*nodeHy)[idx00] + (*nodeHy)[idx10] + (*nodeHy)[idx01] +
                                              (*nodeHy)[idx11]);
                cellHx[idxCell] = hxAvg;
                cellHy[idxCell] = hyAvg;
                cellHmag[idxCell] = std::hypot(hxAvg, hyAvg);
                if (includeEnergyDensity) {
                    cellEnergy[idxCell] = 0.5 * (bxAvg * hxAvg + byAvg * hyAvg);
                }
            }
        }
    }

    std::vector<DataArrayView> arrays;
    arrays.push_back({"Bx", &cellBx});
    arrays.push_back({"By", &cellBy});
    arrays.push_back({"|B|", &cellBmag});
    if (includeH) {
        arrays.push_back({"Hx", &cellHx});
        arrays.push_back({"Hy", &cellHy});
        arrays.push_back({"|H|", &cellHmag});
    }
    if (includeEnergyDensity) {
        arrays.push_back({"energy_density", &cellEnergy});
    }

    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open VTK output: " + path);
    }

    const bool littleEndian = isLittleEndian();
    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\""
        << (littleEndian ? "LittleEndian" : "BigEndian") << "\" header_type=\"UInt64\"\">\n";
    ofs << "  <ImageData WholeExtent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 0\""
        << " Origin=\"" << originX << ' ' << originY << " 0\""
        << " Spacing=\"" << dx << ' ' << dy << " 1\">\n";
    ofs << "    <Piece Extent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 0\">\n";
    ofs << "      <CellData Scalars=\"|B|\">\n";

    std::uint64_t offset = 0;
    for (const auto& array : arrays) {
        const std::size_t count = array.values->size();
        const std::uint64_t bytes = static_cast<std::uint64_t>(count) * sizeof(double);
        ofs << "        <DataArray type=\"Float64\" Name=\"" << array.name
            << "\" format=\"appended\" offset=\"" << offset << "\"/>\n";
        offset += sizeof(std::uint64_t) + bytes;
    }

    ofs << "      </CellData>\n";
    ofs << "    </Piece>\n";
    ofs << "  </ImageData>\n";
    ofs << "  <AppendedData encoding=\"raw\">\n";
    ofs << '_';

    for (const auto& array : arrays) {
        const std::size_t count = array.values->size();
        const std::uint64_t bytes = static_cast<std::uint64_t>(count) * sizeof(double);
        ofs.write(reinterpret_cast<const char*>(&bytes), sizeof(bytes));
        ofs.write(reinterpret_cast<const char*>(array.values->data()), static_cast<std::streamsize>(bytes));
    }

    ofs << "\n";
    ofs << "  </AppendedData>\n";
    ofs << "</VTKFile>\n";
    ofs.flush();
    if (!ofs) {
        throw std::runtime_error("Failed while writing VTK output: " + path);
    }
}

}  // namespace motorsim

