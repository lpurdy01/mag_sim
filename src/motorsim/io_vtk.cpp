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
    int components{1};
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
    std::vector<double> cellBvec;
    std::vector<double> cellHx;
    std::vector<double> cellHy;
    std::vector<double> cellHmag;
    std::vector<double> cellHvec;
    std::vector<double> cellEnergy;
    cellBvec.resize(cells * 3, 0.0);
    if (includeH) {
        cellHx.resize(cells);
        cellHy.resize(cells);
        cellHmag.resize(cells);
        cellHvec.resize(cells * 3, 0.0);
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
            cellBvec[idxCell * 3 + 0] = bxAvg;
            cellBvec[idxCell * 3 + 1] = byAvg;
            cellBvec[idxCell * 3 + 2] = 0.0;

            if (includeH) {
                const double hxAvg = 0.25 * ((*nodeHx)[idx00] + (*nodeHx)[idx10] + (*nodeHx)[idx01] +
                                              (*nodeHx)[idx11]);
                const double hyAvg = 0.25 * ((*nodeHy)[idx00] + (*nodeHy)[idx10] + (*nodeHy)[idx01] +
                                              (*nodeHy)[idx11]);
                cellHx[idxCell] = hxAvg;
                cellHy[idxCell] = hyAvg;
                cellHmag[idxCell] = std::hypot(hxAvg, hyAvg);
                cellHvec[idxCell * 3 + 0] = hxAvg;
                cellHvec[idxCell * 3 + 1] = hyAvg;
                cellHvec[idxCell * 3 + 2] = 0.0;
                if (includeEnergyDensity) {
                    cellEnergy[idxCell] = 0.5 * (bxAvg * hxAvg + byAvg * hyAvg);
                }
            }
        }
    }

    std::vector<DataArrayView> arrays;
    arrays.push_back({"B", &cellBvec, 3});
    arrays.push_back({"Bx", &cellBx, 1});
    arrays.push_back({"By", &cellBy, 1});
    arrays.push_back({"|B|", &cellBmag, 1});
    if (includeH) {
        arrays.push_back({"H", &cellHvec, 3});
        arrays.push_back({"Hx", &cellHx, 1});
        arrays.push_back({"Hy", &cellHy, 1});
        arrays.push_back({"|H|", &cellHmag, 1});
    }
    if (includeEnergyDensity) {
        arrays.push_back({"energy_density", &cellEnergy, 1});
    }

    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open VTK output: " + path);
    }

    const bool littleEndian = isLittleEndian();
    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\""
        << (littleEndian ? "LittleEndian" : "BigEndian") << "\" header_type=\"UInt64\">\n";
    ofs << "  <ImageData WholeExtent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 0\""
        << " Origin=\"" << originX << ' ' << originY << " 0\""
        << " Spacing=\"" << dx << ' ' << dy << " 1\">\n";
    ofs << "    <Piece Extent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 0\">\n";
    ofs << "      <CellData Scalars=\"|B|\" Vectors=\"B\">\n";

    std::uint64_t offset = 0;
    for (const auto& array : arrays) {
        const std::size_t valueCount = array.values->size();
        if (array.components <= 0) {
            throw std::runtime_error("VTK array components must be positive");
        }
        if (valueCount % static_cast<std::size_t>(array.components) != 0) {
            throw std::runtime_error("VTK array value count is not divisible by components for " + array.name);
        }
        const std::uint64_t bytes = static_cast<std::uint64_t>(valueCount) * sizeof(double);
        ofs << "        <DataArray type=\"Float64\" Name=\"" << array.name << "\"";
        if (array.components > 1) {
            ofs << " NumberOfComponents=\"" << array.components << "\"";
        }
        ofs << " format=\"appended\" offset=\"" << offset << "\"/>\n";
        offset += sizeof(std::uint64_t) + bytes;
    }

    ofs << "      </CellData>\n";
    ofs << "    </Piece>\n";
    ofs << "  </ImageData>\n";
    ofs << "  <AppendedData encoding=\"raw\">\n";
    ofs << '_';

    for (const auto& array : arrays) {
        const std::size_t valueCount = array.values->size();
        const std::uint64_t bytes = static_cast<std::uint64_t>(valueCount) * sizeof(double);
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

void write_vtp_outlines(const std::string& path, const std::vector<VtkOutlineLoop>& loops) {
    std::size_t lineCount = 0;
    std::size_t pointTotal = 0;
    for (const auto& loop : loops) {
        if (loop.xs.size() != loop.ys.size()) {
            throw std::invalid_argument("Outline loop coordinate size mismatch");
        }
        if (loop.xs.size() < 2) {
            continue;
        }
        const bool isClosed = std::abs(loop.xs.front() - loop.xs.back()) < 1e-12 &&
                              std::abs(loop.ys.front() - loop.ys.back()) < 1e-12;
        const std::size_t vertices = loop.xs.size() + (isClosed ? 0 : 1);
        pointTotal += vertices;
        ++lineCount;
    }

    std::vector<double> points;
    points.reserve(pointTotal * 3);
    std::vector<long long> connectivity;
    connectivity.reserve(pointTotal);
    std::vector<long long> offsets;
    offsets.reserve(lineCount);
    std::vector<int> kinds;
    kinds.reserve(lineCount);
    std::vector<std::string> labels;
    labels.reserve(lineCount);

    std::size_t pointOffset = 0;
    for (const auto& loop : loops) {
        if (loop.xs.size() != loop.ys.size() || loop.xs.size() < 2) {
            continue;
        }
        const bool isClosed = std::abs(loop.xs.front() - loop.xs.back()) < 1e-12 &&
                              std::abs(loop.ys.front() - loop.ys.back()) < 1e-12;
        const std::size_t baseVertices = loop.xs.size();
        const std::size_t vertices = baseVertices + (isClosed ? 0 : 1);
        for (std::size_t i = 0; i < baseVertices; ++i) {
            points.push_back(loop.xs[i]);
            points.push_back(loop.ys[i]);
            points.push_back(0.0);
            connectivity.push_back(static_cast<long long>(pointOffset + i));
        }
        if (!isClosed) {
            points.push_back(loop.xs.front());
            points.push_back(loop.ys.front());
            points.push_back(0.0);
            connectivity.push_back(static_cast<long long>(pointOffset + baseVertices));
        }
        pointOffset += vertices;
        offsets.push_back(static_cast<long long>(connectivity.size()));
        kinds.push_back(static_cast<int>(loop.kind));
        labels.push_back(loop.label);
    }

    std::ofstream ofs(path);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open VTK outline output: " + path);
    }

    const bool littleEndian = isLittleEndian();
    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\""
        << (littleEndian ? "LittleEndian" : "BigEndian") << "\">\n";
    ofs << "  <PolyData>\n";
    ofs << "    <Piece NumberOfPoints=\"" << points.size() / 3 << "\" NumberOfLines=\""
        << lineCount << "\">\n";
    ofs << "      <Points>\n";
    ofs << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < points.size(); i += 3) {
        ofs << "          " << points[i] << ' ' << points[i + 1] << ' ' << points[i + 2] << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "      </Points>\n";
    ofs << "      <Lines>\n";
    ofs << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    std::size_t cursor = 0;
    for (std::size_t l = 0; l < lineCount; ++l) {
        ofs << "          ";
        const std::size_t end = static_cast<std::size_t>(offsets[l]);
        for (; cursor < end; ++cursor) {
            ofs << connectivity[cursor];
            if (cursor + 1 < end) {
                ofs << ' ';
            }
        }
        ofs << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t l = 0; l < lineCount; ++l) {
        ofs << "          " << offsets[l] << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "      </Lines>\n";
    ofs << "      <CellData Scalars=\"kind\">\n";
    ofs << "        <DataArray type=\"Int32\" Name=\"kind\" format=\"ascii\">\n";
    for (int kind : kinds) {
        ofs << "          " << kind << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "        <DataArray type=\"String\" Name=\"label\" format=\"ascii\">\n";
    for (const auto& label : labels) {
        ofs << "          " << label << "\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "      </CellData>\n";
    ofs << "    </Piece>\n";
    ofs << "  </PolyData>\n";
    ofs << "</VTKFile>\n";
    ofs.flush();
    if (!ofs) {
        throw std::runtime_error("Failed while writing VTK outline output: " + path);
    }
}

}  // namespace motorsim

