#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace motorsim {

struct VtkOutlineLoop {
    enum class Kind : int {
        Domain = 0,
        Material = 1,
        Magnet = 2,
        Wire = 3,
        CurrentRegion = 4,
    };

    Kind kind{Kind::Domain};
    std::string label;
    std::string groupLabel;
    std::vector<double> xs;
    std::vector<double> ys;
};

// Writes a VTK ImageData (.vti) field map using cell-centred values.
// nodeBx/nodeBy should contain nx * ny samples laid out in row-major order.
// nodeHx/nodeHy may be null when H data is not required.
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
                         bool includeEnergyDensity);

// Emits a VTK PolyData (.vtp) file containing closed polyline outlines that
// represent scenario geometry (domain boundary, material polygons, magnets,
// wires, etc.). Each loop is tagged with a categorical `kind` and a `label`
// string stored as cell data for convenient filtering in ParaView. A
// companion CSV stores the label together with any `groupLabel` metadata so
// groupings (e.g. rotor assemblies) can be reconstructed without relying on
// VTK string arrays.
void write_vtp_outlines(const std::string& path, const std::vector<VtkOutlineLoop>& loops);

struct PvdDataSet {
    double time{0.0};
    std::string file;
};

void write_pvd_series(const std::string& path, const std::vector<PvdDataSet>& datasets);

}  // namespace motorsim

