#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace motorsim {

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

}  // namespace motorsim

