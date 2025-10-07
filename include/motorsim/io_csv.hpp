// filename: io_csv.hpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#pragma once

#include <string>
#include <vector>

namespace motorsim {

void write_csv_line_profile(const std::string& path,
                            const std::vector<double>& x,
                            const std::vector<double>& y,
                            const std::vector<double>& v);

void write_csv_field_map(const std::string& path,
                         const std::vector<double>& x,
                         const std::vector<double>& y,
                         const std::vector<double>& bx,
                         const std::vector<double>& by,
                         const std::vector<double>& bmag);

}  // namespace motorsim
