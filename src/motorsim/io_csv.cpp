// filename: io_csv.cpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#include "motorsim/io_csv.hpp"

#include <fstream>
#include <stdexcept>

namespace motorsim {

void write_csv_line_profile(const std::string& path,
                            const std::vector<double>& x,
                            const std::vector<double>& y,
                            const std::vector<double>& v) {
    if (x.size() != y.size() || x.size() != v.size()) {
        throw std::invalid_argument("write_csv_line_profile: mismatched vector sizes");
    }

    std::ofstream ofs(path);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open CSV output: " + path);
    }

    ofs << "x,y,value\n";
    for (std::size_t i = 0; i < x.size(); ++i) {
        ofs << x[i] << ',' << y[i] << ',' << v[i] << '\n';
    }
}

void write_csv_field_map(const std::string& path,
                         const std::vector<double>& x,
                         const std::vector<double>& y,
                         const std::vector<double>& bx,
                         const std::vector<double>& by,
                         const std::vector<double>& bmag) {
    const std::size_t n = x.size();
    if (y.size() != n || bx.size() != n || by.size() != n || bmag.size() != n) {
        throw std::invalid_argument("write_csv_field_map: mismatched vector sizes");
    }

    std::ofstream ofs(path);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open CSV output: " + path);
    }

    ofs << "x,y,Bx,By,Bmag\n";
    for (std::size_t i = 0; i < n; ++i) {
        ofs << x[i] << ',' << y[i] << ',' << bx[i] << ',' << by[i] << ',' << bmag[i]
            << '\n';
    }
}

}  // namespace motorsim
