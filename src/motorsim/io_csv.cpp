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
                         const std::vector<FieldColumnView>& columns) {
    const std::size_t n = x.size();
    if (y.size() != n) {
        throw std::invalid_argument("write_csv_field_map: mismatched vector sizes");
    }
    for (const auto& column : columns) {
        if (column.values == nullptr) {
            throw std::invalid_argument("write_csv_field_map: null column pointer for '" +
                                        column.header + "'");
        }
        if (column.values->size() != n) {
            throw std::invalid_argument(
                "write_csv_field_map: column size mismatch for '" + column.header + "'");
        }
    }

    std::ofstream ofs(path);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open CSV output: " + path);
    }

    ofs << "x,y";
    for (const auto& column : columns) {
        ofs << ',' << column.header;
    }
    ofs << '\n';

    for (std::size_t i = 0; i < n; ++i) {
        ofs << x[i] << ',' << y[i];
        for (const auto& column : columns) {
            ofs << ',' << (*(column.values))[i];
        }
        ofs << '\n';
    }
}

}  // namespace motorsim
