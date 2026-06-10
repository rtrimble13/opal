// Output formatting (aligned tables, JSON, CSV) for the opal CLI.
#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace opal::cli {

enum class OutputFormat { Table, Json, Csv };

inline OutputFormat parse_output(const std::string& s) {
    if (s == "json") return OutputFormat::Json;
    if (s == "csv") return OutputFormat::Csv;
    return OutputFormat::Table;
}

inline std::string fmt_num(double x, int prec = 6) {
    if (std::isnan(x)) return "-";
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec) << x;
    return os.str();
}

/// Key-value report: aligned two-column table or a JSON object.
class Report {
public:
    void add(const std::string& key, double value, int prec = 6) {
        rows_.emplace_back(key, fmt_num(value, prec));
        numeric_.emplace_back(key, value);
    }
    void add(const std::string& key, const std::string& value) {
        rows_.emplace_back(key, value);
        strings_.emplace_back(key, value);
    }

    void print(OutputFormat f, const std::string& title = "") const {
        if (f == OutputFormat::Json) {
            print_json();
            return;
        }
        if (!title.empty()) {
            std::cout << title << "\n";
            std::cout << std::string(title.size(), '-') << "\n";
        }
        std::size_t w = 0;
        for (auto& [k, v] : rows_) w = std::max(w, k.size());
        for (auto& [k, v] : rows_)
            std::cout << "  " << std::left << std::setw(static_cast<int>(w) + 2)
                      << k << v << "\n";
    }

private:
    void print_json() const {
        std::cout << "{";
        bool first = true;
        for (auto& [k, v] : strings_) {
            std::cout << (first ? "" : ", ") << "\"" << k << "\": \"" << v << "\"";
            first = false;
        }
        for (auto& [k, v] : numeric_) {
            std::cout << (first ? "" : ", ") << "\"" << k << "\": ";
            if (std::isnan(v))
                std::cout << "null";
            else
                std::cout << v;
            first = false;
        }
        std::cout << "}\n";
    }

    std::vector<std::pair<std::string, std::string>> rows_;
    std::vector<std::pair<std::string, double>> numeric_;
    std::vector<std::pair<std::string, std::string>> strings_;
};

/// Matrix/grid output: aligned table, CSV, or JSON array-of-rows.
class Table {
public:
    explicit Table(std::vector<std::string> headers) : headers_(std::move(headers)) {}

    void add_row(const std::vector<std::string>& row) { rows_.push_back(row); }

    void print(OutputFormat f) const {
        if (f == OutputFormat::Csv) {
            print_sep(",");
            return;
        }
        if (f == OutputFormat::Json) {
            print_json();
            return;
        }
        std::vector<std::size_t> w(headers_.size(), 0);
        for (std::size_t i = 0; i < headers_.size(); ++i) w[i] = headers_[i].size();
        for (auto& r : rows_)
            for (std::size_t i = 0; i < r.size() && i < w.size(); ++i)
                w[i] = std::max(w[i], r[i].size());
        auto line = [&](const std::vector<std::string>& r) {
            for (std::size_t i = 0; i < r.size(); ++i)
                std::cout << (i ? "  " : "") << std::right
                          << std::setw(static_cast<int>(w[i])) << r[i];
            std::cout << "\n";
        };
        line(headers_);
        std::size_t total = 0;
        for (auto x : w) total += x + 2;
        std::cout << std::string(total > 2 ? total - 2 : total, '-') << "\n";
        for (auto& r : rows_) line(r);
    }

private:
    void print_sep(const std::string& sep) const {
        auto line = [&](const std::vector<std::string>& r) {
            for (std::size_t i = 0; i < r.size(); ++i)
                std::cout << (i ? sep : "") << r[i];
            std::cout << "\n";
        };
        line(headers_);
        for (auto& r : rows_) line(r);
    }

    void print_json() const {
        std::cout << "[";
        for (std::size_t j = 0; j < rows_.size(); ++j) {
            std::cout << (j ? ",\n " : "") << "{";
            for (std::size_t i = 0; i < rows_[j].size() && i < headers_.size(); ++i) {
                std::cout << (i ? ", " : "") << "\"" << headers_[i] << "\": ";
                // Emit numbers unquoted when they parse as such.
                char* end = nullptr;
                std::strtod(rows_[j][i].c_str(), &end);
                bool is_num = end && *end == '\0' && !rows_[j][i].empty();
                if (is_num)
                    std::cout << rows_[j][i];
                else
                    std::cout << "\"" << rows_[j][i] << "\"";
            }
            std::cout << "}";
        }
        std::cout << "]\n";
    }

    std::vector<std::string> headers_;
    std::vector<std::vector<std::string>> rows_;
};

}  // namespace opal::cli
