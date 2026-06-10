// Command line argument parsing for the opal CLI.
#pragma once

#include <cctype>
#include <chrono>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace opal::cli {

class Args {
public:
    Args(int argc, char** argv) {
        for (int i = 1; i < argc; ++i) tokens_.emplace_back(argv[i]);
        if (!tokens_.empty() && tokens_[0][0] != '-') {
            command_ = tokens_[0];
            tokens_.erase(tokens_.begin());
        }
        parse();
    }

    const std::string& command() const { return command_; }

    bool has(const std::string& key) const { return kv_.count(key) > 0; }

    std::string get_str(const std::string& key, const std::string& def = "") const {
        auto it = kv_.find(key);
        return it == kv_.end() ? def : it->second;
    }

    double get_num(const std::string& key, double def) const {
        auto it = kv_.find(key);
        if (it == kv_.end()) return def;
        return parse_num(key, it->second);
    }

    double require_num(const std::string& key) const {
        auto it = kv_.find(key);
        if (it == kv_.end())
            throw std::runtime_error("missing required option --" + key);
        return parse_num(key, it->second);
    }

    long get_int(const std::string& key, long def) const {
        auto it = kv_.find(key);
        if (it == kv_.end()) return def;
        try {
            return std::stol(it->second);
        } catch (...) {
            throw std::runtime_error("option --" + key + ": expected an integer, got '" +
                                     it->second + "'");
        }
    }

    /// Expiry: plain year fraction ("0.5") or a date "YYYY-MM-DD" converted
    /// to an ACT/365F year fraction from today.
    double get_expiry(const std::string& key = "expiry") const {
        auto it = kv_.find(key);
        if (it == kv_.end())
            throw std::runtime_error("missing required option --" + key);
        const std::string& v = it->second;
        if (v.size() == 10 && v[4] == '-' && v[7] == '-') {
            std::tm tm = {};
            tm.tm_year = std::stoi(v.substr(0, 4)) - 1900;
            tm.tm_mon = std::stoi(v.substr(5, 2)) - 1;
            tm.tm_mday = std::stoi(v.substr(8, 2));
            tm.tm_hour = 12;
            std::time_t target = timegm(&tm);
            std::time_t now = std::time(nullptr);
            double days = double(target - now) / 86400.0;
            if (days <= 0.0)
                throw std::runtime_error("--" + key + ": date " + v +
                                         " is not in the future");
            return days / 365.0;
        }
        return parse_num(key, v);
    }

    /// Range option "lo:hi:step" -> vector of values (inclusive of hi).
    std::vector<double> get_range(const std::string& key) const {
        auto it = kv_.find(key);
        if (it == kv_.end())
            throw std::runtime_error("missing required option --" + key);
        const std::string& v = it->second;
        auto p1 = v.find(':');
        auto p2 = v.find(':', p1 + 1);
        if (p1 == std::string::npos || p2 == std::string::npos)
            throw std::runtime_error("option --" + key +
                                     ": expected format lo:hi:step, got '" + v + "'");
        double lo = std::stod(v.substr(0, p1));
        double hi = std::stod(v.substr(p1 + 1, p2 - p1 - 1));
        double step = std::stod(v.substr(p2 + 1));
        if (step <= 0.0 || hi < lo)
            throw std::runtime_error("option --" + key + ": invalid range '" + v + "'");
        std::vector<double> out;
        for (double x = lo; x <= hi + 1e-9; x += step) out.push_back(x);
        return out;
    }

    bool flag(const std::string& key) const {
        auto it = kv_.find(key);
        return it != kv_.end() && (it->second.empty() || it->second == "true");
    }

private:
    static double parse_num(const std::string& key, const std::string& v) {
        // Accept "5%" as 0.05 for rate-like inputs.
        try {
            if (!v.empty() && v.back() == '%')
                return std::stod(v.substr(0, v.size() - 1)) / 100.0;
            return std::stod(v);
        } catch (...) {
            throw std::runtime_error("option --" + key + ": expected a number, got '" +
                                     v + "'");
        }
    }

    void parse() {
        static const std::map<std::string, std::string> shorts = {
            {"-S", "spot"},    {"-K", "strike"}, {"-T", "expiry"},
            {"-r", "rate"},    {"-q", "div"},    {"-v", "vol"},
            {"-t", "type"},    {"-i", "instrument"}, {"-m", "method"},
            {"-H", "barrier"}, {"-o", "output"}, {"-n", "notional"},
        };
        for (std::size_t i = 0; i < tokens_.size(); ++i) {
            std::string tok = tokens_[i];
            std::string key;
            if (tok.rfind("--", 0) == 0) {
                key = tok.substr(2);
            } else if (tok.rfind("-", 0) == 0 && tok.size() >= 2 &&
                       !std::isdigit(static_cast<unsigned char>(tok[1])) && tok[1] != '.') {
                auto it = shorts.find(tok);
                if (it == shorts.end())
                    throw std::runtime_error("unknown option " + tok);
                key = it->second;
            } else {
                throw std::runtime_error("unexpected argument '" + tok + "'");
            }
            // "--key=value" or "--key value" or bare flag
            auto eq = key.find('=');
            if (eq != std::string::npos) {
                kv_[key.substr(0, eq)] = key.substr(eq + 1);
            } else if (i + 1 < tokens_.size() &&
                       (tokens_[i + 1][0] != '-' ||
                        (tokens_[i + 1].size() > 1 &&
                         (std::isdigit(static_cast<unsigned char>(tokens_[i + 1][1])) ||
                          tokens_[i + 1][1] == '.')))) {
                kv_[key] = tokens_[++i];
            } else {
                kv_[key] = "";
            }
        }
    }

    std::string command_;
    std::vector<std::string> tokens_;
    std::map<std::string, std::string> kv_;
};

}  // namespace opal::cli
