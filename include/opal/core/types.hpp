// Core enumerations and result types shared across the library.
#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

namespace opal {

enum class OptionType { Call, Put };
enum class ExerciseStyle { European, American };
enum class BarrierType { DownIn, DownOut, UpIn, UpOut };
enum class AverageType { Arithmetic, Geometric };
enum class StrikeStyle { Fixed, Floating };

inline double type_sign(OptionType t) { return t == OptionType::Call ? 1.0 : -1.0; }

inline std::string to_string(OptionType t) {
    return t == OptionType::Call ? "call" : "put";
}

inline std::string to_string(BarrierType b) {
    switch (b) {
        case BarrierType::DownIn: return "down-in";
        case BarrierType::DownOut: return "down-out";
        case BarrierType::UpIn: return "up-in";
        case BarrierType::UpOut: return "up-out";
    }
    return "?";
}

inline bool is_knock_in(BarrierType b) {
    return b == BarrierType::DownIn || b == BarrierType::UpIn;
}
inline bool is_down_barrier(BarrierType b) {
    return b == BarrierType::DownIn || b == BarrierType::DownOut;
}

/// First and second order sensitivities of an option price.
/// Conventions: vega per 1.00 of vol (divide by 100 for per vol point),
/// theta per year (divide by 365 for per calendar day),
/// rho per 1.00 of rate.
struct Greeks {
    double price = 0.0;
    double delta = 0.0;
    double gamma = 0.0;
    double vega = 0.0;
    double theta = 0.0;
    double rho = 0.0;
    // Higher order (populated by analytic models).
    double vanna = std::nan("");  // d delta / d vol
    double volga = std::nan("");  // d vega / d vol
    double charm = std::nan("");  // d delta / d time
};

/// Result of a Monte Carlo pricing run.
struct McResult {
    double price = 0.0;
    double std_error = 0.0;
    std::size_t paths = 0;
};

inline void require(bool cond, const std::string& msg) {
    if (!cond) throw std::invalid_argument(msg);
}

}  // namespace opal
