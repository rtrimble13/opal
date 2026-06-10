// Discrete cash dividend schedules and the escrowed-spot adjustment.
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "opal/core/types.hpp"

namespace opal {

struct Dividend {
    double time;    // payment time in years (> 0)
    double amount;  // cash amount per share (> 0)
};

using DividendSchedule = std::vector<Dividend>;

inline void validate(const DividendSchedule& divs, double T) {
    for (const auto& d : divs) {
        require(d.time > 0.0, "dividend: time must be positive");
        require(d.amount >= 0.0, "dividend: amount must be non-negative");
        (void)T;
    }
}

/// PV at `at` of dividends paid strictly after `at` and no later than `until`.
inline double pv_dividends(const DividendSchedule& divs, double r, double at,
                           double until) {
    double pv = 0.0;
    for (const auto& d : divs)
        if (d.time > at + 1e-12 && d.time <= until + 1e-12)
            pv += d.amount * std::exp(-r * (d.time - at));
    return pv;
}

}  // namespace opal
