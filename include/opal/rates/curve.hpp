// Discount curve abstraction: flat continuously compounded rate or a
// log-linearly interpolated set of zero rates.
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "opal/core/types.hpp"

namespace opal {

/// Zero curve with continuously compounded zero rates and log-linear
/// discount factor interpolation (flat extrapolation of zero rates).
class DiscountCurve {
public:
    /// Flat curve at continuously compounded rate r.
    explicit DiscountCurve(double r) : flat_rate_(r), flat_(true) {}

    /// Zero curve from (time, zero rate) pillars. Times must be increasing.
    DiscountCurve(std::vector<double> times, std::vector<double> zeros)
        : times_(std::move(times)), zeros_(std::move(zeros)), flat_(false) {
        require(times_.size() == zeros_.size() && !times_.empty(),
                "curve: times and zeros must be same non-zero length");
        for (std::size_t i = 1; i < times_.size(); ++i)
            require(times_[i] > times_[i - 1], "curve: times must be increasing");
        require(times_.front() > 0.0, "curve: times must be positive");
    }

    double zero_rate(double t) const {
        if (flat_) return flat_rate_;
        if (t <= times_.front()) return zeros_.front();
        if (t >= times_.back()) return zeros_.back();
        auto it = std::upper_bound(times_.begin(), times_.end(), t);
        std::size_t i = static_cast<std::size_t>(it - times_.begin());
        double t0 = times_[i - 1], t1 = times_[i];
        double w = (t - t0) / (t1 - t0);
        // Linear in r*t == log-linear in discount factor.
        double rt = (1.0 - w) * zeros_[i - 1] * t0 + w * zeros_[i] * t1;
        return rt / t;
    }

    double discount(double t) const {
        if (t <= 0.0) return 1.0;
        return std::exp(-zero_rate(t) * t);
    }

    /// Simply compounded forward rate over [t1, t2].
    double forward_rate(double t1, double t2) const {
        require(t2 > t1 && t1 >= 0.0, "curve: need t2 > t1 >= 0");
        double tau = t2 - t1;
        return (discount(t1) / discount(t2) - 1.0) / tau;
    }

private:
    double flat_rate_ = 0.0;
    std::vector<double> times_;
    std::vector<double> zeros_;
    bool flat_ = true;
};

}  // namespace opal
