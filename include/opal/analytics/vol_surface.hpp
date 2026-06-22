// Volatility smile / surface containers built on the SABR model. These hold a
// calibrated parameter set and serve an implied (Black/lognormal) volatility
// for any strike — and, for the surface, any expiry — so a whole strike ladder
// or scenario grid can be priced off one calibrated smile rather than a single
// flat vol. Pairs with the SABR/Heston calibrators.
#pragma once

#include <algorithm>
#include <utility>
#include <vector>

#include "opal/core/types.hpp"
#include "opal/models/sabr.hpp"

namespace opal {

/// A single-expiry SABR smile: a forward, an expiry and the SABR parameters.
/// `vol(strike)` returns the Hagan lognormal implied vol at that strike.
struct SabrSmile {
    double forward;
    double expiry;
    SabrParams params;

    double vol(double strike) const {
        return sabr_lognormal_vol(forward, strike, expiry, params);
    }
};

/// A volatility surface: a stack of SABR smiles at increasing expiries. For a
/// queried (strike, expiry) it interpolates the *total variance* (vol^2 * T)
/// linearly in expiry between the bracketing smiles — the standard
/// arbitrage-friendlier interpolation — and extrapolates flat (in vol) past the
/// first and last expiries.
class VolSurface {
public:
    explicit VolSurface(std::vector<SabrSmile> smiles)
        : smiles_(std::move(smiles)) {
        require(!smiles_.empty(), "vol_surface: at least one smile required");
        std::sort(smiles_.begin(), smiles_.end(),
                  [](const SabrSmile& a, const SabrSmile& b) {
                      return a.expiry < b.expiry;
                  });
    }

    /// Implied lognormal vol at (strike, expiry).
    double vol(double strike, double expiry) const {
        require(expiry > 0.0, "vol_surface: expiry must be positive");
        if (expiry <= smiles_.front().expiry) return smiles_.front().vol(strike);
        if (expiry >= smiles_.back().expiry) return smiles_.back().vol(strike);
        std::size_t hi = 1;
        while (smiles_[hi].expiry < expiry) ++hi;
        const SabrSmile& lo = smiles_[hi - 1];
        const SabrSmile& up = smiles_[hi];
        double vlo = lo.vol(strike), vup = up.vol(strike);
        double w_lo = vlo * vlo * lo.expiry;  // total variance at the lower expiry
        double w_up = vup * vup * up.expiry;
        double frac = (expiry - lo.expiry) / (up.expiry - lo.expiry);
        double w = w_lo + frac * (w_up - w_lo);
        return std::sqrt(w / expiry);
    }

    const std::vector<SabrSmile>& smiles() const { return smiles_; }

private:
    std::vector<SabrSmile> smiles_;
};

}  // namespace opal
