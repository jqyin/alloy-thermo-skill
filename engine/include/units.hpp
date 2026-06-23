// units.hpp -- physical constants and unit conversions.
//
// The pair-interaction energies J_ij are read in some energy unit (Ry, eV, or
// meV; selected at runtime via the config) and converted internally to meV.
// Temperatures are read in Kelvin and converted to meV via Boltzmann's
// constant.  Working in a single internal unit (meV) keeps the Metropolis and
// replica-exchange acceptance expressions dimensionless and unambiguous.
#ifndef ALLOY_MC_UNITS_HPP
#define ALLOY_MC_UNITS_HPP

#include <string>
#include <stdexcept>

namespace alloy {

// Boltzmann constant in meV / K (CODATA 2018: 8.617333262e-2 meV/K).
inline constexpr double kB_meV_per_K = 8.617333262e-2;

// Energy-unit -> meV conversion factors.
inline constexpr double Ry_to_meV  = 13605.693122994;  // 1 Rydberg in meV
inline constexpr double eV_to_meV  = 1000.0;
inline constexpr double meV_to_meV = 1.0;

// Resolve a human-readable energy-unit name ("ry", "ev", "mev") to its
// conversion factor into meV.  Used to interpret the J_ij coupling values.
inline double energy_unit_to_meV(const std::string& unit) {
    if (unit == "ry" || unit == "Ry" || unit == "rydberg") return Ry_to_meV;
    if (unit == "ev" || unit == "eV")                        return eV_to_meV;
    if (unit == "mev" || unit == "meV")                      return meV_to_meV;
    throw std::invalid_argument("unknown energy_unit '" + unit +
                                "' (expected one of: ry, ev, mev)");
}

}  // namespace alloy

#endif  // ALLOY_MC_UNITS_HPP
