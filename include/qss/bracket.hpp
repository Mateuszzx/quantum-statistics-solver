#pragma once
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "qss/solver_types.hpp"

namespace qss
{

/**
 * @brief Strategy for finding a bracket.
 */
enum class BracketerType
{
    Expansion,
    Bisection
};

/**
 * @brief Configuration options for the bracketer.
 */
struct BracketerOptions
{
    double expand_factor = 1.6;
    int max_expands = 50;
    double initial_step = 1.0;
};


// ---------- Bracketer ----------
/**
 * @brief Helper class to find a valid bracket for root finding.
 */
class Bracketer
{
public:
    /**
     * @brief Construct a new Bracketer object.
     * 
     * @param type The bracketing strategy to use.
     * @param opts Configuration options.
     */
    explicit Bracketer(BracketerType type, BracketerOptions opts = {}) noexcept
        : type_(type), opts_(opts)
    {
    }

    /**
     * @brief Set the options object.
     * 
     * @param o New options.
     */
    void set_options(const BracketerOptions &o) noexcept
    {
        opts_ = o;
    }
    
    /**
     * @brief Get the current options.
     * 
     * @return const BracketerOptions& 
     */
    const BracketerOptions &options() const noexcept
    {
        return opts_;
    }

    /**
     * @brief Get the bracketing strategy type.
     * 
     * @return BracketerType 
     */
    BracketerType type() const noexcept
    {
        return type_;
    }

    /**
     * @brief Attempt to find a valid bracket for function f starting near [a0, b0].
     * 
     * @param f The function to bracket.
     * @param a0 Initial lower bound guess.
     * @param b0 Initial upper bound guess.
     * @return Bracket The resulting bracket (check .valid flag).
     */
    Bracket bracket(const Function1D &f, double a0, double b0) const;

private:
    BracketerType type_;
    BracketerOptions opts_;
};

} // namespace qss
