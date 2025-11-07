#pragma once
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "solver_types.hpp"

namespace qss
{

enum class BracketerType
{
    Expansion,
    Bisection
};

struct BracketerOptions
{
    double expand_factor = 1.6;
    int max_expands = 50;
    double initial_step = 1.0;
};


// ---------- Bracketer ----------
class Bracketer
{
public:
    explicit Bracketer(BracketerType type, BracketerOptions opts = {}) noexcept
        : type_(type), opts_(opts)
    {
    }

    void set_options(const BracketerOptions &o) noexcept
    {
        opts_ = o;
    }
    const BracketerOptions &options() const noexcept
    {
        return opts_;
    }
    BracketerType type() const noexcept
    {
        return type_;
    }

    Bracket bracket(const Function1D &f, double a0, double b0) const;

private:
    BracketerType type_;
    BracketerOptions opts_;
};

} // namespace qss
