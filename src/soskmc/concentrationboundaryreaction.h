#pragma once

#include "sossolver.h"
#include "../kmcsolver/reaction.h"

#define AXTRANS(callable, n, ...) dim() == 0 ? callable(n, location(), ##__VA_ARGS__) : callable(location(), n, ##__VA_ARGS__)

class concentrationBoundaryReaction : public kMC::Reaction
{
public:
    concentrationBoundaryReaction(const uint dim, const uint orientation, SOSSolver &solver);
    ~concentrationBoundaryReaction();

    double freeBoundaryArea() const;

    double dh(const uint n) const;

    const int &base(const uint n) const;

    const uint &dim() const
    {
        return m_dim;
    }

    const uint &orientation() const
    {
        return m_orientation;
    }

    const uint &location() const
    {
        return m_location;
    }

    const SOSSolver &solver() const
    {
        return m_solver;
    }

    const uint &span() const
    {
        return m_span;
    }

    template<class Callable, class InputType>
    auto axisTransform(Callable &&callable, InputType && input) const -> decltype(callable(input, input))
    {
        if (dim() == 0)
        {
            return callable(std::forward<InputType>(input), location());
        }

        else
        {
            return callable(location(), std::forward<InputType>(input));
        }
    }

private:
    const uint m_dim;
    const uint m_orientation;
    const uint m_location;
    const uint m_span;

    SOSSolver &m_solver;

    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    double rateExpression();
};

