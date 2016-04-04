#pragma once


class InsertStep : public LatticeEvent
{
public:
    InsertStep(SOSSolver &solver, const uint interval, const uint size) :
        LatticeEvent("InsertStep"),
        m_solver(solver),
        m_interval(interval),
        m_size(size)
    {

    }

private:

    SOSSolver &m_solver;
    const uint m_interval;
    const uint m_size;
    uint m_currentLevel;

    void insertStep()
    {
        m_currentLevel += 1;

        for (uint x = 0; x < m_size; ++x)
        {
            const uint yr = round(sqrt(m_size*m_size - x*x));

            for (uint y = 0; y < yr; ++y)
            {
                solver().setHeight(x, y, m_currentLevel, false);
            }

        }

        m_solver.calculateHeightDependentValues();
    }

    SOSSolver &solver() const
    {
        return m_solver;
    }

    // Event interface
public:
    void execute()
    {

    }

    void initialize()
    {
        m_currentLevel = 0;
        insertStep();
    }

    void reset()
    {
        if (cycle() == 0)
        {
            return;
        }

        if (cycle() % m_interval == 0)
        {
            insertStep();
        }
    }
};

