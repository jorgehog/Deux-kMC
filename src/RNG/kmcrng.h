#pragma once

namespace kMC
{

template<typename seedType>
class KMCRNG
{
public:
    typedef seedType type;

    void initialize(const seedType initialSeed)
    {
        m_initialSeed = initialSeed;

        onInitialize();
    }

    virtual double normal() = 0;
    virtual double uniform() = 0;

    const seedType &initialSeed() const
    {
        return m_initialSeed;
    }

private:
    seedType m_initialSeed;

    virtual void onInitialize() = 0;

};

}
