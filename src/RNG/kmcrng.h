#pragma once

namespace kMC
{

template<typename seedType>
class KMCRNG
{
public:
    seedType m_initialSeed;

    void initialize(const seedType initialSeed)
    {
        m_initialSeed = initialSeed;

        onInitialize();
    }

    virtual double normal() = 0;
    virtual double uniform() = 0;

private:
    virtual void onInitialize() = 0;

};

}
