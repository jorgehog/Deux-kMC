#include "reaction.h"

#include <BADAss/badass.h>

using namespace kMC;

Reaction::Reaction():
    m_updateFlag(UNSET_UPDATE_FLAG)
{

}

template<typename T>
void Reaction::registerUpdateFlag(T flag)
{
    static_assert(std::is_scalar<T>::value, "invalid update flag.");

    int iflag = static_cast<int>(flag);

    if (iflag < m_updateFlag)
    {
        m_updateFlag = iflag;
    }
}

template<typename T>
void Reaction::forceUpdateFlag(T flag)
{
    static_assert(std::is_scalar<T>::value, "invalid update flag.");
    m_updateFlag = static_cast<int>(flag);
}
