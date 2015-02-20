#pragma once

#include <sys/types.h>
#include <limits>
#include <sstream>
#include <set>
#include <functional>

#include <libconfig_utils/libconfig_utils.h>


namespace kMC
{

class Reaction
{
public:

    Reaction();

    virtual ~Reaction() {}

    virtual bool isAllowed() const = 0;

    virtual void execute() = 0;

    template<typename T>
    void registerUpdateFlag(T flag);

    const int & updateFlag() const
    {
        return m_updateFlag;
    }

    void resetUpdateFlag()
    {
        m_updateFlag = UNSET_UPDATE_FLAG;
    }

    template<typename T>
    void forceUpdateFlag(const T flag);

    //! Update flags are given in the order such that the minimum of the flag set is the
    //! triumphant flag.
    enum AllUpdateFlags
    {
        UNSET_UPDATE_FLAG = std::numeric_limits<int>::max(),
        defaultUpdateFlag = 0
    };

private:

    int m_updateFlag;

};

}

