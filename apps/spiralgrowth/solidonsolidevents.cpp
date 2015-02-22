#include "solidonsolidevents.h"

SolidOnSolidEvent::~SolidOnSolidEvent()
{

}

void SurfaceSize::execute()
{
    m_localValue = 0;

    for (uint site = 0; site < solver()->length(); ++site)
    {
        m_localValue += 0.5*abs(solver()->height(site) - solver()->height(solver()->rightSite(site)));
    }

    m_localValue /= solver()->length();

    m_sum += solver()->currentTimeStep()*m_localValue;

    setValue(m_sum/(solver()->currentTime() - m_T0));
}


void DumpHeights::execute()
{
    setValue(mean(conv_to<vec>::from(solver()->heights())));

    if (cycle()%100 == 0)
    {
        solver()->heights().save(m_filename);
    }

}
