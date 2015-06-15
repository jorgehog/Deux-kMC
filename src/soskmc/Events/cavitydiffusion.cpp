#include "cavitydiffusion.h"

#include "pressurewall.h"

#include "lammpswriter/lammpswriter.h"

CavityDiffusion::CavityDiffusion(SolidOnSolidSolver &solver,
                                 const double D) :
    SolidOnSolidEvent(solver, "Diffusion"),
    m_D(D)
{
    solver.setDiffusionEvent(*this);
    setDependency(solver.pressureWallEvent());
}

void CavityDiffusion::setupInitialConditions()
{
    const double V = volume();
    const double &h = solver().pressureWallEvent().height();

    uint N = V*exp(solver().gamma()-2*solver().alpha());

    m_particlePositions.set_size(3, N);

    double zMin = solver().heights().min();

    double x0;
    double y0;
    double z0;

    uint n = 0;
    while (n < N)
    {
        do
        {
            x0 = rng.uniform()*(solver().length() - 1);
            y0 = rng.uniform()*(solver().width() - 1);
            z0 = zMin + rng.uniform()*(h - zMin);

        } while(isBlockedPosition(x0, y0, z0));

        m_particlePositions(0, n) = x0;
        m_particlePositions(1, n) = y0;
        m_particlePositions(2, n) = z0;

        n++;
    }

}

void CavityDiffusion::initialize()
{

}

void CavityDiffusion::execute()
{
    _dump(cycle());
}

void CavityDiffusion::reset()
{

    //Has not yet been updated: Represents time spent in current state.
    //    const double &dt = solver().currentTimeStep();
    const double dt = 1.0;

    mat F = zeros(3, nParticles());

    double x1, y1, z1;
    double x0, y0, z0;

    for (uint n = 0; n < nParticles(); ++n)
    {
        x0 = m_particlePositions(0, n);
        y0 = m_particlePositions(1, n);
        z0 = m_particlePositions(2, n);

        BADAssBool(!isBlockedPosition(x0, y0, z0));

        do
        {
            x1 = x0 + sqrt(2*m_D*dt)*rng.normal();// + m_D*F(0, n)*dt;
            y1 = y0 + sqrt(2*m_D*dt)*rng.normal();// + m_D*F(1, n)*dt;
            z1 = z0 + sqrt(2*m_D*dt)*rng.normal();// + m_D*F(2, n)*dt;

            if (solver().dim() == 1)
            {
                y1 = y0;
            }

        } while(isBlockedPosition(x1, y1, z1));

        m_particlePositions(0, n) = x1;
        m_particlePositions(1, n) = y1;
        m_particlePositions(2, n) = z1;

    }

}

void CavityDiffusion::registerHeightChange(const uint x, const uint y, const int value)
{
    if (!hasStarted())
    {
        return;
    }

    (void) x;
    (void) y;
    (void) value;

    //Remove / add the correct particle.

    const double &h = solver().pressureWallEvent().height();
    double zMin = solver().heights().min();

    for (uint n = 0; n < nParticles(); ++n)
    {
        double &x0 = m_particlePositions(0, n);
        double &y0 = m_particlePositions(1, n);
        double &z0 = m_particlePositions(2, n);

        while (isBlockedPosition(x0, y0, z0))
        {
            x0 = rng.uniform()*solver().length();

            if (solver().dim() == 3)
            {
                y0 = rng.uniform()*solver().width();
            }

            z0 = zMin + rng.uniform()*(h - zMin);
        }
    }

}

double CavityDiffusion::localSurfaceSupersaturation(const uint x, const uint y)
{
    (void) x;
    (void) y;

    //do stuff

    return 1.0;
}

double CavityDiffusion::volume() const
{
    return solver().pressureWallEvent().height()*solver().area() - arma::accu(solver().heights());
}

bool CavityDiffusion::isBlockedPosition(const double x, const double y, const double z) const
{
    //DERP: Check periodicity

    bool isOutSideBox_x = (x < 0) || (x > solver().length());

    bool isOutSideBox_y = (y < 0) || (y > solver().width());

    bool isOutSideBox_z = (z > solver().pressureWallEvent().height());

    if (isOutSideBox_x || isOutSideBox_y || isOutSideBox_z)
    {
        return true;
    }

    uint X = uint(round(x));
    uint Y = uint(round(y));

    if (X == solver().length())
    {
        X = solver().length() - 1;
    }

    if (Y == solver().width())
    {
        Y = solver().width() - 1;
    }

    //DERP: collide in x-y plane
    return z < solver().height(X, Y);

}

void CavityDiffusion::_dump(const uint frameNumber) const
{
    const double &h = solver().pressureWallEvent().height();
    const int zMin = solver().heights().min();

    lammpswriter writer(4, "cavitydiff", "/tmp");
    writer.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    writer.initializeNewFile(frameNumber);

    lammpswriter surfacewriter(5, "surfaces", "/tmp");
    surfacewriter.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    surfacewriter.initializeNewFile(frameNumber);


    for (uint n = 0; n < nParticles(); ++n)
    {
        const double &x = m_particlePositions(0, n);
        const double &y = m_particlePositions(1, n);
        const double &z = m_particlePositions(2, n);

        BADAssBool(!isBlockedPosition(x, y, z));

        writer << 0
               << x
               << y
               << z;
    }

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            surfacewriter << 1
                          << x
                          << y
                          << h
                          << 5;

            for (int zSurface = zMin; zSurface <= solver().height(x, y) - 1; ++zSurface)
            {
                surfacewriter << 2
                              << x
                              << y
                              << zSurface
                              << 6;
            }

            surfacewriter << 2
                          << x
                          << y
                          << solver().height(x, y)
                          << solver().nNeighbors(x, y);
        }
    }

    writer.finalize();
    surfacewriter.finalize();



}

