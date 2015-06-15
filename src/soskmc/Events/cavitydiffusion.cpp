#include "cavitydiffusion.h"

#include "pressurewall.h"

#include "lammpswriter/lammpswriter.h"

#include "../solidonsolidreaction.h"

CavityDiffusion::CavityDiffusion(SolidOnSolidSolver &solver,
                                 const double D,
                                 const double dt) :
    SolidOnSolidEvent(solver, "Diffusion"),
    m_D0(D),
    m_D(D*exp(-solver.gamma())/(2*solver.dim())),
    m_dt0(dt),
    m_dt(dt*exp(solver.gamma())/solver.area())
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

    m_F = zeros(3, N);

    m_localRates.set_size(nParticles(), solver().length(), solver().width());

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
    const double &dtFull = solver().currentTimeStep();

    const uint N = dtFull/m_dt;

    cout << dtFull << " " << m_dt << " " << N << endl;

    for (uint i = 0; i < N; ++i)
    {
        diffuse(m_dt);
    }

    diffuse(dtFull - N*m_dt);
}

void CavityDiffusion::registerHeightChange(const uint x, const uint y, const int value)
{
    if (!hasStarted())
    {
        return;
    }


    //Remove / add the correct particle.

    if (value == -1)
    {

        vector<double> localRatesForSite(nParticles());

        double Rtot = 0;
        for (uint n = 0; n < nParticles(); ++n)
        {
            localRatesForSite.at(n) = m_localRates(n, x, y);
            Rtot += localRatesForSite.at(n);
        }

        double R = rng.uniform()*Rtot;

        uint N = binarySearchForInterval(R, localRatesForSite);

        cout << "should remove particle " << N << endl;

    }

    else
    {
        uint ns = 6 - solver().nNeighbors(x, y);

        uint R = (uint)(rng.uniform()*ns);



    }

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

    DiffusionDeposition *r;
    for (uint _x = 0; _x < solver().length(); ++_x)
    {
        for (uint _y = 0; _y < solver().width(); ++_y)
        {
            r = &solver().reaction(_x, _y);
            r->setDepositionRate(r->calculateDepositionRate());
        }
    }
}

double CavityDiffusion::localSurfaceSupersaturation(const uint x, const uint y)
{

    double R = 0;

    const int z = solver().height(x, y) + 1;

    for (uint n = 0; n < nParticles(); ++n)
    {
        const double dxSquared = pow((double)x - m_particlePositions(0, n), 2);
        const double dySquared = pow((double)y - m_particlePositions(1, n), 2);
        const double dzSquared = pow((double)z - m_particlePositions(2, n), 2);

        const double Rn = 1./sqrt(dxSquared + dySquared + dzSquared);
        R += Rn;

        m_localRates(n, x, y) = Rn;
    }

    return R*exp(-solver().gamma());

    return 1.0;
}

double CavityDiffusion::volume() const
{
    //sum (h_l - h_i) - 1
    return (solver().pressureWallEvent().height() - 1)*solver().area() - arma::accu(solver().heights());
}

bool CavityDiffusion::isBlockedPosition(const double x, const double y, const double z) const
{
    //DERP: Check periodicity

    bool isOutSideBox_x = (x < 0) || (x > solver().length());

    bool isOutSideBox_y = (y < 0) || (y > solver().width());

    bool isOutSideBox_z = (z > solver().pressureWallEvent().height() - 0.5);

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
    return z < solver().height(X, Y) + 0.5;

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

void CavityDiffusion::diffuse(const double dt)
{
    if (dt < 0 || fabs(dt) < 1E-15)
    {
        return;
    }

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
            x1 = x0 + sqrt(2*m_D*dt)*rng.normal() + m_D*m_F(0, n)*dt;
            y1 = y0 + sqrt(2*m_D*dt)*rng.normal() + m_D*m_F(1, n)*dt;
            z1 = z0 + sqrt(2*m_D*dt)*rng.normal() + m_D*m_F(2, n)*dt;

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



