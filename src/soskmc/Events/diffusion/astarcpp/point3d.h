#pragma once

#include <sstream>

using std::string;
using std::stringstream;

namespace Tests
{
/// <summary>
/// Author: Roy Triesscheijn (http://www.roy-t.nl)
/// Point3D class mimics some of the Microsoft.Xna.Framework.Vector3
/// but uses Int32's instead of floats.
/// </summary>
class Point3D
{
public:
    int X;
    int Y;
    int Z;

    Point3D(int X, int Y, int Z) :
        X(X),
        Y(Y),
        Z(Z)
    {
        refCounter++;
    }

    ~Point3D()
    {
        refCounter--;
    }

    int GetDistanceSquared(const Point3D* point) const
    {
        return GetDistanceSquared(point->X, point->Y, point->Z);
    }

    int GetDistanceSquared(const int x1, const int y1, const int z1) const
    {
        int dx = X - x1;
        int dy = Y - y1;
        int dz = Z - z1;
        return (dx * dx) + (dy * dy) + (dz * dz);
    }

    bool EqualsSS(Point3D* p) const
    {
        return p->X == X && p->Z == Z && p->Y == Y;
    }

//    int GetHashCode()
//    {
//        return (X + " " + Y + " " + Z).GetHashCode();
//    }

    string ToString() const
    {
        stringstream ss;

        ss << X << ", " << Y << ", " << Z;

        return ss.str();
    }

    static int refCounter;
};

}
