#pragma once

#include <BADAss/badass.h>

#include <armadillo>

using std::vector;
using arma::vec;

namespace kMC
{

inline uint binarySearchForInterval(const double target, const double* intervals, uint size)
{

    BADAss(size, !=, 0u, "Number of intervals cannot be zero.");

    uint imax = size - 1;
    uint MAX = imax;

    uint imin = 0;
    uint imid;

    // continue searching while imax != imin + 1
    do
    {

        // calculate the midpoint for roughly equal partition
        imid = imin + (imax - imin)/2;

        //Is the upper limit above mid?
        if (target > intervals[imid])
        {

            //This means that the target is the last interval.
            if (imid == MAX)
            {
                return MAX;
            }

            //Are we just infront of the limit?
            else if (target < intervals[imid + 1])
            {
                //yes we were! Returning current mid + 1.
                //If item i in accuAllrates > R, then reaction i is selected.
                //This because there is no zero at the start of accuAllrates.

                return imid + 1;
            }

            //No we're not there yet, so we search above us.
            else
            {
                imin = imid + 1;
            }
        }

        //No it's not. Starting new search below mid!
        else
        {

            //This means that the target is the first inteval.
            if (imid == 0)
            {
                return 0;
            }

            imax = imid;
        }


    } while (imid != imin);

    //If we get here, imid = imin, which means that imax = imid + 1 (deduced by integer division).
    //We choose the max value as out match.
    return imid + 1;

}

inline uint binarySearchForInterval(const double target, const vector<double> &intervals)
{
    return binarySearchForInterval(target, &intervals.front(), intervals.size());
}

inline uint binarySearchForInterval(const double target, const vec &intervals)
{
    return binarySearchForInterval(target, intervals.memptr(), intervals.size());
}


} //end of namespace kMC
