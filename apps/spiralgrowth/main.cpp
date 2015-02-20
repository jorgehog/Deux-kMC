#include <iostream>

#include <armadillo>
#include <utils.h>

using namespace arma;
using namespace std;

class test : public ignis::LatticeEvent
{
public:

    test() : ignis::LatticeEvent("test", "!", true) {}

    void execute()
    {
        setValue(cycle());
    }
};

int main()
{
    ignis::Lattice lattice;

    lattice.addEvent(new test);

    lattice.eventLoop(1000);

    return 0;
}

