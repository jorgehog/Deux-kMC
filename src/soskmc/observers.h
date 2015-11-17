#pragma once

#include <vector>
#include <sys/types.h>

class Observer
{
public:
    Observer();
    virtual ~Observer();

    virtual void initializeObserver() = 0;

    virtual void notifyObserver() = 0;
};

//!Makes a connection between a height change and the inherited class
class HeightObserver: public Observer {};

//!Makes a connection between a change in the confining surface and the inherited classc
class ConfiningSurfaceObserver : public Observer {};
