#pragma once

#include <vector>
#include <sys/types.h>
#include <algorithm>

#include <iostream>

namespace kMC
{

template<class IDType = int>
class Observer
{
public:
    Observer() {}
    virtual ~Observer() {}

    virtual void initializeObserver(const IDType &subject) = 0;

    virtual void notifyObserver(const IDType &subject) = 0;
};

template<class IDType = int>
class Subject
{
public:
    using ObserverType = Observer<IDType>;

    Subject() : m_initialized(false) {}

    virtual ~Subject() {}

    void registerObserver(ObserverType *observer)
    {
        if (std::find(m_observers.begin(), m_observers.end(), observer) != m_observers.end())
        {
            std::cerr << "duplicate observer added." << std::endl;
        }

        m_observers.push_back(observer);
    }

    void removeObserver(ObserverType *observer)
    {
        auto &v = m_observers;
        v.erase(std::remove(v.begin(), v.end(), observer), v.end());
    }

    void initializeObservers(const IDType &ID)
    {
        for (ObserverType *observer : m_observers)
        {
            observer->initializeObserver(ID);
        }

        m_initialized = true;
    }

    void notifyObservers(const IDType &ID)
    {
        if (!m_initialized)
        {
            return;
        }

        for (ObserverType *observer : m_observers)
        {
            observer->notifyObserver(ID);
        }
    }

private:
    std::vector<ObserverType* > m_observers;
    bool m_initialized;

};

}
