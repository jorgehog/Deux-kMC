#pragma once

#include <vector>
#include <sys/types.h>
#include <algorithm>

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

    Subject() {}

    virtual ~Subject() {}

    void registerObserver(ObserverType *observer)
    {
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
    }

    void notifyObservers(const IDType &ID)
    {
        for (ObserverType *observer : m_observers)
        {
            observer->notifyObserver(ID);
        }
    }

private:
    std::vector<ObserverType* > m_observers;

};

}
