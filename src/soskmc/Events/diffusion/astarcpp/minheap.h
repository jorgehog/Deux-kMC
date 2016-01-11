#pragma once

#include "searchnode.h"

#include <iostream>
using std::cout;
using std::endl;

///This code is adapted from a C# implementation of Roy Triesscheijn.
/// Original comments are included below.

namespace Tests
{
/// <summary>
/// MinHeap from ZeraldotNet (http://zeraldotnet.codeplex.com/)
/// Modified by Roy Triesscheijn (http://roy-t.nl)
/// -Moved method variables to class variables
/// -Added English Exceptions and comments (instead of Chinese)
/// </summary>
class MinHeap
{
public:
    MinHeap() :
        listHead(nullptr)
    {

    }

    bool HasNext()
    {
        return listHead != nullptr;
    }

    void Add(SearchNode* item)
    {
        if (listHead == nullptr)
        {
            listHead = item;
        }

        else if (listHead->next == nullptr && item->cost <= listHead->cost)
        {
            item->nextListElem = listHead;
            listHead = item;
        }
        else
        {
            SearchNode* ptr = listHead;
            while (ptr->nextListElem != nullptr && ptr->nextListElem->cost < item->cost)
                ptr = ptr->nextListElem;
            item->nextListElem = ptr->nextListElem;
            ptr->nextListElem = item;
        }
    }

    SearchNode* extractFirst()
    {
        SearchNode* result = listHead;
        listHead = listHead->nextListElem;

        return result;
    }

private:
    SearchNode* listHead;
};

}
