#include "error.h"

void DebugStopImpl(const char *fileName, const std::size_t lineN)
{
    std::cerr << "Your chance to put a breakpoint at " << fileName<< ":"<< lineN <<  "\n";
    throw std::bad_exception();
}
