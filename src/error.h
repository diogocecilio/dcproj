#include <iostream>
#include <stdlib.h>
#define DebugStop() DebugStopImpl(__FILE__, __LINE__)
void DebugStopImpl(const char *fileName, const std::size_t lineN);

