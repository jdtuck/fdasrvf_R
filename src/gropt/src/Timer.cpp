
#include "Timer.h"

unsigned long getTickCount(void)
{
//	return clock();

	unsigned long currentTime = 0;
    
    
#ifdef _WIN64
	currentTime = GetTickCount();
#elif _WIN32
	currentTime = GetTickCount();
#elif __APPLE__
    #include "TargetConditionals.h"
    #if TARGET_OS_MAC
	struct timeval current;
	gettimeofday(&current, NULL);
	currentTime = current.tv_sec * 1000 + current.tv_usec / 1000;
    #endif
#elif __linux
struct timeval current;
gettimeofday(&current, NULL);
currentTime = current.tv_sec * 1000 + current.tv_usec / 1000;
#endif // end of checking platforms

	return currentTime;
}
