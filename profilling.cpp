#include <iostream>

#define PROFILLING 1
#include "profiling.h"

int main()
{
	{
		ProfillingTimer t("name");
		for(int i=0;i<1e9;i++);
	}
	
	std::cout << Profiler::Get().LastTime() << "\n";

	{
		ProfillingTimer t("name");
		for(int i=0;i<1e9;i++);
	}
	
	std::cout << Profiler::Get().LastTime() << "\n";
	
	return 0;
}

void profiling(){main();}