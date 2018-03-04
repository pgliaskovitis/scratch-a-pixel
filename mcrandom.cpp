#include <cstdio>
#include <cstdlib>
#include <cmath>

int main(int argc, char **argv)
{
	float seed = 0.4872f, rand = seed;
	int seqlength = 10;
	while (seqlength--) {
		printf("%0.4f\n", rand);
		rand = sqrtf(rand);
		rand = ((int)((rand * 100 - (int)(rand * 100)) * 10000)) / 10000.f;
	}
	return 0;
} 