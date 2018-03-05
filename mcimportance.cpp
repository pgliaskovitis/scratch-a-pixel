#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "drand48.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846f)
#endif

int main(int argc, char **argv)
{
	srand48(13);
	int N = 16;
	for (int n = 0; n < 10; ++n) {
		float sumUniform = 0.f, sumImportance = 0.f;
		for (int i = 0; i < N; ++i) {
			float r = drand48();
			sumUniform += sin(r * M_PI * .5f);
			float xi = sqrtf(r) * M_PI * .5f; // this is our X_i
			sumImportance += sin(xi) / ((8 * xi) / (M_PI * M_PI));
		}
		sumUniform *= (M_PI * .5f) / N;
		sumImportance *= 1.f / N;
		printf("%f %f\n", sumUniform, sumImportance);
	}
	return 0;
} 