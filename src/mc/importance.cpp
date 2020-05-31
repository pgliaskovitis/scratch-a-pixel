/*
 * Copyright (C) 2012  www.scratchapixel.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>

#ifdef _WIN32
#include "drand48.h"
#endif

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