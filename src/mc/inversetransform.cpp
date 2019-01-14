/*
 * Copyright (C) 2016  www.scratchapixel.com
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
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <algorithm>

#ifdef _WIN32
#include "drand48.h"
#endif

#ifndef M_PI
#define M_PI (3.14159265358979323846f)
#endif

// standard normal distribution function
float pdf(const float &x)
{
	return 1 / sqrtf(2 * M_PI) * exp(-x * x * 0.5f);
}

float sample(float *cdf, const uint32_t &nbins, const float &minBound, const float &maxBound)
{
	float r = drand48();
	float *ptr = std::lower_bound(cdf, cdf + nbins + 1, r);
	int off = std::max(0, (int)(ptr - cdf - 1));
	float t = (r - cdf[off]) / (cdf[off + 1] - cdf[off]);
	float x = (off + t) / (float)(nbins);
	return minBound + (maxBound - minBound) * x;
}

int main(int argc, char ** argv)
{
	srand48(13); // seed random generator

	// create CDF
	int nbins = 32;
	float minBound = -5.0f, maxBound = 5.0f;
	float cdf[nbins + 1], dx = (maxBound - minBound) / nbins, sum = 0.0f;
	cdf[0] = 0.0f;

	for (int n = 1; n < nbins; ++n) {
		float x = minBound + (maxBound - minBound) * (n / (float)(nbins));
		float pdf_x = pdf(x) * dx;
		cdf[n] = cdf[n - 1] + pdf_x;
		sum += pdf_x;
	}

	cdf[nbins] = 1;

	// our simulation
	int numSims = 100000;
	int numBins = 100; // to collect data on our sim
	int bins[numBins]; // to collect data on our sim
	memset(bins, 0x0, sizeof(int) * numBins); // set all the bins to 0
	const float dist = 10.0f; // 10 km

	for (int i = 0; i < numSims; ++i) {
		float diff = sample(cdf, nbins, minBound, maxBound); // random var between -5 and 5
		int whichBin = (int)(numBins * (diff - minBound) / (maxBound - minBound));
		bins[whichBin]++;
		float time = 30.0f + diff;
		float speed = 60.0f * dist / time;
	}

	sum = 0;

	for (int i = 0; i < numBins; ++i) {
		float r = bins[i] / (float)numSims;
		printf("%f %f\n", 5 * (2 * (i /(float)(numBins)) -1), r);
		sum += r;
	}

	fprintf(stderr, "sum %f\n", sum);
	return 0;
}