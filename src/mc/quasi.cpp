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
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <fstream>

#ifdef _WIN32
#include "drand48.h"
#endif

//#ifndef REGULAR
//#define REGULAR
//#endif

#ifndef RANDOM
#define RANDOM
#endif

#ifndef M_PI
#define M_PI (3.14159265358979323846f)
#endif

float evalFunc(const float &x, const float &y, const float &xmax, const float &ymax)
{
	return 1. / 2. + 1. / 2. * powf(1. - y / ymax, 3.) * sin( 2. * M_PI * ( x / xmax) * exp(8. * (x / xmax)));
}

int main(int argc, char **argv)
{
	uint32_t width = 512, height = 512;
	uint32_t nsamples = 1;
	unsigned char *pixels = new unsigned char[width * height];
	for (uint32_t y = 0; y < height; ++y) {
		for (uint32_t x = 0; x < width; ++x) {
			float sum = 0.f;
			for (uint32_t ny = 0; ny < nsamples; ++ny) {
				for (uint32_t nx = 0; nx < nsamples; ++nx) {
#ifdef REGULAR
					sum += evalFunc(x + (nx + 0.5) / nsamples, y + (ny + 0.5) / nsamples, width, height);
#endif
#ifdef RANDOM
					sum += evalFunc(x + drand48(), y + drand48(), width, height);
#endif
				}
			}
			pixels[y * width + x] = (unsigned char)(255 * sum / (nsamples * nsamples));
		}
	}

	std::ofstream ofs;
	ofs.open("./mc_quasi.ppm", std::ios::out | std::ios::binary);
	ofs << "P5\n" << width << " " << height << "\n255\n";
	ofs.write((char*)pixels, width * height);
	ofs.close();
	delete [] pixels;
	return 1;
}