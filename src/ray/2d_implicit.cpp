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

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <cmath>

float evalImplicitFunc(const float &x, const float &y)
{
	// to draw an implicit square
	return std::max(std::fabs(x), std::fabs(y)) - 1;
	// to draw an implicit sphere
	// return x * x + y * y - 1;
}

int main(int argc, char **argv)
{
	unsigned int height = 512, width = 512;
	unsigned char *image = new unsigned char[width * height];
	const float scale = 2.f;

	for (unsigned int j = 0; j < width; ++j) {
		for (unsigned i = 0; i < height; ++i) {
			float x = (2 * i / (float)width - 1) * scale;
			float y = (2 * j / (float)height - 1) * scale;
			float d = std::min(1.f, powf(std::fabs(evalImplicitFunc(x, y)), 0.3f));
			image[j * width + i] = static_cast<unsigned char>(d * 255);
		}
	}

	std::ofstream ofs;
	ofs.open("./ray_2d_implicit.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << width << " " << height << "\n255\n";
	for (unsigned int i = 0; i < width * height; ++i)
	ofs << image[i] << image[i] << image[i];
	ofs.close();

	delete [] image;

	return 0;
}