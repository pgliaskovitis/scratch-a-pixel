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

#include <cstring>
#include <fstream>

#include "utils.h"

typedef float Vec2[2];
typedef float Vec3[3];
typedef unsigned char Rgb[3];

int main(int argc, char **argv)
{
	// vertices
	Vec2 v0 = {491.407f, 411.407f};
	Vec2 v1 = {148.593f, 68.5928f};
	Vec2 v2 = {148.593f, 411.407f};

	// colors
	Vec3 c0 = {1.0f, 0.0f, 0.0f};
	Vec3 c1 = {0.0f, 1.0f, 0.0f};
	Vec3 c2 = {0.0f, 0.0f, 1.0f};

	const uint32_t width = 1920;
	const uint32_t height = 1080;

	Rgb *framebuffer = new Rgb[width * height];
	memset(framebuffer, 0x0, width * height * 3);

	float area = scratch::utils::edgeFunction(v0, v1, v2);

	for (uint32_t j = 0; j < height; ++j) {
		for (uint32_t i = 0; i < width; ++i) {
			Vec2 p = {i + 0.5f, j + 0.5f};
			float w0 = scratch::utils::edgeFunction(v1, v2, p);
			float w1 = scratch::utils::edgeFunction(v2, v0, p);
			float w2 = scratch::utils::edgeFunction(v0, v1, p);
			if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
				w0 /= area;
				w1 /= area;
				w2 /= area;
				float r = w0 * c0[0] + w1 * c1[0] + w2 * c2[0];
				float g = w0 * c0[1] + w1 * c1[1] + w2 * c2[1];
				float b = w0 * c0[2] + w1 * c1[2] + w2 * c2[2];
				framebuffer[j * width + i][0] = (unsigned char)(r * 255.0f);
				framebuffer[j * width + i][1] = (unsigned char)(g * 255.0f);
				framebuffer[j * width + i][2] = (unsigned char)(b * 255.0f);
			}
		}
	}

	std::ofstream ofs;
	ofs.open("./raster2d_triangle.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << width << " " << height << "\n255\n";
	ofs.write((char*)framebuffer, width * height * 3);
	ofs.close();

	delete [] framebuffer;

	return 0;
}