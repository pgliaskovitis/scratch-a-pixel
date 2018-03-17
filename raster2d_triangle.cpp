// c++ -o raster2d raster2d.cpp
// (c) www.scratchapixel.com

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <algorithm>

typedef float Vec2[2];
typedef float Vec3[3];
typedef unsigned char Rgb[3];

inline
float edgeFunction(const Vec2 &a, const Vec2 &b, const Vec2 &c)
{
	return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]);
}

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

	const uint32_t width = 640;
	const uint32_t height = 480;

	Rgb *framebuffer = new Rgb[width * height];
	memset(framebuffer, 0x0, width * height * 3);

	float area = edgeFunction(v0, v1, v2);

	for (uint32_t j = 0; j < height; ++j) {
		for (uint32_t i = 0; i < width; ++i) {
			Vec2 p = {i + 0.5f, j + 0.5f};
			float w0 = edgeFunction(v1, v2, p);
			float w1 = edgeFunction(v2, v0, p);
			float w2 = edgeFunction(v0, v1, p);
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