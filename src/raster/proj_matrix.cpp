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

//[header]
// A simple program to demonstrate how to build and use a simple perspective projection matrix
//[/header]

#include <cstring>
#include <fstream>

#include "geometry.h"
#include "vertex_data.h"

//[comment]
// Set the basic perspective projection matrix
//[/comment]
void setProjectionMatrix(const float &angleOfView, const float &near, const float &far, Matrix44f &M)
{
	// do some work here
	float scale = 1 / tan(angleOfView * 0.5f * M_PI / 180);
	M[0][0] = scale;
	M[1][1] = scale;
	M[2][2] = -far / (far - near);
	M[3][2] = -far * near / (far - near);
	M[2][3] = -1;
	M[3][3] = 0;
}

//[comment]
// Point-Matrix multiplication. The input point is assumed to have Cartesian
// coordinates but because we will multiply this point by a 4x4 matrix we actually
// assume it is a point with homogeneous coordinatess (x, y, z, w = 1).
// The point-matrix results in another point with homogeneous coordinates (x', y', z', w').
// To get back to Cartesian coordinates we need to normalize these coordinates: (x'/w', y'/w', z'/w').
//[/comment]
void multPointMatrix(const Vec3f &in, Vec3f &out, const Matrix44f &M)
{
	//out = in * Mproj;
	out.x   = in.x * M[0][0] + in.y * M[1][0] + in.z * M[2][0] + /* in.z = 1 */ M[3][0];
	out.y   = in.x * M[0][1] + in.y * M[1][1] + in.z * M[2][1] + /* in.z = 1 */ M[3][1];
	out.z   = in.x * M[0][2] + in.y * M[1][2] + in.z * M[2][2] + /* in.z = 1 */ M[3][2];
	float w = in.x * M[0][3] + in.y * M[1][3] + in.z * M[2][3] + /* in.z = 1 */ M[3][3];

	// normalize if w is different than 1 (convert from homogeneous to Cartesian coordinates)
	if (w != 1) {
		out.x /= w;
		out.y /= w;
		out.z /= w;
	}
}

int main(int argc, char **argv)
{
	uint32_t imageWidth = 512, imageHeight = 512;
	Matrix44f Mproj;
	Matrix44f worldToCamera;
	worldToCamera[3][1] = -10.0f;
	worldToCamera[3][2] = -20.0f;
	float angleOfView = 90.0f;
	float near = 0.1f;
	float far = 100.0f;

	//[comment]
	// Set the basic perspective projection matrix
	//[/comment]
	setProjectionMatrix(angleOfView, near, far, Mproj);
	unsigned char *buffer = new unsigned char[imageWidth * imageHeight];
	memset(buffer, 0x0, imageWidth * imageHeight);

	//[comment]
	// Loop over all points
	//[/comment]
	for (uint32_t i = 0; i < numVertices; ++i) {
		Vec3f vertCamera, projectedVert;

		//[comment]
		// Transform to camera space
		//[/comment]
		multPointMatrix(vertices[i], vertCamera, worldToCamera);

		//[comment]
		// Project
		//[/comment]
		multPointMatrix(vertCamera, projectedVert, Mproj);
		if (projectedVert.x < -1 || projectedVert.x > 1 || projectedVert.y < -1 || projectedVert.y > 1) continue;
		// convert to raster space and mark the position of the vertex in the image with a simple dot
		uint32_t x = std::min(imageWidth - 1, (uint32_t)((projectedVert.x + 1) * 0.5 * imageWidth));
		uint32_t y = std::min(imageHeight - 1, (uint32_t)((1 - (projectedVert.y + 1) * 0.5f) * imageHeight));
		buffer[y * imageWidth + x] = 255;
	}
	// export to image
	std::ofstream ofs;
	ofs.open("./raster_projmatrix.ppm", std::ios::out | std::ios::binary);
	ofs << "P5\n" << imageWidth << " " << imageHeight << "\n255\n";
	ofs.write((char*)buffer, imageWidth * imageHeight);
	ofs.close();
	delete [] buffer;

	return 0;
}
