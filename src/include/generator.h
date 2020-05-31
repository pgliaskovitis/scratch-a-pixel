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

#pragma once

#include <sstream>

#include "objects.h"

namespace scratch
{
namespace generator
{
	TriangleMesh* generatePolySphere(float rad, uint32_t divs)
	{
		// generate points
		uint32_t numVertices = (divs - 1) * divs + 2;
		std::unique_ptr<Vec3f []> P(new Vec3f[numVertices]);
		std::unique_ptr<Vec3f []> N(new Vec3f[numVertices]);
		std::unique_ptr<Vec2f []> st(new Vec2f[numVertices]);

		float u = -M_PI_2;
		float v = -M_PI;
		float du = M_PI / divs;
		float dv = 2 * M_PI / divs;

		P[0] = N[0] = Vec3f(0, -rad, 0);
		uint32_t k = 1;
		for (uint32_t i = 0; i < divs - 1; i++) {
			u += du;
			v = -M_PI;
			for (uint32_t j = 0; j < divs; j++) {
				float x = rad * cos(u) * cos(v);
				float y = rad * sin(u);
				float z = rad * cos(u) * sin(v) ;
				P[k] = N[k] = Vec3f(x, y, z);
				st[k].x = u / M_PI + 0.5;
				st[k].y = v * 0.5 / M_PI + 0.5;
				v += dv, k++;
			}
		}
		P[k] = N[k] = Vec3f(0, rad, 0);

		uint32_t npolys = divs * divs;
		std::unique_ptr<uint32_t []> faceIndex(new uint32_t[npolys]);
		std::unique_ptr<uint32_t []> vertsIndex(new uint32_t[(6 + (divs - 1) * 4) * divs]);

		// create the connectivity lists
		uint32_t vid = 1, numV = 0, l = 0;
		k = 0;
		for (uint32_t i = 0; i < divs; i++) {
			for (uint32_t j = 0; j < divs; j++) {
				if (i == 0) {
					faceIndex[k++] = 3;
					vertsIndex[l] = 0;
					vertsIndex[l + 1] = j + vid;
					vertsIndex[l + 2] = (j == (divs - 1)) ? vid : j + vid + 1;
					l += 3;
				}
				else if (i == (divs - 1)) {
					faceIndex[k++] = 3;
					vertsIndex[l] = j + vid + 1 - divs;
					vertsIndex[l + 1] = vid + 1;
					vertsIndex[l + 2] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs;
					l += 3;
				}
				else {
					faceIndex[k++] = 4;
					vertsIndex[l] = j + vid + 1 - divs;
					vertsIndex[l + 1] = j + vid + 1;
					vertsIndex[l + 2] = (j == (divs - 1)) ? vid + 1 : j + vid + 2;
					vertsIndex[l + 3] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs;
					l += 4;
				}
				numV++;
			}
			vid = numV;
		}

		return new TriangleMesh(npolys, faceIndex, vertsIndex, P, N, st);
	}
}
}
