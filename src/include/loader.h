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

#include "geometry.h"

namespace
{
	void loadGeoFile(
		const char *file,
		uint32_t &numFaces,
		std::unique_ptr<Vec3f []> &verts,
		std::unique_ptr<Vec2f []> &st,
		std::unique_ptr<uint32_t []> &vertsIndex)
	{
		std::ifstream ifs;
		try {
			ifs.open(file);
			if (ifs.fail()) throw;
			std::stringstream ss;
			ss << ifs.rdbuf();
			ss >> numFaces;

			uint32_t vertsIndexArraySize = 0;
			// reading face index array
			for (uint32_t i = 0; i < numFaces; ++i) {
				uint32_t tmp;
				ss >> tmp; //faceIndex[i];
				vertsIndexArraySize += tmp; //faceIndex[i];
			}
			vertsIndex = std::unique_ptr<uint32_t []>(new uint32_t[vertsIndexArraySize]);

			uint32_t vertsArraySize = 0;
			// reading vertex index array
			for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
				ss >> vertsIndex[i];
				if (vertsIndex[i] > vertsArraySize) vertsArraySize = vertsIndex[i];
			}
			vertsArraySize += 1;

			// reading vertices
			verts = std::unique_ptr<Vec3f []>(new Vec3f[vertsArraySize]);
			for (uint32_t i = 0; i < vertsArraySize; ++i) {
				ss >> verts[i].x >> verts[i].y >> verts[i].z;
			}

			// reading normals
			for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
				Vec3f normal;
				ss >> normal.x >> normal.y >> normal.z;
			}

			// reading st coordinates
			st = std::unique_ptr<Vec2f []>(new Vec2f[vertsIndexArraySize]);
			for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
				ss >> st[i].x >> st[i].y;
			}
		}
		catch (...) {
			ifs.close();
		}
		ifs.close();
	}
}
