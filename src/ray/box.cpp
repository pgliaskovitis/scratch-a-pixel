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
// Ray trace a bounding box
//[/header]

#include <cstdlib>
#include <random>

#include "objects.h"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

int main(int argc, char **argv)
{
	BBox<float> box(Vec3f(-1), Vec3f(1));
	gen.seed(0);
	for (uint32_t i = 0; i < 16; ++i) {
		Vec3f randDir(2 * dis(gen) - 1, 2 * dis(gen) - 1, 2 * dis(gen) - 1);
		randDir.normalize();
		Ray ray(Vec3f(0.f), randDir);
		float t;
		if (box.intersect(ray, t)) {
			Vec3f Phit = ray.orig + ray.dir * t;
			std::cerr << ray.orig << " " << Phit << std::endl;
		}
	}
	return 0;
}
