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
// Simple image manipulations
//[/header]

#include "images.h"

// [comment]
// Simulate the bokeh effect (downalod the test images or create your own)
// [/comment]
int main(int argc, char **argv)
{
	try {
		Image I = readPPM("data/xmas.ppm");
		Image J = readPPM("data/heart.ppm");
		int w = J.w, h = J.h;
		Image K(w, h);
		float total = 0;
		for (int j = 0; j < h; ++j) {
			for (int i = 0; i < w; ++i) {
				if (J(i, j) != Image::kBlack) {
					K += J(i, j) * Image::circshift(I, std::pair<int, int>(i, j));
					total += J(i, j);
				}
			}
		}
		K /= total;
		savePPM(K, "./img_bokeh.ppm");
	}
	catch (const std::exception &e) { // catch general exception (bad_alloc mainly?)
		fprintf(stderr, "Error: %s\n", e.what());
	}

	return 0;
}