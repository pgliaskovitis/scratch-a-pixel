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
// Digital Images
//[/header]

#include "images.h"

// [comment]
// Read/Write an image stored in the PPM format
// [/comment]
int main(int argc, char **argv)
{
	Image I = readPPM("data/xmas.ppm");
	savePPM(I, "./img_readwrite.ppm");

	return 0;
}