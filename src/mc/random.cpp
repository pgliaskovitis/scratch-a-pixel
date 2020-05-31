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
#include <cstdio>
#include <cmath>

int main(int argc, char **argv)
{
	float seed = 0.4872f, rand = seed;
	int seqlength = 10;
	while (seqlength--) {
		printf("%0.4f\n", rand);
		rand = sqrtf(rand);
		rand = ((int)((rand * 100 - (int)(rand * 100)) * 10000)) / 10000.f;
	}
	return 0;
}