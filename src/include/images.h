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

#include <cstring>
#include <cassert>
#include <fstream>
#include <algorithm>

// [comment]
// The main Image class
// [/comment]
class Image
{
public:

	struct Rgb
	{
		Rgb() : r(0), g(0), b(0)  {}
		Rgb(float c) : r(c), g(c), b(c) {}
		Rgb(float _r, float _g, float _b) : r(_r), g(_g), b(_b) {}
		bool operator != (const Rgb &c) const { return c.r != r && c.g != g && c.b != b; }
		Rgb& operator *= (const Rgb &rgb) { r *= rgb.r, g *= rgb.g, b *= rgb.b; return *this; }
		Rgb& operator += (const Rgb &rgb) { r += rgb.r, g += rgb.g, b += rgb.b; return *this; }
		friend float& operator += (float &f, const Rgb rgb)
		{ f += (rgb.r + rgb.g + rgb.b) / 3.f; return f; }
		float r, g, b;
	};

	Image() : w(0), h(0), pixels(nullptr)
	{ /* empty image */ }

	Image(const unsigned int &_w, const unsigned int &_h, const Rgb &c = kBlack) :
		w(_w), h(_h), pixels(nullptr)
	{
		pixels = new Rgb[w * h];
		for (unsigned int i = 0; i < w * h; ++i) pixels[i] = c;
	}
	Image(const Image &img) : w(img.w), h(img.h), pixels(nullptr)
	{
		pixels = new Rgb[w * h];
		memcpy(pixels, img.pixels, sizeof(Rgb) * w * h);
	}
	// move constructor
	Image(Image &&img) : w(0), h(0), pixels(nullptr)
	{
		w = img.w;
		h = img.h;
		pixels = img.pixels;
		img.pixels = nullptr;
		img.w = img.h = 0;
	}
	// move assignment operator
	Image& operator = (Image &&img)
	{
		if (this != &img) {
			if (pixels != nullptr) delete [] pixels;
			w = img.w, h = img.h;
			pixels = img.pixels;
			img.pixels = nullptr;
			img.w = img.h = 0;
		}
		return *this;
	}
	Rgb& operator () (const unsigned &x, const unsigned int &y) const
	{
		assert(x < w && y < h);
		return pixels[y * w + x];
	}
	Image& operator *= (const Rgb &rgb)
	{
		for (unsigned int i = 0; i < w * h; ++i) pixels[i] *= rgb;
		return *this;
	}
	Image& operator += (const Image &img)
	{
		for (unsigned int i = 0; i < w * h; ++i) pixels[i] += img[i];
		return *this;
	}
	Image& operator /= (const float &div)
	{
		float invDiv = 1 / div;
		for (unsigned int i = 0; i < w * h; ++i) pixels[i] *= invDiv;
		return *this;
	}
	friend Image operator * (const Rgb &rgb, const Image &img)
	{
		Image tmp(img);
		tmp *= rgb;
		return tmp;
	}
	Image operator * (const Image &img)
	{
		Image tmp(*this);
		// multiply pixels together
		for (unsigned int i = 0; i < w * h; ++i) tmp[i] *= img[i];
		return tmp;
	}
	static Image circshift(const Image &img, const std::pair<int,int> &shift)
	{
		Image tmp(img.w, img.h);
		unsigned int w = img.w, h = img.h;
		for (unsigned int j = 0; j < h; ++j) {
			int jmod = (j + shift.second) % h;
			for (unsigned int i = 0; i < w; ++i) {
				int imod = (i + shift.first) % w;
				tmp[jmod * w + imod] = img[j * w + i];
			}
		}
		return tmp;
	}
	const Rgb& operator [] (uint8_t i) const { return pixels[i]; }
	Rgb& operator [] (uint8_t i) { return pixels[i]; }
	~Image() { if (pixels != nullptr) delete [] pixels; }
	unsigned int w, h;
	Rgb *pixels;
	static const Rgb kBlack, kWhite, kRed, kGreen, kBlue;
};

// [comment]
// Save an image to PPM image file
// [/comment]
void savePPM(const Image &img, const char *filename)
{
	if (img.w == 0 || img.h == 0) { fprintf(stderr, "Can't save an empty image\n"); return; }
	std::ofstream ofs;
	try {
		ofs.open(filename, std::ios::binary); // need to spec. binary mode for Windows users
		if (ofs.fail()) throw("Can't open output file");
		ofs << "P6\n" << img.w << " " << img.h << "\n255\n";
		unsigned char r, g, b;
		// loop over each pixel in the image, clamp and convert to byte format
		for (unsigned int i = 0; i < img.w * img.h; ++i) {
			r = static_cast<unsigned char>(std::min(1.f, img.pixels[i].r) * 255);
			g = static_cast<unsigned char>(std::min(1.f, img.pixels[i].g) * 255);
			b = static_cast<unsigned char>(std::min(1.f, img.pixels[i].b) * 255);
			ofs << r << g << b;
		}
		ofs.close();
	}
	catch (const char *err) {
		fprintf(stderr, "%s\n", err);
		ofs.close();
	}
}

// [comment]
// Read a PPM image file
// [/comment]
Image readPPM(const char *filename)
{
	std::ifstream ifs;
	ifs.open(filename, std::ios::binary); // need to spec. binary mode for Windows users
	Image img;
	try {
		if (ifs.fail()) { throw("Can't open input file"); }
		std::string header;
		unsigned int w, h, b;
		ifs >> header;
		if (strcmp(header.c_str(), "P6") != 0) throw("Can't read input file");
		ifs >> w >> h >> b;
		img.w = w; img.h = h;
		img.pixels = new Image::Rgb[w * h]; // this is throw an exception if bad_alloc
		ifs.ignore(256, '\n'); // skip empty lines in necessary until we get to the binary data
		unsigned char pix[3];
		// read each pixel one by one and convert bytes to floats
		for (unsigned int i = 0; i < w * h; ++i) {
			ifs.read(reinterpret_cast<char *>(pix), 3);
			img.pixels[i].r = pix[0] / 255.f;
			img.pixels[i].g = pix[1] / 255.f;
			img.pixels[i].b = pix[2] / 255.f;

			// [comment]
			// This is just to make the bokeh effect more visible, book high value pixels brightness
			// [/comment]
			// if (img.pixels[i].r > 0.7) img.pixels[i].r *= 3;
			// if (img.pixels[i].g > 0.7) img.pixels[i].g *= 3;
			// if (img.pixels[i].b > 0.7) img.pixels[i].b *= 3;
		}
		ifs.close();
	}
	catch (const char *err) {
		fprintf(stderr, "%s\n", err);
		ifs.close();
	}

	return img;
}

const Image::Rgb Image::kBlack = Image::Rgb(0);
const Image::Rgb Image::kWhite = Image::Rgb(1);
const Image::Rgb Image::kRed = Image::Rgb(1,0,0);
const Image::Rgb Image::kGreen = Image::Rgb(0,1,0);
const Image::Rgb Image::kBlue = Image::Rgb(0,0,1);