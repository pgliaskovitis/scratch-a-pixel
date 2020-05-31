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
// Simple example that demonstrates how to convert a color defined as a spectrum to RGB
//[/header]

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>

// [comment]
// Color Matching Function: 380-730x5
// [/comment]
const float colorMatchingFunc[3][72] =
{
	{0.001368f, 0.002236f, 0.004243f, 0.007650f, 0.014310f, 0.023190f, 0.043510f, 0.077630f, 0.134380f, 0.214770f, 0.283900f, 0.328500f, 0.348280f, 0.348060f, 0.336200f, 0.318700f, 0.290800f, 0.251100f, 0.195360f, 0.142100f, 0.095640f, 0.057950f, 0.032010f, 0.014700f, 0.004900f, 0.002400f, 0.009300f, 0.029100f, 0.063270f, 0.109600f, 0.165500f, 0.225750f, 0.290400f, 0.359700f, 0.433450f, 0.512050f, 0.594500f, 0.678400f, 0.762100f, 0.842500f, 0.916300f, 0.978600f, 1.026300f, 1.056700f, 1.062200f, 1.045600f, 1.002600f, 0.938400f, 0.854450f, 0.751400f, 0.642400f, 0.541900f, 0.447900f, 0.360800f, 0.283500f, 0.218700f, 0.164900f, 0.121200f, 0.087400f, 0.063600f, 0.046770f, 0.032900f, 0.022700f, 0.015840f, 0.011359f, 0.008111f, 0.005790f, 0.004106f, 0.002899f, 0.002049f, 0.001440f, 0.000000f},
	{0.000039f, 0.000064f, 0.000120f, 0.000217f, 0.000396f, 0.000640f, 0.001210f, 0.002180f, 0.004000f, 0.007300f, 0.011600f, 0.016840f, 0.023000f, 0.029800f, 0.038000f, 0.048000f, 0.060000f, 0.073900f, 0.090980f, 0.112600f, 0.139020f, 0.169300f, 0.208020f, 0.258600f, 0.323000f, 0.407300f, 0.503000f, 0.608200f, 0.710000f, 0.793200f, 0.862000f, 0.914850f, 0.954000f, 0.980300f, 0.994950f, 1.000000f, 0.995000f, 0.978600f, 0.952000f, 0.915400f, 0.870000f, 0.816300f, 0.757000f, 0.694900f, 0.631000f, 0.566800f, 0.503000f, 0.441200f, 0.381000f, 0.321000f, 0.265000f, 0.217000f, 0.175000f, 0.138200f, 0.107000f, 0.081600f, 0.061000f, 0.044580f, 0.032000f, 0.023200f, 0.017000f, 0.011920f, 0.008210f, 0.005723f, 0.004102f, 0.002929f, 0.002091f, 0.001484f, 0.001047f, 0.000740f, 0.000520f, 0.000000f},
	{0.006450f, 0.010550f, 0.020050f, 0.036210f, 0.067850f, 0.110200f, 0.207400f, 0.371300f, 0.645600f, 1.039050f, 1.385600f, 1.622960f, 1.747060f, 1.782600f, 1.772110f, 1.744100f, 1.669200f, 1.528100f, 1.287640f, 1.041900f, 0.812950f, 0.616200f, 0.465180f, 0.353300f, 0.272000f, 0.212300f, 0.158200f, 0.111700f, 0.078250f, 0.057250f, 0.042160f, 0.029840f, 0.020300f, 0.013400f, 0.008750f, 0.005750f, 0.003900f, 0.002750f, 0.002100f, 0.001800f, 0.001650f, 0.001400f, 0.001100f, 0.001000f, 0.000800f, 0.000600f, 0.000340f, 0.000240f, 0.000190f, 0.000100f, 0.000050f, 0.000030f, 0.000020f, 0.000010f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f}
};

// [comment]
// McBeth Spectral data: 380-730x10
// [/comment]
const float spectralData[24][36] = {
	{0.055f, 0.058f, 0.061f, 0.062f, 0.062f, 0.062f, 0.062f, 0.062f, 0.062f, 0.062f, 0.062f, 0.063f, 0.065f, 0.070f, 0.076f, 0.079f, 0.081f, 0.084f, 0.091f, 0.103f, 0.119f, 0.134f, 0.143f, 0.147f, 0.151f, 0.158f, 0.168f, 0.179f, 0.188f, 0.190f, 0.186f, 0.181f, 0.182f, 0.187f, 0.196f, 0.209f},
	{0.117f, 0.143f, 0.175f, 0.191f, 0.196f, 0.199f, 0.204f, 0.213f, 0.228f, 0.251f, 0.280f, 0.309f, 0.329f, 0.333f, 0.315f, 0.286f, 0.273f, 0.276f, 0.277f, 0.289f, 0.339f, 0.420f, 0.488f, 0.525f, 0.546f, 0.562f, 0.578f, 0.595f, 0.612f, 0.625f, 0.638f, 0.656f, 0.678f, 0.700f, 0.717f, 0.734f},
	{0.130f, 0.177f, 0.251f, 0.306f, 0.324f, 0.330f, 0.333f, 0.331f, 0.323f, 0.311f, 0.298f, 0.285f, 0.269f, 0.250f, 0.231f, 0.214f, 0.199f, 0.185f, 0.169f, 0.157f, 0.149f, 0.145f, 0.142f, 0.141f, 0.141f, 0.141f, 0.143f, 0.147f, 0.152f, 0.154f, 0.150f, 0.144f, 0.136f, 0.132f, 0.135f, 0.147f},
	{0.051f, 0.054f, 0.056f, 0.057f, 0.058f, 0.059f, 0.060f, 0.061f, 0.062f, 0.063f, 0.065f, 0.067f, 0.075f, 0.101f, 0.145f, 0.178f, 0.184f, 0.170f, 0.149f, 0.133f, 0.122f, 0.115f, 0.109f, 0.105f, 0.104f, 0.106f, 0.109f, 0.112f, 0.114f, 0.114f, 0.112f, 0.112f, 0.115f, 0.120f, 0.125f, 0.130f},
	{0.144f, 0.198f, 0.294f, 0.375f, 0.408f, 0.421f, 0.426f, 0.426f, 0.419f, 0.403f, 0.379f, 0.346f, 0.311f, 0.281f, 0.254f, 0.229f, 0.214f, 0.208f, 0.202f, 0.194f, 0.193f, 0.200f, 0.214f, 0.230f, 0.241f, 0.254f, 0.279f, 0.313f, 0.348f, 0.366f, 0.366f, 0.359f, 0.358f, 0.365f, 0.377f, 0.398f},
	{0.136f, 0.179f, 0.247f, 0.297f, 0.320f, 0.337f, 0.355f, 0.381f, 0.419f, 0.466f, 0.510f, 0.546f, 0.567f, 0.574f, 0.569f, 0.551f, 0.524f, 0.488f, 0.445f, 0.400f, 0.350f, 0.299f, 0.252f, 0.221f, 0.204f, 0.196f, 0.191f, 0.188f, 0.191f, 0.199f, 0.212f, 0.223f, 0.232f, 0.233f, 0.229f, 0.229f},
	{0.054f, 0.054f, 0.053f, 0.054f, 0.054f, 0.055f, 0.055f, 0.055f, 0.056f, 0.057f, 0.058f, 0.061f, 0.068f, 0.089f, 0.125f, 0.154f, 0.174f, 0.199f, 0.248f, 0.335f, 0.444f, 0.538f, 0.587f, 0.595f, 0.591f, 0.587f, 0.584f, 0.584f, 0.590f, 0.603f, 0.620f, 0.639f, 0.655f, 0.663f, 0.663f, 0.667f},
	{0.122f, 0.164f, 0.229f, 0.286f, 0.327f, 0.361f, 0.388f, 0.400f, 0.392f, 0.362f, 0.316f, 0.260f, 0.209f, 0.168f, 0.138f, 0.117f, 0.104f, 0.096f, 0.090f, 0.086f, 0.084f, 0.084f, 0.084f, 0.084f, 0.084f, 0.085f, 0.090f, 0.098f, 0.109f, 0.123f, 0.143f, 0.169f, 0.205f, 0.244f, 0.287f, 0.332f},
	{0.096f, 0.115f, 0.131f, 0.135f, 0.133f, 0.132f, 0.130f, 0.128f, 0.125f, 0.120f, 0.115f, 0.110f, 0.105f, 0.100f, 0.095f, 0.093f, 0.092f, 0.093f, 0.096f, 0.108f, 0.156f, 0.265f, 0.399f, 0.500f, 0.556f, 0.579f, 0.588f, 0.591f, 0.593f, 0.594f, 0.598f, 0.602f, 0.607f, 0.609f, 0.609f, 0.610f},
	{0.092f, 0.116f, 0.146f, 0.169f, 0.178f, 0.173f, 0.158f, 0.139f, 0.119f, 0.101f, 0.087f, 0.075f, 0.066f, 0.060f, 0.056f, 0.053f, 0.051f, 0.051f, 0.052f, 0.052f, 0.051f, 0.052f, 0.058f, 0.073f, 0.096f, 0.119f, 0.141f, 0.166f, 0.194f, 0.227f, 0.265f, 0.309f, 0.355f, 0.396f, 0.436f, 0.478f},
	{0.061f, 0.061f, 0.062f, 0.063f, 0.064f, 0.066f, 0.069f, 0.075f, 0.085f, 0.105f, 0.139f, 0.192f, 0.271f, 0.376f, 0.476f, 0.531f, 0.549f, 0.546f, 0.528f, 0.504f, 0.471f, 0.428f, 0.381f, 0.347f, 0.327f, 0.318f, 0.312f, 0.310f, 0.314f, 0.327f, 0.345f, 0.363f, 0.376f, 0.381f, 0.378f, 0.379f},
	{0.063f, 0.063f, 0.063f, 0.064f, 0.064f, 0.064f, 0.065f, 0.066f, 0.067f, 0.068f, 0.071f, 0.076f, 0.087f, 0.125f, 0.206f, 0.305f, 0.383f, 0.431f, 0.469f, 0.518f, 0.568f, 0.607f, 0.628f, 0.637f, 0.640f, 0.642f, 0.645f, 0.648f, 0.651f, 0.653f, 0.657f, 0.664f, 0.673f, 0.680f, 0.684f, 0.688f},
	{0.066f, 0.079f, 0.102f, 0.146f, 0.200f, 0.244f, 0.282f, 0.309f, 0.308f, 0.278f, 0.231f, 0.178f, 0.130f, 0.094f, 0.070f, 0.054f, 0.046f, 0.042f, 0.039f, 0.038f, 0.038f, 0.038f, 0.038f, 0.039f, 0.039f, 0.040f, 0.041f, 0.042f, 0.044f, 0.045f, 0.046f, 0.046f, 0.048f, 0.052f, 0.057f, 0.065f},
	{0.052f, 0.053f, 0.054f, 0.055f, 0.057f, 0.059f, 0.061f, 0.066f, 0.075f, 0.093f, 0.125f, 0.178f, 0.246f, 0.307f, 0.337f, 0.334f, 0.317f, 0.293f, 0.262f, 0.230f, 0.198f, 0.165f, 0.135f, 0.115f, 0.104f, 0.098f, 0.094f, 0.092f, 0.093f, 0.097f, 0.102f, 0.108f, 0.113f, 0.115f, 0.114f, 0.114f},
	{0.050f, 0.049f, 0.048f, 0.047f, 0.047f, 0.047f, 0.047f, 0.047f, 0.046f, 0.045f, 0.044f, 0.044f, 0.045f, 0.046f, 0.047f, 0.048f, 0.049f, 0.050f, 0.054f, 0.060f, 0.072f, 0.104f, 0.178f, 0.312f, 0.467f, 0.581f, 0.644f, 0.675f, 0.690f, 0.698f, 0.706f, 0.715f, 0.724f, 0.730f, 0.734f, 0.738f},
	{0.058f, 0.054f, 0.052f, 0.052f, 0.053f, 0.054f, 0.056f, 0.059f, 0.067f, 0.081f, 0.107f, 0.152f, 0.225f, 0.336f, 0.462f, 0.559f, 0.616f, 0.650f, 0.672f, 0.694f, 0.710f, 0.723f, 0.731f, 0.739f, 0.746f, 0.752f, 0.758f, 0.764f, 0.769f, 0.771f, 0.776f, 0.782f, 0.790f, 0.796f, 0.799f, 0.804f},
	{0.145f, 0.195f, 0.283f, 0.346f, 0.362f, 0.354f, 0.334f, 0.306f, 0.276f, 0.248f, 0.218f, 0.190f, 0.168f, 0.149f, 0.127f, 0.107f, 0.100f, 0.102f, 0.104f, 0.109f, 0.137f, 0.200f, 0.290f, 0.400f, 0.516f, 0.615f, 0.687f, 0.732f, 0.760f, 0.774f, 0.783f, 0.793f, 0.803f, 0.812f, 0.817f, 0.825f},
	{0.108f, 0.141f, 0.192f, 0.236f, 0.261f, 0.286f, 0.317f, 0.353f, 0.390f, 0.426f, 0.446f, 0.444f, 0.423f, 0.385f, 0.337f, 0.283f, 0.231f, 0.185f, 0.146f, 0.118f, 0.101f, 0.090f, 0.082f, 0.076f, 0.074f, 0.073f, 0.073f, 0.074f, 0.076f, 0.077f, 0.076f, 0.075f, 0.073f, 0.072f, 0.074f, 0.079f},
	{0.189f, 0.255f, 0.423f, 0.660f, 0.811f, 0.862f, 0.877f, 0.884f, 0.891f, 0.896f, 0.899f, 0.904f, 0.907f, 0.909f, 0.911f, 0.910f, 0.911f, 0.914f, 0.913f, 0.916f, 0.915f, 0.916f, 0.914f, 0.915f, 0.918f, 0.919f, 0.921f, 0.923f, 0.924f, 0.922f, 0.922f, 0.925f, 0.927f, 0.930f, 0.930f, 0.933f},
	{0.171f, 0.232f, 0.365f, 0.507f, 0.567f, 0.583f, 0.588f, 0.590f, 0.591f, 0.590f, 0.588f, 0.588f, 0.589f, 0.589f, 0.591f, 0.590f, 0.590f, 0.590f, 0.589f, 0.591f, 0.590f, 0.590f, 0.587f, 0.585f, 0.583f, 0.580f, 0.578f, 0.576f, 0.574f, 0.572f, 0.571f, 0.569f, 0.568f, 0.568f, 0.566f, 0.566f},
	{0.144f, 0.192f, 0.272f, 0.331f, 0.350f, 0.357f, 0.361f, 0.363f, 0.363f, 0.361f, 0.359f, 0.358f, 0.358f, 0.359f, 0.360f, 0.360f, 0.361f, 0.361f, 0.360f, 0.362f, 0.362f, 0.361f, 0.359f, 0.358f, 0.355f, 0.352f, 0.350f, 0.348f, 0.345f, 0.343f, 0.340f, 0.338f, 0.335f, 0.334f, 0.332f, 0.331f},
	{0.105f, 0.131f, 0.163f, 0.180f, 0.186f, 0.190f, 0.193f, 0.194f, 0.194f, 0.192f, 0.191f, 0.191f, 0.191f, 0.192f, 0.192f, 0.192f, 0.192f, 0.192f, 0.192f, 0.193f, 0.192f, 0.192f, 0.191f, 0.189f, 0.188f, 0.186f, 0.184f, 0.182f, 0.181f, 0.179f, 0.178f, 0.176f, 0.174f, 0.173f, 0.172f, 0.171f},
	{0.068f, 0.077f, 0.084f, 0.087f, 0.089f, 0.090f, 0.092f, 0.092f, 0.091f, 0.090f, 0.090f, 0.090f, 0.090f, 0.090f, 0.090f, 0.090f, 0.090f, 0.090f, 0.090f, 0.090f, 0.090f, 0.089f, 0.089f, 0.088f, 0.087f, 0.086f, 0.086f, 0.085f, 0.084f, 0.084f, 0.083f, 0.083f, 0.082f, 0.081f, 0.081f, 0.081f},
	{0.031f, 0.032f, 0.032f, 0.033f, 0.033f, 0.033f, 0.033f, 0.033f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.032f, 0.033f}
};

// [comment]
// CIE XYZ to RGB transformation matrix
// [/comment]
const double XYZ_to_RGB[][3] = {
	{ 2.3706743, -0.9000405, -0.4706338},
	{-0.5138850,  1.4253036,  0.0885814},
	{ 0.0052982, -0.0146949,  1.0093968}
};

// [comment]
// Convert XYZ color to RGB color space
// [/comment]
void XYZtoRGB(const float &X, const float &Y, const float &Z, float &r, float &g, float &b)
{
	r = std::max(0., X * XYZ_to_RGB[0][0] + Y * XYZ_to_RGB[0][1] + Z * XYZ_to_RGB[0][2]);
	g = std::max(0., X * XYZ_to_RGB[1][0] + Y * XYZ_to_RGB[1][1] + Z * XYZ_to_RGB[1][2]);
	b = std::max(0., X * XYZ_to_RGB[2][0] + Y * XYZ_to_RGB[2][1] + Z * XYZ_to_RGB[2][2]);
}

// [comment]
// Convert a spectrum to a XYZ color
// [/comment]
void spectrumToXYZ(int colorIndex, float& X, float& Y, float& Z)
{
	float S = 0;
	for (int i = 0; i < 36; ++i) {
		X += colorMatchingFunc[0][i * 2] * spectralData[colorIndex][i];
		Y += colorMatchingFunc[1][i * 2] * spectralData[colorIndex][i];
		Z += colorMatchingFunc[2][i * 2] * spectralData[colorIndex][i];
		S += colorMatchingFunc[1][i * 2];
	}
	X /= S;
	Y /= S;
	Z /= S;
}

class Image
{
public:
	Image(const int &w, const int &h) : width(w), height(h)
	{
		imageData = new float [w * h * 3];
	}
	~Image()
	{
		delete [] imageData;
	}
	void setPixel(const float *pixelValues, const int &x, const int &y)
	{
		imageData[(y * width + x) * 3] = pixelValues[0];
		imageData[(y * width + x) * 3 + 1] = pixelValues[1];
		imageData[(y * width + x) * 3 + 2] = pixelValues[2];
	}
	void saveToPpm(const char *filename)
	{
		float gamma = 1;
		std::ofstream ofs;
		ofs.open(filename, std::ios::out | std::ios::binary);
		ofs << "P6\n" << width << " " << height << "\n255\n";
		float *pixel = imageData;
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width; ++i) {
				unsigned char r = (unsigned char)(std::max(0.f, std::min(255.f, powf(pixel[0], 1 / gamma) * 255 + 0.5f)));
				unsigned char g = (unsigned char)(std::max(0.f, std::min(255.f, powf(pixel[1], 1 / gamma) * 255 + 0.5f)));
				unsigned char b = (unsigned char)(std::max(0.f, std::min(255.f, powf(pixel[2], 1 / gamma) * 255 + 0.5f)));
				ofs << r << g << b;
				pixel += 3;
			}
		}
		ofs.close();
	}
	float *imageData;
	int width, height;
};

int main(int argc, char **argv)
{
	// [comment]
	// Convert the spectrum data of each bucket from the McBeth chart to a XYZ then RGB color
	// [/comment]
	float rgb[24][3];
	for (int i = 0; i < 24; ++i) {
		float X(0), Y(0), Z(0);
		spectrumToXYZ(i, X, Y, Z);
		XYZtoRGB(X, Y, Z, rgb[i][0], rgb[i][1], rgb[i][2]);
		std::cerr << rgb[i][0] << " " << rgb[i][1] << " " << rgb[i][2] << std::endl;
		fprintf(stderr, "%d RGB %d %d %d\n", i + 1,
			(unsigned char)(255 * rgb[i][0]),
			(unsigned char)(255 * rgb[i][1]),
			(unsigned char)(255 * rgb[i][2]));
	}

	// [comment]
	// Store result to an image (fill a small bucket with the bucket's color)
	// [/comment]
	int patchSize = 64;
	int width = patchSize * 6;
	int height = patchSize * 4;
	Image image(width, height);
	for (int j = 0; j < 4; ++j) {
		int offsetj = j * patchSize;
		for (int i = 0; i < 6; ++i) {
			int offseti = i * patchSize;
			for (int jj = 0; jj < patchSize; ++jj) {
				for (int ii = 0; ii < patchSize; ++ii) {
					image.setPixel(rgb[j * 6 + i], offseti + ii, offsetj + jj);
				}
			}
		}
	}

	image.saveToPpm("./img_mcbeth.ppm");

	return 0;
}