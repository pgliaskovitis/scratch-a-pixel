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
// This program renders a wireframe image of a 3D object whose description is stored
// in the program itself. The result is stored to a SVG file. To draw an image
// of each triangle making up that object, we project the vertices of each triangle
// onto the screen using perspective projection (effectively transforming the vertices
// world coordinates to 2D pixel coordinates). Triangles are stored in the SVG
// file by connecting their respective vertices to each other with lines.
//[/header]

#include <fstream>

#include "geometry.h"

//[comment]
// List of vertices making up the object
//[/comment]
const Vec3f verts[146] = {
	{        0.0f,    39.034f,         0.0f}, {  0.76212f,    36.843f,         0.0f},
	{        3.0f,    36.604f,         0.0f}, {        1.0f,    35.604f,         0.0f},
	{   2.0162f,    33.382f,         0.0f}, {        0.0f,    34.541f,         0.0f},
	{  -2.0162f,    33.382f,         0.0f}, {       -1.0f,    35.604f,         0.0f},
	{       -3.0f,    36.604f,         0.0f}, { -0.76212f,    36.843f,         0.0f},
	{-0.040181f,     34.31f,         0.0f}, {   3.2778f,    30.464f,         0.0f},
	{-0.040181f,    30.464f,         0.0f}, {-0.028749f,    30.464f,         0.0f},
	{   3.2778f,    30.464f,         0.0f}, {   1.2722f,    29.197f,         0.0f},
	{   1.2722f,    29.197f,         0.0f}, {-0.028703f,    29.197f,         0.0f},
	{   1.2722f,    29.197f,         0.0f}, {   5.2778f,    25.398f,         0.0f},
	{ -0.02865f,    25.398f,         0.0f}, {   1.2722f,    29.197f,         0.0f},
	{   5.2778f,    25.398f,         0.0f}, {   3.3322f,    24.099f,         0.0f},
	{-0.028683f,    24.099f,         0.0f}, {   7.1957f,    20.299f,         0.0f},
	{ -0.02861f,    20.299f,         0.0f}, {   5.2778f,    19.065f,         0.0f},
	{-0.028663f,    18.984f,         0.0f}, {   9.2778f,    15.265f,         0.0f},
	{-0.028571f,    15.185f,         0.0f}, {   9.2778f,    15.265f,         0.0f},
	{   7.3772f,    13.999f,         0.0f}, {-0.028625f,    13.901f,         0.0f},
	{   9.2778f,    15.265f,         0.0f}, {   12.278f,    8.9323f,         0.0f},
	{-0.028771f,    8.9742f,         0.0f}, {   12.278f,    8.9323f,         0.0f},
	{   10.278f,    7.6657f,         0.0f}, {-0.028592f,    7.6552f,         0.0f},
	{   15.278f,    2.5994f,         0.0f}, {-0.028775f,    2.6077f,         0.0f},
	{   15.278f,    2.5994f,         0.0f}, {   13.278f,    1.3329f,         0.0f},
	{-0.028727f,    1.2617f,         0.0f}, {   18.278f,   -3.7334f,         0.0f},
	{   18.278f,   -3.7334f,         0.0f}, {   2.2722f,   -1.2003f,         0.0f},
	{-0.028727f,   -1.3098f,         0.0f}, {   4.2722f,        -5.0f,         0.0f},
	{   4.2722f,        -5.0f,         0.0f}, {-0.028727f,        -5.0f,         0.0f},
	{  -3.3582f,    30.464f,         0.0f}, {  -3.3582f,    30.464f,         0.0f},
	{  -1.3526f,    29.197f,         0.0f}, {  -1.3526f,    29.197f,         0.0f},
	{  -1.3526f,    29.197f,         0.0f}, {  -5.3582f,    25.398f,         0.0f},
	{  -1.3526f,    29.197f,         0.0f}, {  -5.3582f,    25.398f,         0.0f},
	{  -3.4126f,    24.099f,         0.0f}, {   -7.276f,    20.299f,         0.0f},
	{  -5.3582f,    19.065f,         0.0f}, {  -9.3582f,    15.265f,         0.0f},
	{  -9.3582f,    15.265f,         0.0f}, {  -7.4575f,    13.999f,         0.0f},
	{  -9.3582f,    15.265f,         0.0f}, {  -12.358f,    8.9323f,         0.0f},
	{  -12.358f,    8.9323f,         0.0f}, {  -10.358f,    7.6657f,         0.0f},
	{  -15.358f,    2.5994f,         0.0f}, {  -15.358f,    2.5994f,         0.0f},
	{  -13.358f,    1.3329f,         0.0f}, {  -18.358f,   -3.7334f,         0.0f},
	{  -18.358f,   -3.7334f,         0.0f}, {  -2.3526f,   -1.2003f,         0.0f},
	{  -4.3526f,        -5.0f,         0.0f}, {  -4.3526f,        -5.0f,         0.0f},
	{        0.0f,     34.31f,  0.040181f}, {        0.0f,    30.464f,   -3.2778f},
	{        0.0f,    30.464f,  0.040181f}, {        0.0f,    30.464f,  0.028749f},
	{        0.0f,    30.464f,   -3.2778f}, {        0.0f,    29.197f,   -1.2722f},
	{        0.0f,    29.197f,   -1.2722f}, {        0.0f,    29.197f,  0.028703f},
	{        0.0f,    29.197f,   -1.2722f}, {        0.0f,    25.398f,   -5.2778f},
	{        0.0f,    25.398f,   0.02865f}, {        0.0f,    29.197f,   -1.2722f},
	{        0.0f,    25.398f,   -5.2778f}, {        0.0f,    24.099f,   -3.3322f},
	{        0.0f,    24.099f,  0.028683f}, {        0.0f,    20.299f,   -7.1957f},
	{        0.0f,    20.299f,   0.02861f}, {        0.0f,    19.065f,   -5.2778f},
	{        0.0f,    18.984f,  0.028663f}, {        0.0f,    15.265f,   -9.2778f},
	{        0.0f,    15.185f,  0.028571f}, {        0.0f,    15.265f,   -9.2778f},
	{        0.0f,    13.999f,   -7.3772f}, {        0.0f,    13.901f,  0.028625f},
	{        0.0f,    15.265f,   -9.2778f}, {        0.0f,    8.9323f,   -12.278f},
	{        0.0f,    8.9742f,  0.028771f}, {        0.0f,    8.9323f,   -12.278f},
	{        0.0f,    7.6657f,   -10.278f}, {        0.0f,    7.6552f,  0.028592f},
	{        0.0f,    2.5994f,   -15.278f}, {        0.0f,    2.6077f,  0.028775f},
	{        0.0f,    2.5994f,   -15.278f}, {        0.0f,    1.3329f,   -13.278f},
	{        0.0f,    1.2617f,  0.028727f}, {        0.0f,   -3.7334f,   -18.278f},
	{        0.0f,   -3.7334f,   -18.278f}, {        0.0f,   -1.2003f,   -2.2722f},
	{        0.0f,   -1.3098f,  0.028727f}, {        0.0f,        -5.0f,   -4.2722f},
	{        0.0f,        -5.0f,   -4.2722f}, {        0.0f,        -5.0f,  0.028727f},
	{        0.0f,    30.464f,    3.3582f}, {        0.0f,    30.464f,    3.3582f},
	{        0.0f,    29.197f,    1.3526f}, {        0.0f,    29.197f,    1.3526f},
	{        0.0f,    29.197f,    1.3526f}, {        0.0f,    25.398f,    5.3582f},
	{        0.0f,    29.197f,    1.3526f}, {        0.0f,    25.398f,    5.3582f},
	{        0.0f,    24.099f,    3.4126f}, {        0.0f,    20.299f,     7.276f},
	{        0.0f,    19.065f,    5.3582f}, {        0.0f,    15.265f,    9.3582f},
	{        0.0f,    15.265f,    9.3582f}, {        0.0f,    13.999f,    7.4575f},
	{        0.0f,    15.265f,    9.3582f}, {        0.0f,    8.9323f,    12.358f},
	{        0.0f,    8.9323f,    12.358f}, {        0.0f,    7.6657f,    10.358f},
	{        0.0f,    2.5994f,    15.358f}, {        0.0f,    2.5994f,    15.358f},
	{        0.0f,    1.3329f,    13.358f}, {        0.0f,   -3.7334f,    18.358f},
	{        0.0f,   -3.7334f,    18.358f}, {        0.0f,   -1.2003f,    2.3526f},
	{        0.0f,        -5.0f,    4.3526f}, {        0.0f,        -5.0f,    4.3526f}
};

const uint32_t numTris = 128;

//[comment]
// Triangle index array. A triangle has 3 vertices. Each successive group of 3
// integers in this array represent the positions of the vertices in the vertex
// array making up one triangle of that object. For example, the first 3 integers
// from this array, 8/7/9 represent the positions of the vertices making up the
// the first triangle. You can access these vertices with the following code:
//
//     verts[8]; /* first vertex  */
//
//     verts[7]; /* second vertex */
//
//     verts[9]; /* third vertex  */
//
// 6/5/5 are the positions of the vertices in the vertex array making up the second
// triangle, and so on.
// To find the indices of the n-th triangle, use the following code:
//
//     tris[n * 3];     /* index of the first vertex in the verts array */
//
//     tris[n * 3 + 1]; /* index of the second vertexin the verts array */
//
//     tris[n * 3 + 2]; /* index of the third vertex in the verts array */
//[/comment]
const uint32_t tris[numTris * 3] = {
	  8,   7,   9,   6,   5,   7,   4,   3,   5,   2,   1,   3,   0,   9,   1,
	  5,   3,   7,   7,   3,   9,   9,   3,   1,  10,  12,  11,  13,  15,  14,
	 15,  13,  16,  13,  17,  16,  18,  20,  19,  17,  20,  21,  20,  23,  22,
	 20,  24,  23,  23,  26,  25,  24,  26,  23,  26,  27,  25,  26,  28,  27,
	 27,  30,  29,  28,  30,  27,  30,  32,  31,  30,  33,  32,  27,  30,  34,
	 32,  36,  35,  33,  36,  32,  36,  38,  37,  36,  39,  38,  38,  41,  40,
	 39,  41,  38,  41,  43,  42,  41,  44,  43,  44,  45,  43,  44,  47,  46,
	 44,  48,  47,  48,  49,  47,  48,  51,  50,  10,  52,  12,  13,  53,  54,
	 55,  17,  54,  13,  54,  17,  56,  57,  20,  17,  58,  20,  20,  59,  60,
	 20,  60,  24,  60,  61,  26,  24,  60,  26,  26,  61,  62,  26,  62,  28,
	 62,  63,  30,  28,  62,  30,  30,  64,  65,  30,  65,  33,  62,  66,  30,
	 65,  67,  36,  33,  65,  36,  36,  68,  69,  36,  69,  39,  69,  70,  41,
	 39,  69,  41,  41,  71,  72,  41,  72,  44,  44,  72,  73,  44,  74,  75,
	 44,  75,  48,  48,  75,  76,  48,  77,  51,  78,  80,  79,  81,  83,  82,
	 83,  81,  84,  81,  85,  84,  86,  88,  87,  85,  88,  89,  88,  91,  90,
	 88,  92,  91,  91,  94,  93,  92,  94,  91,  94,  95,  93,  94,  96,  95,
	 95,  98,  97,  96,  98,  95,  98, 100,  99,  98, 101, 100,  95,  98, 102,
	100, 104, 103, 101, 104, 100, 104, 106, 105, 104, 107, 106, 106, 109, 108,
	107, 109, 106, 109, 111, 110, 109, 112, 111, 112, 113, 111, 112, 115, 114,
	112, 116, 115, 116, 117, 115, 116, 119, 118,  78, 120,  80,  81, 121, 122,
	123,  85, 122,  81, 122,  85, 124, 125,  88,  85, 126,  88,  88, 127, 128,
	 88, 128,  92, 128, 129,  94,  92, 128,  94,  94, 129, 130,  94, 130,  96,
	130, 131,  98,  96, 130,  98,  98, 132, 133,  98, 133, 101, 130, 134,  98,
	133, 135, 104, 101, 133, 104, 104, 136, 137, 104, 137, 107, 137, 138, 109,
	107, 137, 109, 109, 139, 140, 109, 140, 112, 112, 140, 141, 112, 142, 143,
	112, 143, 116, 116, 143, 144, 116, 145, 119
};

//[comment]
// Compute the 2D pixel coordinates of a point defined in world space. This function
// requires the point original world coordinates of course, the world-to-camera
// matrix (which you can get from computing the inverse of the camera-to-world matrix,
// the matrix transforming the camera), the canvas dimension and the image width and
// height in pixels.
//
// Note that we don't check in this version of the function if the point is visible
// or not. If the absolute value of of any of the screen coordinates are greater
// that their respective maximum value (the canvas width for the x-coordinate,
// and the canvas height for the y-coordinate) then the point is not visible.
// When a SVG file is displayed to the screen, the application displaying the content
// of the file clips points and lines which are not visible. Thus we can store lines or point
// in the file whether they are visible or not. The final clipping will be done when the
// image is displayed to the screen.
//[/comment]
void computePixelCoordinates(
	const Vec3f pWorld,
	Vec2i &pRaster,
	const Matrix44f &worldToCamera,
	const float &canvasWidth,
	const float &canvasHeight,
	const uint32_t &imageWidth,
	const uint32_t &imageHeight)
{
	Vec3f pCamera;
	worldToCamera.multVecMatrix(pWorld, pCamera);
	Vec2f pScreen;
	pScreen.x = pCamera.x / -pCamera.z;
	pScreen.y = pCamera.y / -pCamera.z;
	Vec2f pNDC;
	pNDC.x = (pScreen.x + canvasWidth * 0.5f) / canvasWidth;
	pNDC.y = (pScreen.y + canvasHeight * 0.5f) / canvasHeight;
	pRaster.x = (int)(pNDC.x * imageWidth);
	pRaster.y = (int)((1 - pNDC.y) * imageHeight);
}

int main(int argc, char **argv)
{
	std::ofstream ofs;
	//[comment]
	// We will export the result to a SVG file.
	//[/comment]
	ofs.open("./persproj.svg");
	ofs << "<svg version=\"1.1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns=\"http://www.w3.org/2000/svg\" height=\"512\" width=\"512\">" << std::endl;
	//[comment]
	// We exported the camera matrix from Maya. We need to compute its inverse,
	// which is the matrix used in the computePixelCoordinates() function.
	//[/comment]
	Matrix44f cameraToWorld(0.871214f, 0.0f, -0.490904f, 0.0f, -0.192902f, 0.919559f, -0.342346f, 0.0f, 0.451415f, 0.392953f, 0.801132f, 0.0f, 14.777467f, 29.361945f, 27.993464f, 1.0f);
	Matrix44f worldToCamera = cameraToWorld.inverse();
	std::cerr << worldToCamera << std::endl;
	float canvasWidth = 2, canvasHeight = 2;
	uint32_t imageWidth = 512, imageHeight = 512;
	//[comment]
	// Loop over all the triangles making up the object, get the 3 vertices
	// making the current triangle, convert the 3 vertices wolrd coordinates
	// to 2D pixel coordinates using the computePixelCoordinates function.
	// Then finally store the result as 3 lines connecting the projected
	// vertices of the current triangle to each other.
	//[/comment]
	for (uint32_t i = 0; i < numTris; ++i) {
		const Vec3f &v0World = verts[tris[i * 3]];
		const Vec3f &v1World = verts[tris[i * 3 + 1]];
		const Vec3f &v2World = verts[tris[i * 3 + 2]];
		Vec2i v0Raster, v1Raster, v2Raster;
		computePixelCoordinates(v0World, v0Raster, worldToCamera, canvasWidth, canvasHeight, imageWidth, imageHeight);
		computePixelCoordinates(v1World, v1Raster, worldToCamera, canvasWidth, canvasHeight, imageWidth, imageHeight);
		computePixelCoordinates(v2World, v2Raster, worldToCamera, canvasWidth, canvasHeight, imageWidth, imageHeight);
		std::cerr << v0Raster << ", " << v1Raster << ", " << v2Raster << std::endl;
		ofs << "<line x1=\"" << v0Raster.x << "\" y1=\"" << v0Raster.y << "\" x2=\"" << v1Raster.x << "\" y2=\"" << v1Raster.y << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />\n";
		ofs << "<line x1=\"" << v1Raster.x << "\" y1=\"" << v1Raster.y << "\" x2=\"" << v2Raster.x << "\" y2=\"" << v2Raster.y << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />\n";
		ofs << "<line x1=\"" << v2Raster.x << "\" y1=\"" << v2Raster.y << "\" x2=\"" << v0Raster.x << "\" y2=\"" << v0Raster.y << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />\n";
	}
	ofs << "</svg>\n";
	ofs.close();

	return 0;
}
