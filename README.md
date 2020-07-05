# scratch-a-pixel
A repository of the code contained in the amazing Scratchapixel lessons https://www.scratchapixel.com/

This is for personal use and experimentation first and foremost. I do not own the original copyright to this code, https://www.scratchapixel.com/ does.

An initial goal is to support both GCC and MSVC.

Changes compared to the original code:

- Header files utils.h, geometry.h, geometry_utils.h, objects.h now have consolidated classes, reused by various source files.

- Explicit floating point literals, constant definitions and headers have been added as expected by MSVC.

- For Windows machines output file streams are always opened with the flags:

`std::ios::out | std::ios::binary`

- Original raster3d.cpp needs a change to this line when writing to output file stream:

`ofs.write((char*)frameBuffer, imageWidth * imageHeight * (sizeof *frameBuffer));`

- Original loadPolyMeshFromFile function needs to close the file when returning a valid mesh object.

- Original acceleration.cpp has a few compilation issues and needs more extensive consolidation with the common header files.

- Render timers have been added where applicable.

- Resolution of generated images has been increased to 1920 * 1080 pixels (1080p).

Suggestions on where to go from here
====================================
There are numerous internet resources touching and expanding upon the topics discussed by the first scratch-a-pixel lessons, although not as thoroughly and cohesively, one has to dig their way through. I am going to mention some of them, surely there are others I am not aware of.

1. Inigo Quilez's site, this is a _treasure trove_ of practical techniques, e.g.:
- Parallelizing ray tracing for CPUs: https://www.iquilezles.org/www/articles/cputiles/cputiles.htm
- A simple Monte Carlo path tracing architecture: https://www.iquilezles.org/www/articles/simplepathtracing/simplepathtracing.htm
- Analytical bounding boxes for Bezier curves: https://www.iquilezles.org/www/articles/bezierbbox/bezierbbox.htm
- Smooth blending for constructive solid geometry: http://www.iquilezles.org/www/articles/smin/smin.htm
- Analytical derivatives for value noise: https://www.iquilezles.org/www/articles/morenoise/morenoise.htm

2. Peter Shirley's Ray Tracing series of books: https://github.com/petershirley/raytracinginoneweekend

3. UC Davis graduate course on ray tracing for global illumination: https://www.youtube.com/watch?v=wENIThh7XWo&list=PLslgisHe5tBPckSYyKoU3jEA4bqiFmNBJ&index=1

4. Charles University of Prague graduate course on computer graphics: https://cgg.mff.cuni.cz/~jaroslav/teaching/2015-npgr010/

5. Total Compedium (of global illumination): https://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf

6. Physically-based Rendering, this is the definitive reference: http://www.pbr-book.org/3ed-2018/contents.html

7. Real-time Rendering, this has many further links: http://www.realtimerendering.com/
