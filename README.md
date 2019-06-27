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

- Render timers have been added where applicable.

- Resolution of generated images has been increased to 1920 * 1080 pixels (1080p).

Suggestions on where to go from here
====================================
There are numerous internet resources touching and expanding upon the topics discussed by the first scratch-a-pixel lessons, although not as thoroughly and cohesively, one has to dig their way through. I am going to mention some of them, surely there are others I am not aware of.

1. Inigo Quilez's site, this is a treasure trove of practical techniques:
- Parallelizing ray tracing for CPUs: https://www.iquilezles.org/www/articles/cputiles/cputiles.htm
- A simple Monte Carlo path tracing architecture: https://www.iquilezles.org/www/articles/simplepathtracing/simplepathtracing.htm

2. Path tracing with SIMD instructions and GPUs: https://aras-p.info/blog/2018/03/28/Daily-Pathtracer-Part-0-Intro/

3. Peter Shirley's Ray Tracing series of books: https://github.com/petershirley/raytracinginoneweekend

4. Pixar's optimization for shooting path tracing rays: https://graphics.pixar.com/library/OrthonormalB/paper.pdf

5. Total Compedium (of global illumination): https://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf

6. Physically-based Rendering, this is the definitive reference: http://www.pbr-book.org/3ed-2018/contents.html

7. Real-time Rendering, this has many further links: http://www.realtimerendering.com/
