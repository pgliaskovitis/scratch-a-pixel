# scratch-a-pixel
A repository of the code contained in the amazing Scratchapixel lessons https://www.scratchapixel.com/

This is for personal use and experimentation first and foremost. I do not own the original copyright to this code, https://www.scratchapixel.com/ does.

An initial goal is to support both GCC and MSVC.

Changes compared to the original code:

- Header files util.h, geometry.h, objects.h now have consolidated classes, reused by various source files.

- Explicit floating point literals, constant definitions and headers have been added as expected by MSVC.

- For Windows machines output file streams are always opened with the flags:

`std::ios::out | std::ios::binary`

- Original raster3d.cpp needs a change to this line when writing to output file stream:

`ofs.write((char*)frameBuffer, imageWidth * imageHeight * (sizeof *frameBuffer));`

- Original loadPolyMeshFromFile function needs to close the file when returning a valid mesh object.

- Render timers have been added where applicable.

- Resolution of generated images has been increased to 1920 * 1080 pixels (1080p).
