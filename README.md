# scratch-a-pixel
A repository of the code contained in the amazing Scratchapixel lessons https://www.scratchapixel.com/

This is for personal use and experimentation first and foremost. I do not own the copyright to this code, https://www.scratchapixel.com/ does. 

An initial goal is to support both GCC and MSVC. Minor issues from the original code such as floating point literals, constant definitions and occasionally missing headers as expected by MSVC have been corrected.

Specific issues found in the original code:

-For Windows machines output file streams must always be opened with the flags: 

`std::ios::out | std::ios::binary`

-Original raster3d.cpp needs a change to this line when writing to output file stream:

`ofs.write((char*)frameBuffer, imageWidth * imageHeight * (sizeof *frameBuffer));`

