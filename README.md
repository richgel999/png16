# png16
Compiled with Visual Studio 2022.

Lossless half float PNG support prototype, with built-in tone mapping. Note the FOURCC will be changed to "haLf". The images written by this code are currently for test purposes only.

This is a very early prototype, but it's already extremely promising. The tone map operator is lossless and crams the entire dynamic range of the input half float image into a viewable SDR image. It does this in a very simple way (see the code in hdrpng.cpp). We are contributing this technology to the Public Domain.
