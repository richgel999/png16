# png16
Compiled with Visual Studio 2022. Working on cmake/Linux support.

**Losslessly tone mapped 48-bit PNG test images are here:**
https://github.com/richgel999/png16/tree/main/bin

Lossless half float PNG support prototype, with built-in tone mapping. Note the FOURCC will be changed to "haLf". The images written by this code are currently for test purposes only.

This is a very early prototype, but it's already extremely promising. The tone map operator is lossless and crams the entire dynamic range of the input half float image into a viewable SDR image. It does this in a very simple way (see the code in hdrpng.cpp). We are contributing this technology to the Public Domain.

Here's an example image processed using this tool, with a very high dynamic range. Note you're viewing a 48bit PNG, but the actual data downloaded by your browser (the compressed 48bpp pixel data) is only displayed as SDR (typically by displaying just the high byte of the 16-bit components). *PNG is already ready for lossless half float images, we just need that "last mile" extension ("haLf") and the viewing software.*

Tone mapping is subjective, and this tonemap operator is lossless to half floats. We're hoping this triggers more research in this area, because we need the ability to interchange half float images without using enormous libraries such as EXR. More soon.

![HDR PNG](https://github.com/richgel999/png16/blob/main/bin/Tree.png)


