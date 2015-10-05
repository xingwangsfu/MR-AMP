/**
 * @file imageio.h
 * @brief Implements ReadImage and WriteImage functions
 * @author Pascal Getreuer <getreuer@gmail.com>
 */
#ifndef _IMAGEIO_H_
#define _IMAGEIO_H_

#include <stdio.h>
#include <stdlib.h>

/* Portable integer types
 * 
 * Types uint8_t, uint16_t, uint32_t should be defined as
 * unsigned integer types such that
 *    uint8_t  is 8-bit,  range 0 to 255
 *    uint16_t is 16-bit, range 0 to 65535
 *    uint32_t is 32-bit, range 0 to 4294967295
 *
 * Similarly, int8_t, int16_t, int32_t should be defined as
 * signed integer types such that
 *    int8_t  is  8-bit, range        -128 to +127
 *    int16_t is 16-bit, range      -32768 to +32767
 *    int32_t is 32-bit, range -2147483648 to +2147483647
 *
 * These definitions are implemented with types __int8, __int16, and
 * __int32 under Windows and by including stdint.h under UNIX.
 */
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)

    /* Windows system: Use __intN types to define uint8_t, etc. */
    typedef unsigned __int8 uint8_t;
    typedef unsigned __int16 uint16_t;
    typedef unsigned __int32 uint32_t;
    typedef __int8 int8_t;
    typedef __int16 int16_t;
    typedef __int32 int32_t;
    
#else

    /* UNIX system: Use stdint to define uint8_t, etc. */
    #include <stdint.h>

#endif

/** @brief Limit on the maximum allowed image width or height (security). */
#define MAX_IMAGE_SIZE 10000

/* Build string macros listing the supported formats */
#ifdef LIBJPEG_SUPPORT
#define SUPPORTEDSTRING_JPEG	"/JPEG"
#else
#define SUPPORTEDSTRING_JPEG	""
#endif
#ifdef LIBPNG_SUPPORT
#define SUPPORTEDSTRING_PNG		"/PNG"
#else
#define SUPPORTEDSTRING_PNG		""
#endif
#ifdef LIBTIFF_SUPPORT
#define SUPPORTEDSTRING_TIFF	"/TIFF"
#else
#define SUPPORTEDSTRING_TIFF	""
#endif

/** String macro listing supported formats for \c ReadImage.  For example,
@code
    printf("The supported formats for reading are " READIMAGE_FORMATS_SUPPORTED ".\n"); 
@endcode
*/
#define READIMAGE_FORMATS_SUPPORTED	\
    "BMP" SUPPORTEDSTRING_JPEG SUPPORTEDSTRING_PNG SUPPORTEDSTRING_TIFF
    
/** String macro listing supported formats for \c WriteImage */
#define WRITEIMAGE_FORMATS_SUPPORTED	\
    "BMP" SUPPORTEDSTRING_JPEG SUPPORTEDSTRING_PNG SUPPORTEDSTRING_TIFF


/* Definitions for specifying image formats */
#define IMAGEIO_U8            0x0000
#define IMAGEIO_SINGLE        0x0001
#define IMAGEIO_DOUBLE        0x0002
#define IMAGEIO_STRIP_ALPHA   0x0010
#define IMAGEIO_BGRFLIP       0x0020
#define IMAGEIO_AFLIP         0x0040
#define IMAGEIO_GRAYSCALE     0x0080
#define IMAGEIO_PLANAR        0x0100
#define IMAGEIO_COLUMNMAJOR   0x0200
#define IMAGEIO_RGB           (IMAGEIO_STRIP_ALPHA)
#define IMAGEIO_BGR           (IMAGEIO_STRIP_ALPHA | IMAGEIO_BGRFLIP)
#define IMAGEIO_RGBA          0x0000
#define IMAGEIO_BGRA          (IMAGEIO_BGRFLIP)
#define IMAGEIO_ARGB          (IMAGEIO_AFLIP)
#define IMAGEIO_ABGR          (IMAGEIO_BGRFLIP | IMAGEIO_AFLIP)

#ifndef _CRT_SECURE_NO_WARNINGS
/** @brief Avoid MSVC warnings on using fopen */
#define _CRT_SECURE_NO_WARNINGS
#endif

int IdentifyImageType(char *Type, const char *FileName);

void *ReadImage(int *Width, int *Height, 
    const char *FileName, unsigned Format);

int WriteImage(void *Image, int Width, int Height, 
    const char *FileName, unsigned Format, int Quality);
    
#endif /* _IMAGEIO_H_ */
