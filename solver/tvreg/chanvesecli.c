/**
* @file chanvesecli.c
* @brief Chan-Vese image segmentation command line interface
* @author Pascal Getreuer <getreuer@gmail.com>
*
* Pascal Getreuer 2010 
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cliio.h"
#include "chanvese.h"

#define ROUNDCLAMP(x)	((x < 0.0) ? 0 : \
    ((x > 1.0) ? 255 : (uint8_t)floor(255.0*(x) + 0.5)))


typedef enum
{
    SEGMENTDISPLAY_COMPOSITE,
    SEGMENTDISPLAY_BINARY,
    SEGMENTDISPLAY_CURVE,
    SEGMENTDISPLAY_INSIDE,
    SEGMENTDISPLAY_OUTSIDE
} segmentdisplay;


/** @brief Program parameters struct */
typedef struct
{
    /** @brief Input file name */
    const char *InputFile;
    /** @brief Output file name */
    const char *OutputFile;
    /** @brief Quality for saving JPEG images (0 to 100) */
    int JpegQuality;
    /** @brief Display style of the segmentation output */
    segmentdisplay Display;
    
    /** @brief Level set */
    image Phi;
    /** @brief ChanVese options object */
    chanveseopt *Opt;    
} programparams;


int WriteOutput(const image Phi, const image f, segmentdisplay Style,
    const char *FileName, int JpegQuality);
int ParseParams(programparams *Params, int argc, const char *argv[]);
int PhiRescale(image *Phi);
int ReadDisplayStyle(segmentdisplay *Display, const char *String);

int main(int argc, char *argv[])
{
    programparams Params;
    image f = NullImage;
    num c1[3], c2[3];
    int Status = 1;

    
    if(!ParseParams(&Params, argc, (const char **)argv))
        goto Catch;
    
    /* Read the input image */
    if(!ReadImageObj(&f, Params.InputFile))
        goto Catch;
    
    if(Params.Phi.Data && (f.Width != Params.Phi.Width || f.Height != Params.Phi.Height))
    {
        fprintf(stderr, "Size mismatch: phi0 (%dx%d) does not match image size (%dx%d).\n",
            Params.Phi.Width, Params.Phi.Height, f.Width, f.Height);
        goto Catch;
    }
            
    printf("Segmentation parameters\n");
    printf("f         : [%d x %d %s]\n", 
        f.Width, f.Height, (f.NumChannels == 1) ? "grayscale" : "RGB");
    printf("phi0      : %s\n", (Params.Phi.Data) ? "custom" : "default");
    ChanVesePrintOpt(Params.Opt);
#ifdef NUM_SINGLE    
    printf("datatype  : single precision float\n");
#else
    printf("datatype  : double precision float\n");
#endif    
    printf("\n");
    
    if(!Params.Phi.Data)
    {
        if(!AllocImageObj(&Params.Phi, f.Width, f.Height, 1))
        {
            fprintf(stderr, "Out of memory.");
            goto Catch;
        }
        
        ChanVeseInitPhi(Params.Phi.Data, Params.Phi.Width, Params.Phi.Height);
    }

    /* Perform the segmentation */
    if(!ChanVese(Params.Phi.Data, f.Data, 
        f.Width, f.Height, f.NumChannels, Params.Opt))
    {
        fprintf(stderr, "Error in ChanVese.");
        goto Catch;
    }    
    
    /* Compute the final region averages */
    RegionAverages(c1, c2, Params.Phi.Data, f.Data, 
        f.Width, f.Height, f.NumChannels);
    
    printf("\nRegion averages\n");
    
    if(f.NumChannels == 1)
        printf("c1        : %.4f\nc2        : %.4f\n\n", c1[0], c2[0]);
    else if(f.NumChannels == 3)
        printf("c1        : (%.4f, %.4f, %.4f)\nc2        : (%.4f, %.4f, %.4f)\n\n", 
            c1[0], c1[1], c1[2], c2[0], c2[1], c2[2]);
        
    /* Write the output image */
    if(!WriteOutput(Params.Phi, f, Params.Display, 
        Params.OutputFile, Params.JpegQuality))
    {
        fprintf(stderr, "Error writing \"%s\".\n", Params.OutputFile);
        goto Catch;
    }
    else
        printf("Output written to \"%s\".\n", Params.OutputFile);
    
    Status = 0;
Catch:
    FreeImageObj(Params.Phi);
    FreeImageObj(f);
    ChanVeseFreeOpt(Params.Opt);
    return Status;
}


int WriteOutput(const image Phi, const image f, segmentdisplay Display,
    const char *FileName, int JpegQuality)
{
    const int NumPixels = f.Width*f.Height;
    const num *Red = f.Data;
    const num *Green = f.Data + NumPixels;
    const num *Blue = f.Data + 2*NumPixels;
    num OutRed, OutGreen, OutBlue;
    uint32_t *Data = NULL;
    uint8_t *OutPtr, OutAlpha = 255;
    int x, y, n, Inside, Edge, Success = 0;
    
    
    if(!(Data = (uint32_t *)malloc(sizeof(uint32_t)*NumPixels)))
    {
        fprintf(stderr, "Out of memory.\n");
        goto Catch;
    }
    
    OutPtr = (uint8_t *)Data;
    
    for(y = 0, n = 0; y < f.Height; y++)
        for(x = 0; x < f.Width; x++, n++)
        {
            if(Phi.Data[n] >= 0)     /* Inside the curve */
            {
                Inside = 1;
                
                if(    (x > 0            && Phi.Data[n - 1] < 0)
                    || (x + 1 < f.Width  && Phi.Data[n + 1] < 0)
                    || (y > 0            && Phi.Data[n - f.Width] < 0)
                    || (y + 1 < f.Height && Phi.Data[n + f.Width] < 0))
                    Edge = 1;       /* Inside the curve, on the edge */
                else
                    Edge = 0;                
            }
            else                    /* Outside the curve */
                Inside = Edge = 0;
            
            switch(Display)
            {
            case SEGMENTDISPLAY_BINARY:
                OutRed = OutGreen = OutBlue = (num)Inside;
                break;
            case SEGMENTDISPLAY_CURVE:
                OutRed = OutGreen = OutBlue = (num)Edge;
                break;
            case SEGMENTDISPLAY_INSIDE:
                if(Inside)
                {
                    OutRed = Red[n];
                    OutGreen = Green[n];
                    OutBlue = Blue[n];
                    OutAlpha = 255;
                }
                else
                    OutRed = OutGreen = OutBlue = OutAlpha = 0;
                break;
            case SEGMENTDISPLAY_OUTSIDE:
                if(!Inside)
                {
                    OutRed = Red[n];
                    OutGreen = Green[n];
                    OutBlue = Blue[n];
                    OutAlpha = 255;
                }
                else
                    OutRed = OutGreen = OutBlue = OutAlpha = 0;
                break;
            default:
                OutRed = OutGreen = OutBlue = (num)0.8*((num)0.2989 * Red[n]
                                            + (num)0.5870 * Green[n] 
                                            + (num)0.1140 * Blue[n]);

                if(Inside)
                {
                    /* if(Edge)
                    {
                        OutRed /= 2;
                        OutGreen /= 2;
                        OutBlue = 1;
                    }
                    else */
                        OutBlue += (num)0.6;
                }
                else
                    OutRed += (num)0.6;
                break;
            }
            
            OutPtr[4*n + 0] = ROUNDCLAMP(OutRed);
            OutPtr[4*n + 1] = ROUNDCLAMP(OutGreen);
            OutPtr[4*n + 2] = ROUNDCLAMP(OutBlue);
            OutPtr[4*n + 3] = OutAlpha;
        }
        
    if(!WriteImage(Data, f.Width, f.Height, FileName,
        IMAGEIO_U8 | IMAGEIO_RGBA, JpegQuality))
        goto Catch;
    
    Success = 1;
Catch:    
    if(Data)
        free(Data);
    return Success;
}


void PrintHelpMessage()
{
    printf("chanvese, P. Getreuer 2007-2010\n"
        "Chan-Vese active contours without edges image segmentation\n\n");
    printf("Usage: chanvese [param:value ...] input output \n\n");
    printf("where \"input\" and \"output\" are " READIMAGE_FORMATS_SUPPORTED " files.\n");
    printf("Please see the included PDF file for details.\n\n");
    printf("Parameters\n\n");
    printf("   mu:<number>           length penalty (default 0.25)\n");
    printf("   nu:<number>           area penalty (default 0.0)\n");
    printf("   lambda1:<number>      fit weight inside the cuve (default 1.0)\n");
    printf("   lambda2:<number>      fit weight outside the curve (default 1.0)\n");
    printf("   phi0:<file>           read initial level set from an image or text file\n");
    printf("   tol:<number>          convergence tolerance (default 1e-4)\n");
    printf("   maxiter:<number>      maximum number of iterations (default 500)\n");    
    printf("   dt:<number>           time step (default 0.5)\n\n");
    printf("   display:<style>       display style of the output segmentation\n");
    printf("      display:composite     as a colorful overlay composited with the input\n");
    printf("      display:binary        as a binary image (white = inside, black = outside)\n");
    printf("      display:curve         as a binary image of only the curve\n");
    printf("      display:inside        as a cutout of the input image inside the curve\n");
    printf("      display:outside       as a cutout of the input image outside the curve\n\n");
#ifdef LIBJPEG_SUPPORT
    printf("   jpegquality:<number>  Quality for saving JPEG images (0 to 100)\n\n");
#endif
    printf("Example:\n");
    printf("   chanvese tol:1e-5 mu:0.5 input.bmp segmented.bmp\n");
}


int ParseParams(programparams *Params, int argc, const char *argv[])
{
    static const char *DefaultOutputFile = (char *)"out.bmp";
    const char *Param, *Value;
    num NumValue;
    char TokenBuf[256];
    int k, kread, Skip;
    
    
    /* Set parameter defaults */
    Params->InputFile = NULL;
    Params->OutputFile = DefaultOutputFile;
    Params->JpegQuality = 85;
    Params->Display = SEGMENTDISPLAY_COMPOSITE;
    Params->Phi = NullImage;
    Params->Opt = NULL;
    
    if(!(Params->Opt = ChanVeseNewOpt()))
    {
        fprintf(stderr, "Out of memory.\n");
        return 0;
    }
        
    if(argc < 2)
    {
        PrintHelpMessage();
        return 0;
    }    
        
    k = 1;
    
    while(k < argc)
    {
        Skip = (argv[k][0] == '-') ? 1 : 0;        
        kread = CliParseArglist(&Param, &Value, TokenBuf, sizeof(TokenBuf),
            k, &argv[k][Skip], argc, argv, ":");        
       
        if(!Param)
        {
            if(!Params->InputFile)
                Param = (char *)"f";
            else
                Param = (char *)"u";
        }
        
        if(Param[0] == '-')     /* Argument begins with two dashes "--" */
        {
            PrintHelpMessage();
            return 0;
        }

        if(!strcmp(Param, "f") || !strcmp(Param, "input"))
        {
            if(!Value)
            {
                fprintf(stderr, "Expected a value for option %s.\n", Param);
                return 0;
            }
            Params->InputFile = Value;
        }
        else if(!strcmp(Param, "u") || !strcmp(Param, "output"))
        {
            if(!Value)
            {
                fprintf(stderr, "Expected a value for option %s.\n", Param);
                return 0;
            }
            Params->OutputFile = Value;
        }       
        else if(!strcmp(Param, "tol"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                ChanVeseSetTol(Params->Opt, NumValue);
            else
                return 0;
        }
        else if(!strcmp(Param, "mu"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                ChanVeseSetMu(Params->Opt, NumValue);
            else
                return 0;
        }
        else if(!strcmp(Param, "nu"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                ChanVeseSetNu(Params->Opt, NumValue);
            else
                return 0;
        }
        else if(!strcmp(Param, "lambda1"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                ChanVeseSetLambda1(Params->Opt, NumValue);
            else
                return 0;
        }
        else if(!strcmp(Param, "lambda2"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                ChanVeseSetLambda2(Params->Opt, NumValue);
            else 
                return 0;
        }
        else if(!strcmp(Param, "dt"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                ChanVeseSetDt(Params->Opt, NumValue);
            else
                return 0;
        }
        else if(!strcmp(Param, "maxiter"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                ChanVeseSetMaxIter(Params->Opt, (int)NumValue);
            else
                return 0;
        }
        else if(!strcmp(Param, "phi0"))
        {
            if(!Value)
            {
                fprintf(stderr, "Expected a value for option %s.\n", Param);
                return 0;
            }
            
            if(Params->Phi.Data)
                FreeImageObj(Params->Phi);
            
            if(!(ReadMatrixFromFile(&Params->Phi, Value, PhiRescale)))
                return 0;
        }
        else if(!strcmp(Param, "display"))
        {
            if(!ReadDisplayStyle(&Params->Display, Value))
                return 0;
        }
        else if(!strcmp(Param, "jpegquality"))
        {
            if(!CliGetNum(&NumValue, Value, Param))
                return 0;
            else if(NumValue < 0 || 100 < NumValue)
            {
                fprintf(stderr, "JPEG quality must be between 0 and 100.\n");
                return 0;
            } 
            else
                Params->JpegQuality = (int)NumValue;
        }
        else if(Skip)
        {
            fprintf(stderr, "Unknown option \"%s\".\n", Param);
            return 0;
        }
        else
        {
            if(!Params->InputFile)
                Params->InputFile = argv[k];
            else
                Params->OutputFile = argv[k];
            
            kread = k;
        }

        k = kread + 1;
    }

    if(!Params->InputFile)
    {
        PrintHelpMessage();
        return 0;
    }

    return 1;
}


/* If phi is read from an image file, this function is called to rescale
   it from the range [0,1] to [-1,1].  */
int PhiRescale(image *Phi)
{    
    const int NumEl = Phi->Width * Phi->Height;
    num *Data = Phi->Data;
    int n;

    for(n = 0; n < NumEl; n++)
        Data[n] = 2*Data[n] - 1;

    return 1;
}


int ReadDisplayStyle(segmentdisplay *Display, const char *String)
{
    if(!String)
    {
        fprintf(stderr, "Expected a value for option style.\n");
        return 0;
    }
    else if(!strcmp(String, "composite"))
        *Display = SEGMENTDISPLAY_COMPOSITE;
    else if(!strcmp(String, "binary"))
        *Display = SEGMENTDISPLAY_BINARY;
    else if(!strcmp(String, "curve"))
        *Display = SEGMENTDISPLAY_CURVE;
    else if(!strcmp(String, "inside"))
        *Display = SEGMENTDISPLAY_INSIDE;
    else if(!strcmp(String, "outside"))
        *Display = SEGMENTDISPLAY_OUTSIDE;
    else
    {
        fprintf(stderr, "Unknown style \"%s\".\n", String);
        return 0;
    }
     
    return 1;
}
