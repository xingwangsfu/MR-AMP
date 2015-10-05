/**
 * @file tvrestorecli.c 
 * @brief TV-regularized image restoration command line interface
 * @author Pascal Getreuer <getreuer@gmail.com>
 */
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "cliio.h"
#include "tvreg.h"

#ifndef M_SQRT1_2
/** @brief The constant sqrt(1/2) */
#define M_SQRT1_2   0.70710678118654752440084436210484904
#endif


/** @brief Program parameters struct */
typedef struct
{
    /** @brief Input file name */
    const char *InputFile;
    /** @brief Output file name */
    const char *OutputFile;
    /** @brief Quality for saving JPEG images (0 to 100) */
    int JpegQuality;
    
    /** @brief tvregopt options object */
    tvregopt *Opt;
    /** @brief Lambda fidelity weight parameter */
    num Lambda;
    /** @brief Spatially-varying Lambda */
    image VaryingLambda;
    /** @brief Inpainting domain */
    image Domain;
    /** @brief Blur kernel */
    image Kernel;
} programparams;


int ParseParams(programparams *Params, int argc, const char *argv[]);
int SetLambda(programparams *Params, int ImageWidth, int ImageHeight);
int ReadLambda(programparams *Params, const char *String);
int ReadDomain(programparams *Params, const char *String);
int ReadKernel(programparams *Params, const char *OptionString);


int main(int argc, char *argv[])
{
    programparams Params;
    image f = NullImage, u = NullImage;
    
    
    if(!ParseParams(&Params, argc, (const char **)argv))
        goto Catch;
    
    /* Read the input image */
    if(!ReadImageObj(&f, Params.InputFile))
        goto Catch;
    
    if(!SetLambda(&Params, f.Width, f.Height))
        goto Catch;
    
    printf("Restoration parameters\n");
    printf("f         : [%d x %d %s]\n", 
        f.Width, f.Height, (f.NumChannels == 1) ? "grayscale" : "RGB");
    TvRegPrintOpt(Params.Opt);
#ifdef NUM_SINGLE    
    printf("datatype  : single precision float\n");
#else
    printf("datatype  : double precision float\n");
#endif
    printf("\n");
    
    if(!AllocImageObj(&u, f.Width, f.Height, f.NumChannels))
    {
        fprintf(stderr, "Out of memory.\n");
        goto Catch;
    }
    
    /* Set initial guess as u = f */
    memcpy(u.Data, f.Data, sizeof(num)*f.Width*f.Height*f.NumChannels);
    
    /* Perform restoration */
    if(!TvRestore(u.Data, f.Data, f.Width, f.Height, f.NumChannels, Params.Opt))
    {
        fprintf(stderr, "Error in computation.\n");
        goto Catch;
    }    
    
    printf("\n");
    
    /* Write output */
    if(!WriteImageObj(u, Params.OutputFile, Params.JpegQuality))    
        fprintf(stderr, "Error writing to \"%s\".\n", Params.OutputFile);
    else
        printf("Output written to \"%s\".\n", Params.OutputFile);
    
Catch:
    FreeImageObj(u);
    FreeImageObj(f); 
    FreeImageObj(Params.Domain);
    FreeImageObj(Params.VaryingLambda);
    FreeImageObj(Params.Kernel);
    TvRegFreeOpt(Params.Opt);
    return 0;
}


static void PrintHelpMessage()
{    
    printf("tvrestore, P. Getreuer 2009-2010\n"
        "Total variation regularized image denoising/deconvolution/inpainting\n\n");
    printf("Usage: tvrestore [param:value ...] input output\n\n"
        "where \"input\" and \"output\" are " READIMAGE_FORMATS_SUPPORTED " files.\n");
    printf("Please see the included PDF file for details.\n\n");
    printf("Problem parameters\n\n");
    printf("  lambda:<weight>    fidelity weight\n\n");
    printf("      lambda:<number>              constant\n");
    printf("      lambda:<file>                read spatially varying weight from file\n");
    printf("      lambda:<scalefactor>:<file>  scale weights by a constant factor\n\n");
    printf("  K:<kernel>         blur kernel for deconvolution\n\n");
    printf("      K:disk:<radius>              filled disk kernel\n");
    printf("      K:gaussian:<sigma>           Gaussian kernel\n");
    printf("      K:<file>                     read kernel from file\n\n");
    printf("  D:<file>           read inpainting domain from file\n\n");
    printf("  noise:<model>      noise model\n\n");
    printf("      noise:gaussian (or noise:l2) Gaussian noise model\n");
    printf("      noise:laplace (or noise:l1)  Laplace noise model\n");
    printf("      noise:poisson                Poisson noise model\n\n");
    printf("Algorithm parameters\n");
    printf("   tol:<number>          convergence tolerance\n"); 
    printf("   maxiter:<number>      maximum number of iterations\n");
    printf("   gamma1:<number>       constraint weight on d = grad u\n");
    printf("   gamma2:<number>       constraint weight on z = Ku\n\n");
    printf("File parameters\n");
    printf("   f:<file>              input file (alternative syntax)\n");
    printf("   u:<file>              output file (alternative syntax)\n");
#ifdef LIBJPEG_SUPPORT
    printf("   jpegquality:<number>  quality for saving JPEG images (0 to 100)\n");
#endif
    printf("\nExample: \n"
        "   tvrestore lambda:50 noise:l1 noisy.bmp denoised.bmp\n");
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
    
    Params->Lambda = 25;
    Params->VaryingLambda = NullImage;
    Params->Domain = NullImage;
    Params->Kernel = NullImage;
    Params->Opt = NULL;
        
    if(argc < 2)
    {
        PrintHelpMessage();
        return 0;
    }    
    
    if(!(Params->Opt = TvRegNewOpt()))
        return 0;
        
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
        else if(!strcmp(Param, "lambda"))
        {
            if(!Value)
            {
                fprintf(stderr, "Expected a value for option %s.\n", Param);
                return 0;
            }
            else if(!ReadLambda(Params, Value))
                return 0;
        }
        else if(!strcmp(Param, "K"))
        {
            if(!Value)
            {
                fprintf(stderr, "Expected a value for option %s.\n", Param);
                return 0;
            }
            else if(!ReadKernel(Params, Value))
                return 0;
        }
        else if(!strcmp(Param, "D"))
        {
            if(!Value)
            {
                fprintf(stderr, "Expected a value for option %s.\n", Param);
                return 0;
            }
            
            if(Params->Domain.Data)
                FreeImageObj(Params->Domain);
            
            if(!(ReadMatrixFromFile(&Params->Domain, Value, NULL)))
                return 0;
        }
        else if(!strcmp(Param, "noise"))
        {
            if(!Value)
            {
                fprintf(stderr, "Expected a value for option %s.\n", Param);
                return 0;
            }
            else if(!TvRegSetNoiseModel(Params->Opt, Value))
            {
                fprintf(stderr, "Unknown noise model \"%s\".\n", Value);
                return 0;
            }
        }
        else if(!strcmp(Param, "tol"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                TvRegSetTol(Params->Opt, NumValue);
            else
                return 0;
        }
        else if(!strcmp(Param, "maxiter"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                TvRegSetMaxIter(Params->Opt, (int)NumValue);
            else
                return 0;
        }
        else if(!strcmp(Param, "gamma1"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                TvRegSetGamma1(Params->Opt, NumValue);
            else
                return 0;
        }            
        else if(!strcmp(Param, "gamma2"))
        {
            if(CliGetNum(&NumValue, Value, Param))
                TvRegSetGamma2(Params->Opt, NumValue);
            else
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


int SetLambda(programparams *Params, int ImageWidth, int ImageHeight)
{
    const int NumPixels = ImageWidth * ImageHeight;
    num *VaryingLambda = Params->VaryingLambda.Data;
    num *Domain = Params->Domain.Data;
    num Lambda = Params->Lambda;
    int n;
    
    
    if(VaryingLambda && 
        (ImageWidth != Params->VaryingLambda.Width 
        || ImageHeight != Params->VaryingLambda.Height))
    {
        fprintf(stderr, "Size mismatch: lambda (%dx%d) does not match image size (%dx%d).\n",
            Params->VaryingLambda.Width, Params->VaryingLambda.Height, ImageWidth, ImageHeight);
        return 0;
    }
    
    if(Params->Domain.Data &&
        (ImageWidth != Params->Domain.Width 
        || ImageHeight != Params->Domain.Height))
    {
        fprintf(stderr, "Size mismatch: domain (%dx%d) does not match image size (%dx%d).\n",
            Params->Domain.Width, Params->Domain.Height, ImageWidth, ImageHeight);
        return 0;
    }
    
    if(!VaryingLambda)
    {
        if(!Params->Domain.Data)     /* Constant lambda */
            TvRegSetLambda(Params->Opt, Lambda);
        else
        {                       /* Inpainting */
            for(n = 0; n < NumPixels; n++)
                Domain[n] = (Domain[n]) ? 0 : Lambda;
            
            TvRegSetVaryingLambda(Params->Opt, Domain, 
                ImageWidth, ImageHeight);
        }
    }
    else                        /* Spatially-varying lambda */
    {
        if(Params->Domain.Data)      /* Inpainting domain also specified */
            for(n = 0; n < NumPixels; n++)
                if(Domain[n] > (num)0.5)
                    VaryingLambda[n] = 0;
     
        TvRegSetVaryingLambda(Params->Opt, VaryingLambda,
            ImageWidth, ImageHeight);
    }
    
    return 1;
}


/* "<number>", "<file>", or "<scalefactor>:<file>" */
int ReadLambda(programparams *Params, const char *String)
{
    double DoubleValue;    
    int i;        
    
    if(Params->VaryingLambda.Data)
    {
        FreeImageObj(Params->VaryingLambda);
        Params->VaryingLambda = NullImage;
    }
    
    i = ParseDouble(&DoubleValue, String);
        
    if(!i)                          /* "<file>" */
        return ReadMatrixFromFile(&Params->VaryingLambda, String, NULL);
    else if(String[i] == ':')       /* "<scalefactor>:<file>" */
    {
        if(ReadMatrixFromFile(&Params->VaryingLambda, String + i + 1, NULL))
        {
            num *VaryingLambda = Params->VaryingLambda.Data;
            num ScaleFactor = (num)DoubleValue;
            int n, NumPixels = Params->VaryingLambda.Width 
                * Params->VaryingLambda.Height;
                
            for(n = 0; n < NumPixels; n++)
                VaryingLambda[n] *= ScaleFactor;
            
            return 1;
        }
        else
            return 0;
    }
    else if(!String[i])             /* "<number>" */
    {
        Params->Lambda = (num)DoubleValue;
        return 1;
    }
    else
    {
        fprintf(stderr, "Invalid syntax, \"%s\".\n", String);
        return 0;
    }    
}


int KernelRescale(image *Kernel)
{    
    const int NumEl = Kernel->Width * Kernel->Height;
    num *Data = Kernel->Data;
    num Sum = 0;
    int n;

    for(n = 0; n < NumEl; n++)
        Sum += Data[n];

    if(Sum == 0)
    {
        fprintf(stderr, "Kernel must have nonzero sum.\n");
        return 0;
    }

    for(n = 0; n < NumEl; n++)
        Data[n] /= Sum;

    return 1;
}


int GaussianKernel(image *Kernel, num Sigma)
{
    const int R = (int)ceil(4*Sigma);
    const int W = 2*R + 1;
    const double ExpDenom = 2.0*(double)Sigma*(double)Sigma;
    num *Data;
    double Sum;
    int x, y;
    

    if(!Kernel || Sigma < 0.0 || !(AllocImageObj(Kernel, W, W, 1)))
        return 0;
    
    Data = Kernel->Data;
    
    if(Sigma == 0.0)
        Data[0] = 1;
    else
    {
        for(y = -R, Sum = 0; y <= R; y++)
            for(x = -R; x <= R; x++)
            {
                Data[(R + x) + W*(R + y)] = (num)exp(-(x*x + y*y)/ExpDenom);
                Sum += Data[(R + x) + W*(R + y)];
            }

        for(y = -R; y <= R; y++)
            for(x = -R; x <= R; x++)
                Data[(R + x) + W*(R + y)] = 
                    (num)(Data[(R + x) + W*(R + y)]/Sum);
    }
    
    return 1;
}


int DiskKernel(image *Kernel, num Radius)
{
    const int R = (int)ceil(Radius - 0.5);
    const int W = 2*R + 1;
    const int Res = 8;
    const double RadiusSqr = (double)Radius*(double)Radius;
    const double RadiusInnerSqr = 
        ((double)Radius - M_SQRT1_2)*((double)Radius - M_SQRT1_2);
    const double RadiusOuterSqr = 
        ((double)Radius + M_SQRT1_2)*((double)Radius + M_SQRT1_2);
    num *Data;
    double Sum, xl, yl, Start = -0.5 + 0.5/Res, Step = 1.0/Res;
    int c, x, y, m, n;
    

    if(!Kernel || Radius < 0.0 || !(AllocImageObj(Kernel, W, W, 1)))
        return 0;
    
    Data = Kernel->Data;
    
    if(Radius <= 0.5)
        Data[0] = 1;
    else
    {
        for(y = -R, Sum = 0; y <= R; y++)
            for(x = -R; x <= R; x++)
            {
                if(x*x + y*y <= RadiusInnerSqr)
                    c = Res*Res;
                else if(x*x + y*y > RadiusOuterSqr)
                    c = 0;
                else
                    for(n = 0, yl = y + Start, c = 0; n < Res; n++, yl += Step)
                        for(m = 0, xl = x + Start; m < Res; m++, xl += Step)
                            if(xl*xl + yl*yl <= RadiusSqr)
                                c++;
                            
                Data[(R + x) + W*(R + y)] = (num)c;
                Sum += c;
            }

        for(y = -R; y <= R; y++)
            for(x = -R; x <= R; x++)
                Data[(R + x) + W*(R + y)] =
                    (num)(Data[(R + x) + W*(R + y)]/Sum);
    }

    return 1;
}

    
int MakeNamedKernel(image *Kernel, const char *String)
{
    const char *ColonPtr;    
    num KernelParam;
    int Length;
    char KernelName[16];
    
    if(!Kernel || !(ColonPtr = strchr(String, ':')) 
        || (Length = (int)(ColonPtr - String)) > 9)
        return 0;
                
    strncpy(KernelName, String, Length);
    KernelName[Length] = '\0';
    KernelParam = (num)atof(ColonPtr + 1);
    
    if(!strcmp(KernelName, "Gaussian") || !strcmp(KernelName, "gaussian"))
        return GaussianKernel(Kernel, KernelParam);
    else if(!strcmp(KernelName, "Disk") || !strcmp(KernelName, "disk"))
        return DiskKernel(Kernel, KernelParam);
    else
        return 0;
}


/* Read a kernel: the syntax can be "disk:<radius>" or "gaussian:<sigma>".
   Otherwise, the string is assumed to be a file name. */
int ReadKernel(programparams *Params, const char *String)
{
    if(Params->Kernel.Data)
    {
        FreeImageObj(Params->Kernel);
        Params->Kernel = NullImage;
    }
    
    if(MakeNamedKernel(&Params->Kernel, String)
        || ReadMatrixFromFile(&Params->Kernel, String, KernelRescale))
    {
        TvRegSetKernel(Params->Opt, Params->Kernel.Data, 
            Params->Kernel.Width, Params->Kernel.Height);
        return 1;
    }
    else
        return 0;
}
