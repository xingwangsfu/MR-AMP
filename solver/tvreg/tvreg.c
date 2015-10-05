/**
 * @file tvreg.c
 * @brief TV-regularized image restoration
 * @author Pascal Getreuer <getreuer@gmail.com>
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fftw3.h>
#include "tvreg.h"

#ifdef MATLAB_MEX_FILE
    #include "mex.h"
    #define Malloc(s)    mxMalloc(s)
    #define Free(p)      mxFree(p)
    
    #ifdef MATLAB_CTRL_C
        #ifdef __cplusplus 
            extern "C" bool utIsInterruptPending();
        #else
            extern bool utIsInterruptPending();
        #endif
    #endif    
#else
    #define Malloc(s)    malloc(s)
    #define Free(p)      free(p)
#endif

#ifdef __GNUC__
    #ifndef ATTRIBUTE_UNUSED
    /** @brief Macro for the unused attribue GNU extension */
    #define ATTRIBUTE_UNUSED __attribute__((unused))
    #endif    
#else
    #define ATTRIBUTE_UNUSED
#endif

#ifndef M_PI
/** @brief The constant pi */
#define M_PI        3.14159265358979323846264338327950288
#endif
#ifndef M_2PI
/** @brief The constant 2 pi */
#define M_2PI       6.28318530717958647692528676655900577
#endif

/** @brief Size of the string buffer for holding the algorithm description */
#define ALGSTRING_SIZE  128


/*
 * Local type definitions
 */

/** @brief 2D vector type */
typedef struct
{
    num x;
    num y;
} numvec2;

/** @brief Complex value type */
typedef num numcomplex[2];

/** @brief Enum of the noise models supported by TvRestore */
typedef enum {
    NOISEMODEL_L2,
    NOISEMODEL_L1,
    NOISEMODEL_POISSON
} noisemodel;

/** @brief Options handling for TvRestore */
struct tvregstruct
{
    num Lambda;
    const num *VaryingLambda;
    int LambdaWidth;
    int LambdaHeight;
    const num *Kernel;
    int KernelWidth;
    int KernelHeight;
    num Tol;
    num Gamma1;
    num Gamma2;
    int MaxIter;
    noisemodel NoiseModel;
    int (*PlotFun)(int, int, num, const num*, int, int, int, void*);
    void *PlotParam;
    char *AlgString;
};

#ifdef __GNUC__
int TvRestoreSimplePlot(int State, int Iter, num Delta,
    __attribute__((unused)) const num *u, 
    __attribute__((unused)) int Width, 
    __attribute__((unused)) int Height, 
    __attribute__((unused)) int NumChannels,
    __attribute__((unused)) void *Param);    
#else
int TvRestoreSimplePlot(int State, int Iter, num Delta,
    const num *u, 
    int Width, 
    int Height, 
    int NumChannels,
    void *Param);
#endif

/** @brief Default options struct */
static struct tvregstruct DefaultOpt =
        {25, NULL, 0, 0, NULL, 0, 0, (num)1e-3, 5, 8, 100, 
        NOISEMODEL_L2, TvRestoreSimplePlot, NULL, NULL};

        
/*
 * Subroutines used by TvRestore
 */

/* Algorithm planning function */
static int ChooseAlgorithm(int *UseZ, int *DeconvFlag, int *DctFlag,
        void (**ZSolveFun)(num*, num*, const num*, const num*, const num*,
        int, int, int, int, int, num, num), const tvregopt *Opt);

/* Test if Kernel is symmetric */
static int IsSymmetric(const num *Kernel, int KernelWidth, int KernelHeight);

/*** d-subproblem ***/

static void DShrink(numvec2 *d, numvec2 *dtilde, const num *u, 
    int Width, int Height, int NumChannels, num Gamma1);


/*** u-subproblem ***/

static num UGaussSeidel(num *u, const num *ztilde, 
    const numvec2 *dtilde, const num *VaryingLambda,
    int Width, int Height, int NumChannels, num Lambda, num Gamma1);
static num UTrimUpdate(num *u, const num *B, int Width, int Height,
    int PadWidth, int PadHeight, int NumChannels);

/* Subroutines for DCT-based deconvolution */
static int InitDeconvDct(
    FFT(plan) *TransformA, FFT(plan) *InvTransformA, 
    FFT(plan) *TransformB, FFT(plan) *InvTransformB,
    num *KernelTrans, num *DenomTrans,
    num *A, num *ATrans, num *B, num *BTrans,
    const num *Kernel, int KernelWidth, int KernelHeight, 
    int Width, int Height, int NumChannels, num Alpha);
static void AdjBlurDct(num *ATrans, num *A, FFT(plan) TransformA,
    const num *KernelTrans, const num *ztilde,
    int Width, int Height, int NumChannels, num Alpha);
static void UTransSolveDct(num *BTrans, num *B, FFT(plan) TransformB,
    num *ATrans, const numvec2 *dtilde, 
    const num *DenomTrans, int Width, int Height, int NumChannels);
static num UDeconvDct(    
    num *u, num *B, num *BTrans, 
    FFT(plan) TransformB, FFT(plan) InvTransformB,
    num *ATrans, const numvec2 *dtilde, 
    const num *DenomTrans, int Width, int Height, int NumChannels);
static num UDeconvDctZ(
    num *u, num *A, num *ATrans, num *B, num *BTrans, 
    FFT(plan) TransformA, FFT(plan) InvTransformA,
    FFT(plan) TransformB, FFT(plan) InvTransformB,
    const numvec2 *dtilde, const num *ztilde, const num *KernelTrans,
    const num *DenomTrans, int Width, int Height, int NumChannels, num Alpha);

/* Subroutines for Fourier-based deconvolution */
static int InitDeconvFourier(
    FFT(plan) *TransformA, FFT(plan) *InvTransformA, 
    FFT(plan) *TransformB, FFT(plan) *InvTransformB,
    numcomplex *KernelTrans, num *DenomTrans,
    num *A, numcomplex *ATrans, num *B, numcomplex *BTrans, 
    const num *Kernel, int KernelWidth, int KernelHeight,
    int PadWidth, int PadHeight, int NumChannels, num Alpha);
static void AdjBlurFourier(numcomplex *ATrans, num *A, FFT(plan) TransformA,
    const numcomplex *KernelTrans, const num *ztilde,
    int Width, int Height, int NumChannels, num Alpha);
static void UTransSolveFourier(numcomplex *BTrans, num *B, FFT(plan) TransformB,
    numcomplex *ATrans, const numvec2 *dtilde, 
    const num *DenomTrans, int Width, int Height, int NumChannels);
static num UDeconvFourier(    
    num *u, num *B, numcomplex *BTrans, 
    FFT(plan) TransformB, FFT(plan) InvTransformB,
    numcomplex *ATrans, const numvec2 *dtilde, 
    const num *DenomTrans, int Width, int Height, int NumChannels);
static num UDeconvFourierZ(
    num *u, num *A, numcomplex *ATrans, num *B, numcomplex *BTrans, 
    FFT(plan) TransformA, FFT(plan) InvTransformA,
    FFT(plan) TransformB, FFT(plan) InvTransformB,
    const numvec2 *dtilde, const num *ztilde, const numcomplex *KernelTrans,
    const num *DenomTrans, int Width, int Height, int NumChannels, num Alpha);


/*** z-subproblem ***/

static void ZSolveL2(num *z, num *ztilde, 
    const num *Ku, const num *f, const num *VaryingLambda, 
    int Width, int Height, int PadWidth, int PadHeight, int NumChannels,
    num Lambda, num Gamma2);
static void ZSolveL1(num *z, num *ztilde, 
    const num *Ku, const num *f, const num *VaryingLambda, 
    int Width, int Height, int PadWidth, int PadHeight, int NumChannels,
    num Lambda, num Gamma2);
static void ZSolvePoisson(num *z, num *ztilde, 
    const num *Ku, const num *f, const num *VaryingLambda, 
    int Width, int Height, int PadWidth, int PadHeight, int NumChannels,
    num Lambda, num Gamma2);


/**
 * @brief Total variation based image restoration
 * @param u initial guess, overwritten with restored image
 * @param f input image
 * @param Width, Height, NumChannels dimensions of the input image
 * @param Opt tvregopt options object
 * @return 0 on failure, 1 on success, 2 on maximum iterations exceeded
 * 
 * This routine implements simultaneous denoising, deconvolution, and 
 * inpainting with total variation (TV) regularization, using either the 
 * Gaussian (L2), Laplace (L1), or Poisson noise model, such that Kernel*u is
 * approximately f outside of the inpainting domain.
 * 
 * The input image f should be a contiguous array of size Width by Height by
 * NumChannels in planar row-major order, 
 *    f[x + Width*(y + Height*k)] = kth component of pixel (x,y).
 * 
 * The image intensity values of f should be scaled so that the maximum 
 * intensity range of the true clean image is from 0 to 1.  It is allowed that
 * f have values outside of [0,1] (as spurious noisy pixels can exceed this
 * range), but it should be scaled so that the restored image is in [0,1].
 * This scaling is especially important for the Poisson noise model.
 * 
 * Typically, NumChannels is either 1 (grayscale image) or 3 (color image),
 * but NumChannels is allowed to be any positive integer.  If NumChannels > 1,
 * then vectorial TV (VTV) regularization is used in place of TV.
 * 
 * Image u is both an input and output of the routine.  Image u should be
 * set by the caller to an initial guess, for example a good generic 
 * initialization is to set u as a copy of f.  Image u is overwritten with the
 * restored image.
 * 
 * Other options are specified through the options object Opt.  First use 
 *    tvregopt Opt = TvRegNewOpt()
 * to create a new options object with default options (denoising with the 
 * Gaussian noise model).  Then use the following functions to make settings.
 * 
 *    - TvRegSetLambda:         set fidelity weight
 *    - TvRegSetVaryingLambda:  set spatially varying fidelity weight
 *    - TvRegSetKernel:         Kernel for deconvolution problems
 *    - TvRegSetTol:            convergence tolerance
 *    - TvRegSetMaxIter:        maximum number of iterations
 *    - TvRegSetNoiseModel:     noise model 
 *    - TvRegSetGamma1:         constraint weight on d = grad u
 *    - TvRegSetGamma2:         constraint weight on z = Ku
 *    - TvRegSetPlotFun:        custom plotting function
 * 
 * When done, call TvRegFreeOpt(Opt) to free the options object.  Setting
 * Opt = NULL uses the default options (denoising with Gaussian noise model).
 * 
 * In all filtering operations, boundaries are handled by (conceptually) 
 * extending the image with half-sample even symmetry.   For example, the 
 * sequence "abcde" is extended with half-sample even symmetry by mirroring
 * the sequence to twice the length as "abcdeedcba" and then periodizing,
 * "...cbaabcdeedc...".
 * 
 * The split Bregman method is used to solve the minimization,
 *    T. Goldstein and S. Osher,  "The Split Bregman Algorithm for L1
 *    Regularized Problems", UCLA CAM Report 08-29.
 * 
 * The routine automatically adapts the algorithm according to the inputs.  If
 * no deconvolution is needed, Gauss-Seidel is used to solve the u-subproblem.
 * If the kernel is symmetric, a DCT-based solver is used, otherwise, a 
 * (slower) Fourier-based solver is used.  For the Gaussian noise model, the 
 * routine uses a simpler splitting of the problem with two auxiliary 
 * variables.  For non-Gaussian noise models, a splitting with three auxiliary
 * variables is used.
 */
int TvRestore(num *u, const num *f, int Width, int Height, int NumChannels,
     const tvregopt *Opt)
{
    void (*ZSolveFun)(num*, num*, const num*, const num*, const num*, 
        int, int, int, int, int, num, num);
    int (*PlotFun)(int, int, num, const num*, int, int, int, void*);
    const int NumPixels = Width*Height;    
    const int NumEl = NumPixels*NumChannels;
    numvec2 *d = NULL, *dtilde = NULL;
    num *z = NULL, *ztilde = NULL;
    num *A = NULL, *ATrans = NULL, *B = NULL, *BTrans = NULL;
    num *KernelTrans = NULL, *DenomTrans = NULL;
    FFT(plan) TransformA = NULL, InvTransformA = NULL;
    FFT(plan) TransformB = NULL, InvTransformB = NULL;
    num Lambda, Gamma1, Gamma2, Tol, Alpha, Beta;
    num fNorm, DiffNorm;
    int i, PadWidth, PadHeight, TransWidth, NumTransPixels, Success = 0;
    int UseZ, DeconvFlag, DctFlag, Iter, MaxIter;
    

    if(!u || !f || u == f || Width <= 0 || Height <= 0 || NumChannels <= 0)
        return 0;
    
    /*** Set algorithm flags ***/
    
    if(!Opt)
        Opt = &DefaultOpt;
    
    if(!ChooseAlgorithm(&UseZ, &DeconvFlag, &DctFlag, &ZSolveFun, Opt))
        return 0;    
        
#if VERBOSE > 0
    TvRegPrintOpt(Opt);
#endif
    
    if(Opt->VaryingLambda && (Opt->LambdaWidth != Width || Opt->LambdaHeight != Height))
    {
        fprintf(stderr, "Size mismatch: image is %dx%d but lambda array is %dx%d.\n",
                Width, Height, Opt->LambdaWidth, Opt->LambdaHeight);
        return 0;
    }
    
    Gamma1 = Opt->Gamma1;
    Gamma2 = Opt->Gamma2;
    Tol = Opt->Tol;
    MaxIter = Opt->MaxIter;
    PlotFun = Opt->PlotFun;
    DiffNorm = (Tol > 0) ? 1000*Tol : 1000;
    Lambda = Opt->Lambda;
    Alpha = (!UseZ) ? Lambda/Gamma1 : Gamma2/Gamma1;    
    Beta = Lambda/Gamma2;
    PadWidth = Width;
    PadHeight = Height;
    
    /*** Allocate memory ***/
    
    if(!(d = (numvec2 *)Malloc(sizeof(numvec2)*NumEl))
        || !(dtilde = (numvec2 *)Malloc(sizeof(numvec2)*NumEl)))
        goto Catch;
    
    if(UseZ)
        if(!(z = (num *)Malloc(sizeof(num)*NumEl))        
            || !(ztilde = (num *)Malloc(sizeof(num)*NumEl)))
            goto Catch;
    
    if(DeconvFlag && DctFlag)
    {
        /* Prepare for DCT-based deconvolution */
        if(!(ATrans = (num *)FFT(malloc)(sizeof(num)*NumEl))
            || !(BTrans = (num *)FFT(malloc)(sizeof(num)*NumEl))
            || !(A = (num *)FFT(malloc)(sizeof(num)*NumEl))
            || !(B = (num *)FFT(malloc)(sizeof(num)*NumEl))
            || !(KernelTrans = (num *)FFT(malloc)(sizeof(num)*NumPixels))
            || !(DenomTrans = (num *)Malloc(sizeof(num)*NumPixels)))
            goto Catch;
        
        if(!InitDeconvDct(
            &TransformA, &InvTransformA, &TransformB, &InvTransformB,        
            KernelTrans, DenomTrans, 
            A, ATrans, B, BTrans,
            Opt->Kernel, Opt->KernelWidth, Opt->KernelHeight,
            Width, Height, NumChannels, Alpha))
            goto Catch;
    }
    else if(DeconvFlag)
    {
        /* Prepare for Fourier-based deconvolution */
        PadWidth = 2*Width;
        PadHeight = 2*Height;
        TransWidth = PadWidth/2 + 1;
        NumTransPixels = TransWidth*PadHeight;
        
        if(!(ATrans = (num *)FFT(malloc)(sizeof(numcomplex)*NumTransPixels*NumChannels))
            || !(BTrans = (num *)FFT(malloc)(sizeof(numcomplex)*NumTransPixels*NumChannels))
            || !(A = (num *)FFT(malloc)(sizeof(num)*PadWidth*PadHeight*NumChannels))
            || !(B = (num *)FFT(malloc)(sizeof(num)*PadWidth*PadHeight*NumChannels))
            || !(KernelTrans = (num *)FFT(malloc)(sizeof(numcomplex)*NumTransPixels))
            || !(DenomTrans = (num *)Malloc(sizeof(num)*NumTransPixels)))
            goto Catch;
                
        if(!InitDeconvFourier(
            &TransformA, &InvTransformA, &TransformB, &InvTransformB,        
            (numcomplex *)KernelTrans, DenomTrans, 
            A, (numcomplex *)ATrans, B, (numcomplex *)BTrans, 
            Opt->Kernel, Opt->KernelWidth, Opt->KernelHeight,
            PadWidth, PadHeight, NumChannels, Alpha))
            goto Catch;
    }

    /* Set convergence threshold scaled by norm of f */
    fNorm = 0;
    
    for(i = 0; i < NumEl; i++)
        fNorm += f[i]*f[i];
    
    fNorm = (num)sqrt(fNorm);
    
    if(fNorm == 0)  /* Special case, input image is zero */
    {
        memcpy(u, f, sizeof(num)*NumEl);
        Success = 1;
        goto Catch;
    }
    
    /* Initialize d = dtilde = 0 and u = z = ztilde = f */
    for(i = 0; i < NumEl; i++)
    {
        d[i].x = 0;
        d[i].y = 0;
    }
    
    for(i = 0; i < NumEl; i++)
    {
        dtilde[i].x = 0;
        dtilde[i].y = 0;
    }
    
    if(UseZ)
    {
        memcpy(z, u, sizeof(num)*NumEl);
        memcpy(ztilde, u, sizeof(num)*NumEl);
    }
    else if(DeconvFlag)        
    {    
        if(DctFlag)
            /* Compute ATrans = Alpha . KernelTrans . DCT[f] */
            AdjBlurDct(ATrans, A, TransformA, 
                KernelTrans, f, Width, Height, NumChannels, Alpha);
        else
            /* Compute ATrans = Alpha . conj(KernelTrans) . DFT[f] */
            AdjBlurFourier((numcomplex *)ATrans, A, TransformA, 
                (const numcomplex *)KernelTrans, f, Width, Height, NumChannels, Alpha);
    }    

    Success = 2;

    if(PlotFun)
        if(!PlotFun(0, 0, DiffNorm, u, Width, Height, NumChannels, Opt->PlotParam))
            goto Catch;

    /*** Algorithm main loop: Bregman iterations ***/
    for(Iter = 1; Iter <= MaxIter; Iter++)
    {
        /* Solve d subproblem and update dtilde */
        DShrink(d, dtilde, u, Width, Height, NumChannels, Gamma1);
        
        /* Solve u subproblem */
        if(!DeconvFlag)
        {
            if(!UseZ)
                DiffNorm = UGaussSeidel(u, f, 
                    dtilde, Opt->VaryingLambda, Width, Height, NumChannels,
                    Lambda, Gamma1);
            else
                DiffNorm = UGaussSeidel(u, ztilde, 
                    dtilde, NULL, Width, Height, NumChannels,
                    Gamma2, Gamma1);
        }
        else if(DctFlag)
        {   /* DCT deconvolution */
            if(!UseZ)
                DiffNorm = UDeconvDct(
                    u, B, BTrans, TransformB, InvTransformB,
                    ATrans, dtilde, DenomTrans, 
                    Width, Height, NumChannels);
            else
                /* Update u and compute A = Ku */
                DiffNorm = UDeconvDctZ(
                    u, A, ATrans, B, BTrans,
                    TransformA, InvTransformA, TransformB, InvTransformB,
                    dtilde, ztilde, KernelTrans, 
                    DenomTrans, Width, Height, NumChannels, Alpha);
        }
        else
        {   /* Fourier deconvolution */
            if(!UseZ)
                DiffNorm = UDeconvFourier(
                    u, B, (numcomplex *)BTrans, TransformB, InvTransformB,
                    (numcomplex *)ATrans, dtilde, DenomTrans, 
                    Width, Height, NumChannels);
            else
                /* Update u and compute A = Ku */
                DiffNorm = UDeconvFourierZ(
                    /* Pass a horrendous number of arguments */
                    u, A, (numcomplex *)ATrans, B, (numcomplex *)BTrans,
                    TransformA, InvTransformA, TransformB, InvTransformB,
                    dtilde, ztilde, (const numcomplex *)KernelTrans, 
                    DenomTrans, Width, Height, NumChannels, Alpha);
        }
        
        DiffNorm = (num)sqrt(DiffNorm)/fNorm;
        
        if(Iter >= 2 + UseZ && DiffNorm < Tol)
            break;
        
        /* Solve z subproblem and update ztilde */
        if(UseZ)
            ZSolveFun(z, ztilde, (!DeconvFlag) ? u : A, f, Opt->VaryingLambda,
                Width, Height, PadWidth, PadHeight, NumChannels, 
                Lambda, Gamma2);
            
        if(PlotFun)
            if(!(PlotFun(0, Iter, DiffNorm, u, 
                Width, Height, NumChannels, Opt->PlotParam)))
                goto Catch;
    }

    Success = (Iter <= MaxIter) ? 1 : 2;
    
    if(PlotFun)
        PlotFun(Success, (Iter <= MaxIter) ? Iter : MaxIter, 
            DiffNorm, u, Width, Height, NumChannels, Opt->PlotParam);
Catch:
    /* Release memory */
    if(ztilde)
        Free(ztilde);
    if(z)
        Free(z);
    if(dtilde)
        Free(dtilde);
    if(d)
        Free(d);
    if(DenomTrans)
        Free(DenomTrans);
    if(KernelTrans)
        FFT(free)(KernelTrans);
    if(B)
        FFT(free)(B);
    if(A)
        FFT(free)(A);
    if(BTrans)
        FFT(free)(BTrans);
    if(ATrans)
        FFT(free)(ATrans);
    if(DeconvFlag)
    {
        FFT(destroy_plan)(InvTransformB);
        FFT(destroy_plan)(TransformB);
        FFT(destroy_plan)(InvTransformA);
        FFT(destroy_plan)(TransformA);
        FFT(cleanup)();
    }
    return Success;
}


/** @brief Solves the d-subproblem */
static void DShrink(numvec2 *d, numvec2 *dtilde, const num *u, 
    int Width, int Height, int NumChannels, num Gamma1)
{
    const num Thresh = 1/Gamma1;
    const num ThreshSquared = Thresh*Thresh;
    const int ChannelStride = Width*Height;
    const int NumEl = NumChannels*ChannelStride;
    numvec2 dnew;
    num Magnitude;
    int i, x, y;
    
    
    for(y = 0; y < Height - 1; y++, d++, dtilde++, u++)
    {
        /* Interior points */
        for(x = 0; x < Width - 1; x++, d++, dtilde++, u++)
        {
            Magnitude = 0;
            
            for(i = 0; i < NumEl; i += ChannelStride)
            {
                d[i].x += (u[i + 1]     - u[i]) - dtilde[i].x;
                d[i].y += (u[i + Width] - u[i]) - dtilde[i].y;
                Magnitude += d[i].x*d[i].x + d[i].y*d[i].y;
            }
            
            if(Magnitude <= ThreshSquared)
            {
                for(i = 0; i < NumEl; i += ChannelStride)
                {
                    dtilde[i].x = -d[i].x;
                    dtilde[i].y = -d[i].y;
                    d[i].x = 0;
                    d[i].y = 0;
                }
            }
            else
            {
                Magnitude = 1 - Thresh/(num)sqrt(Magnitude);

                for(i = 0; i < NumEl; i += ChannelStride)
                {  
                    dnew.x = Magnitude*d[i].x;
                    dnew.y = Magnitude*d[i].y;
                    dtilde[i].x = 2*dnew.x - d[i].x;
                    dtilde[i].y = 2*dnew.y - d[i].y;
                    d[i] = dnew;
                }
            }
        }
        
        /* Right edge */
        Magnitude = 0;
            
        for(i = 0; i < NumEl; i += ChannelStride)
        {
            d[i].y += (u[i + Width] - u[i]) - dtilde[i].y;
            Magnitude += d[i].y*d[i].y;
            d[i].x = 0;
            dtilde[i].x = 0;
        }
        
        if(Magnitude <= ThreshSquared)
        {
            for(i = 0; i < NumEl; i += ChannelStride)
            {
                dtilde[i].y = -d[i].y;
                d[i].y = 0;
            }
        }
        else
        {
            Magnitude = 1 - Thresh/(num)sqrt(Magnitude);
            
            for(i = 0; i < NumEl; i += ChannelStride)
            {
                dnew.y = Magnitude*d[i].y;
                dtilde[i].y = 2*dnew.y - d[i].y;
                d[i].y = dnew.y;
            }
        }
    }
    
    /* Bottom edge */
    for(x = 0; x < Width - 1; x++, d++, dtilde++, u++)
    {
        Magnitude = 0;
            
        for(i = 0; i < NumEl; i += ChannelStride)
        {
            d[i].x += (u[i + 1] - u[i]) - dtilde[i].x;
            Magnitude += d[i].x*d[i].x;
            d[i].y = 0;
            dtilde[i].y = 0;
        }
        
        if(Magnitude <= ThreshSquared)
        {
            for(i = 0; i < NumEl; i += ChannelStride)
            {
                dtilde[i].x = -d[i].x;
                d[i].x = 0;
            }
        }
        else
        {
            Magnitude = 1 - Thresh/(num)sqrt(Magnitude);
            
            for(i = 0; i < NumEl; i += ChannelStride)
            {
                dnew.x = Magnitude*d[i].x;
                dtilde[i].x = 2*dnew.x - d[i].x;
                d[i].x = dnew.x;
            }
        }
    }
    
    /* Bottom-right corner */
    for(i = 0; i < NumEl; i += ChannelStride)
    {
        d[i].x = 0;
        d[i].y = 0;
        dtilde[i].x = 0;
        dtilde[i].y = 0;
    }
}


/** 
 * @brief Approximately solve the u-subproblem by Gauss-Seidel 
 * 
 * Performs one Gauss-Seidel iteration on u to improve the solution in the
 * u-subproblem when no deconvolution is being performed (pure denoising or
 * inpainting problems).
 */
static num UGaussSeidel(num *u, const num *ztilde, 
    const numvec2 *dtilde, const num *VaryingLambda,
    int Width, int Height, int NumChannels, num Lambda, num Gamma1)
{
    num unew, DiffNorm = 0, Temp;
    int up, down, x, y, k;
    
    
    if(!VaryingLambda)  /* Fidelity weight is constant */
    {
        const num Alpha = Lambda/Gamma1;
        const num Denom = 4 + Alpha;
        
        for(k = 0; k < NumChannels; k++)
        {
            for(y = 0; y < Height; y++)
            {
                up = (y == 0) ? 0 : -Width;
                down = (y == Height - 1) ? 0 : Width;
                
                unew = ( Alpha*ztilde[0] 
                    - dtilde[0].y + dtilde[up].y
                    + u[0] + u[1] + u[up] + u[down] ) / Denom;
                Temp = unew - u[0];
                DiffNorm += Temp*Temp;
                u[0] = unew;
                
                for(x = 1; x < Width - 1; x++)
                {
                    unew = ( Alpha*ztilde[x] 
                        - dtilde[x].x + dtilde[x - 1].x
                        - dtilde[x].y + dtilde[x + up].y
                        + u[x - 1] + u[x + 1] + u[x + up] + u[x + down] ) 
                        / Denom;
                    Temp = unew - u[x];
                    DiffNorm += Temp*Temp;
                    u[x] = unew;
                }
                
                unew = ( Alpha*ztilde[x] 
                    - dtilde[x].y + dtilde[x + up].y
                    + u[x - 1] + u[x] + u[x + up] + u[x + down] ) / Denom;
                Temp = unew - u[x];
                DiffNorm += Temp*Temp;
                u[x] = unew;
                
                u += Width;
                ztilde += Width;
                dtilde += Width;
            }
        }
    }
    else    /* Fidelity weight is spatially varying */
    {
        const num *LambdaPtr;
        num Alpha;
        
        for(k = 0; k < NumChannels; k++)
        {
            for(y = 0, LambdaPtr = VaryingLambda; y < Height; y++)
            {
                up = (y == 0) ? 0 : -Width;
                down = (y == Height - 1) ? 0 : Width;
                
                Alpha = LambdaPtr[0]/Gamma1;
                unew = ( Alpha*ztilde[0] 
                    - dtilde[0].y + dtilde[up].y
                    + u[0] + u[1] + u[up] + u[down] ) / (4 + Alpha);
                        
                Temp = unew - u[0];
                DiffNorm += Temp*Temp;
                u[0] = unew;
                
                for(x = 1; x < Width - 1; x++)
                {
                    Alpha = LambdaPtr[x]/Gamma1;
                    unew = ( Alpha*ztilde[x] 
                        - dtilde[x].x + dtilde[x - 1].x
                        - dtilde[x].y + dtilde[x + up].y
                        + u[x - 1] + u[x + 1] + u[x + up] + u[x + down] ) 
                        / (4 + Alpha);
                        
                    Temp = unew - u[x];
                    DiffNorm += Temp*Temp;
                    u[x] = unew;
                }
                
                Alpha = LambdaPtr[x]/Gamma1;
                unew = ( Alpha*ztilde[x] 
                    - dtilde[x].y + dtilde[x + up].y
                    + u[x - 1] + u[x] + u[x + up] + u[x + down] ) 
                    / (4 + Alpha);
                        
                Temp = unew - u[x];
                DiffNorm += Temp*Temp;
                u[x] = unew;
                
                u += Width;
                ztilde += Width;
                dtilde += Width;
                LambdaPtr += Width;
            }
        }
    }
        
    return DiffNorm;
}


/** @brief Trims padding, computes ||B - u||, and assigns u = B */
static num UTrimUpdate(num *u, const num *B, int Width, int Height,
    int PadWidth, int PadHeight, int NumChannels)
{
    const int PadJump = PadWidth*(PadHeight - Height);
    num unew, Temp, DiffNorm = 0;
    int x, y, k;
    
        
    for(k = 0; k < NumChannels; k++, B += PadJump)
        for(y = 0; y < Height; y++, u += Width, B += PadWidth)
            for(x = 0; x < Width; x++)
            {
                unew = B[x]; /*CLAMP(B[x], 0, 1);*/
                Temp = unew - u[x];
                DiffNorm += Temp*Temp;
                u[x] = unew;
            }
            
    return DiffNorm;
}


/**
 * @brief Boundary handling function for periodic extension
 * @param N is the data length
 * @param i is an index into the data
 * @return an index that is always between 0 and N - 1
 */
static int PeriodicExtension(int N, int i)
{
    while(1)
    {
        if(i < 0)
            i += N;
        else if(i >= N)
            i -= N;
        else
            return i;
    }
}


/**
 * @brief Boundary handling function for whole-sample symmetric extension
 * @param N is the data length
 * @param i is an index into the data
 * @return an index that is always between 0 and N - 1
 */
static int WSymExtension(int N, int i)
{
    while(1)
    {
        if(i < 0)
            i = -i;
        else if(i >= N)        
            i = (2*N - 2) - i;
        else
            return i;
    }
}


/** @brief Intializations to prepare TvRestore for DCT-based deconvolution */
static int InitDeconvDct(
    FFT(plan) *TransformA, FFT(plan) *InvTransformA, 
    FFT(plan) *TransformB, FFT(plan) *InvTransformB,
    num *KernelTrans, num *DenomTrans,
    num *A, num *ATrans, num *B, num *BTrans,
    const num *Kernel, int KernelWidth, int KernelHeight, 
    int Width, int Height, int NumChannels, num Alpha)
{
    const int NumPixels = Width*Height;
    FFT(plan) Plan = NULL;
    int Size[2];
    FFT(r2r_kind) Kind[2];
    int i, x0, y0, x, y, xi, yi;
    
    
    for(i = 0; i < NumPixels; i++)
        B[i] = 0;
    
    x0 = -KernelWidth/2;
    y0 = -KernelHeight/2;
    
    /* Pad Kernel to size Width by Height.  If Kernel
       happens to be larger, it is folded. */
    for(y = 0; y < y0 + KernelHeight; y++)
    {
        yi = WSymExtension(Height, y);
    
        for(x = 0; x < x0 + KernelWidth; x++)
        {
            xi = WSymExtension(Width, x);
            B[xi + Width*yi] += Kernel[(x - x0) + KernelWidth*(y - y0)];
        }
    }
    
    /* Compute the DCT-I transform of the padded Kernel */
    if(!(Plan = FFT(plan_r2r_2d)(Height, Width, B, KernelTrans,
        FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE | FFTW_DESTROY_INPUT)))
        return 0;
    
    FFT(execute)(Plan);
    FFT(destroy_plan)(Plan);
    
    /* Precompute the denominator that will be used in the u-subproblem. */    
    for(y = i = 0; y < Height; y++)
        for(x = 0; x < Width; x++, i++)
            DenomTrans[i] = 
                (num)(4*NumPixels*(Alpha*KernelTrans[i]*KernelTrans[i]
                + 2*(2 - cos(x*M_PI/Width) - cos(y*M_PI/Height))));

    /* Plan DCT-II transforms */
    Size[1] = Width;
    Size[0] = Height;
    Kind[0] = Kind[1] = FFTW_REDFT10;
    
    if(!(*TransformA = FFT(plan_many_r2r)(2, Size, NumChannels, 
        A, NULL, 1, NumPixels, ATrans, NULL, 1, NumPixels, Kind,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT))
        || !(*TransformB = FFT(plan_many_r2r)(2, Size, NumChannels, 
        B, NULL, 1, NumPixels, BTrans, NULL, 1, NumPixels, Kind,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT)))
        return 0;
    
    /* Plan inverse DCT-II transforms (DCT-III) */
    Kind[0] = Kind[1] = FFTW_REDFT01;
    
    if(!(*InvTransformA = FFT(plan_many_r2r)(2, Size, NumChannels, 
        ATrans, NULL, 1, NumPixels, A, NULL, 1, NumPixels, Kind,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT))
        || !(*InvTransformB = FFT(plan_many_r2r)(2, Size, NumChannels, 
        BTrans, NULL, 1, NumPixels, B, NULL, 1, NumPixels, Kind,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT)))
        return 0;
    
    return 1;
}


/** @brief Compute ATrans = Alpha . KernelTrans . DCT[ztilde] */
static void AdjBlurDct(num *ATrans, num *A, FFT(plan) TransformA,
    const num *KernelTrans, const num *ztilde,
    int Width, int Height, int NumChannels, num Alpha)
{
    const int NumPixels = Width*Height;
    int i, k;
    
    memcpy(A, ztilde, sizeof(num)*NumPixels*NumChannels);
    
    /* Compute ATrans = DCT[A] */
    FFT(execute)(TransformA);
    
    /* Compute ATrans = Alpha . KernelTrans . ATrans */
    for(k = 0; k < NumChannels; k++, ATrans += NumPixels)
        for(i = 0; i < NumPixels; i++)
            ATrans[i] = Alpha * KernelTrans[i] * ATrans[i];
}


/** @brief Compute BTrans = ( ATrans - DCT[div(dtilde)] ) / DenomTrans */
static void UTransSolveDct(num *BTrans, num *B, FFT(plan) TransformB,
    num *ATrans, const numvec2 *dtilde, 
    const num *DenomTrans, int Width, int Height, int NumChannels)
{
    const int NumPixels = Width*Height;
    int up, x, y, k;
    
    /* Compute B = div(dtilde) */
    for(k = 0; k < NumChannels; k++)
        for(y = 0; y < Height; y++, B += Width, dtilde += Width)
        {
            up = (y == 0) ? 0 : -Width;
            B[0] = dtilde[0].y - dtilde[up].y;
            
            for(x = 1; x < Width; x++)
                B[x] = dtilde[x].x - dtilde[x - 1].x
                    + dtilde[x].y - dtilde[x + up].y;
        }
    
    /* Compute BTrans = DCT[B] */
    FFT(execute)(TransformB);

    /* Compute BTrans = ( ATrans - BTrans ) / DenomTrans */
    for(k = 0; k < NumChannels; k++, ATrans += NumPixels, BTrans += NumPixels)
        for(x = 0; x < NumPixels; x++)
            BTrans[x] = (ATrans[x] - BTrans[x]) / DenomTrans[x];
}


/** @brief Compute u = IDCT[ ( ATrans - DCT[div(dtilde)] ) / DenomTrans ] */
static num UDeconvDct(    
    num *u, num *B, num *BTrans, 
    FFT(plan) TransformB, FFT(plan) InvTransformB,
    num *ATrans, const numvec2 *dtilde, 
    const num *DenomTrans, int Width, int Height, int NumChannels)
{
    /* BTrans = ( ATrans - DCT[div(dtilde)] ) / DenomTrans */
    UTransSolveDct(BTrans, B, TransformB, ATrans, dtilde,
        DenomTrans, Width, Height, NumChannels);
    /* B = IDCT[BTrans] */
    FFT(execute)(InvTransformB);
    /* Trim padding, compute ||B - u||, and assign u = B */
    return UTrimUpdate(u, B, Width, Height, Width, Height, NumChannels);
}


/** 
 * @brief Extension of UDeconvDct for when UseZ = 1
 *
 * This extended version of UDeconvDct is used when performing DCT-based 
 * deconvolution with the three-auxiliary variable algorithm (UseZ = 1),
 * that is, in a deconvolution problem with a symmetric kernel and non-
 * Gaussian noise model.
 */
static num UDeconvDctZ(
    num *u, num *A, num *ATrans, num *B, num *BTrans, 
    FFT(plan) TransformA, FFT(plan) InvTransformA,
    FFT(plan) TransformB, FFT(plan) InvTransformB,
    const numvec2 *dtilde, const num *ztilde, const num *KernelTrans,
    const num *DenomTrans, int Width, int Height, int NumChannels, num Alpha)
{
    const int NumPixels = Width*Height;
    int x, k;
    
    
    /* Compute ATrans = Alpha . KernelTrans . DCT[ztilde] */
    AdjBlurDct(ATrans, A, TransformA, KernelTrans, ztilde,
        Width, Height, NumChannels, Alpha);
    /* BTrans = ( ATrans - DCT[div(dtilde)] ) / DenomTrans */
    UTransSolveDct(BTrans, B, TransformB, ATrans, dtilde,
        DenomTrans, Width, Height, NumChannels);
    
    /* Compute ATrans = KernelTrans . BTrans */
    for(k = 0; k < NumChannels; k++, ATrans += NumPixels, BTrans += NumPixels)
        for(x = 0; x < NumPixels; x++)
            ATrans[x] = KernelTrans[x] * BTrans[x];
    
    /* A = IDCT[ATrans] = new Ku */
    FFT(execute)(InvTransformA);
    /* B = IDCT[BTrans] = new u */
    FFT(execute)(InvTransformB);
    
    /* Trim padding, compute ||B - u||, and assign u = B */
    return UTrimUpdate(u, B, Width, Height, Width, Height, NumChannels);
}


/** @brief Intializations to prepare TvRestore for Fourier deconvolution */
static int InitDeconvFourier(
    FFT(plan) *TransformA, FFT(plan) *InvTransformA, 
    FFT(plan) *TransformB, FFT(plan) *InvTransformB,
    numcomplex *KernelTrans, num *DenomTrans,
    num *A, numcomplex *ATrans, num *B, numcomplex *BTrans, 
    const num *Kernel, int KernelWidth, int KernelHeight,
    int PadWidth, int PadHeight, int NumChannels, num Alpha)
{
    const int PadNumPixels = PadWidth*PadHeight;
    const int TransWidth = PadWidth/2 + 1;    
    FFT(plan) Plan = NULL;
    int PadSize[2];
    int i, x0, y0, x, y, xi, yi;
    
    
    for(i = 0; i < PadNumPixels; i++)
        B[i] = 0;
    
    x0 = -KernelWidth/2;
    y0 = -KernelHeight/2;
    
    /* Pad Kernel to size PadWidth by PadHeight.  If Kernel
       happens to be larger, it is wrapped. */
    for(y = y0, i = 0; y < y0 + KernelHeight; y++)
    {
        yi = PeriodicExtension(PadHeight, y);
    
        for(x = x0; x < x0 + KernelWidth; x++, i++)
        {
            xi = PeriodicExtension(PadWidth, x);
            B[xi + PadWidth*yi] += Kernel[i];
        }
    }
    
    /* Compute the Fourier transform of the padded Kernel */
    if(!(Plan = FFT(plan_dft_r2c_2d)(PadHeight, PadWidth, B, 
        KernelTrans, FFTW_ESTIMATE | FFTW_DESTROY_INPUT)))
        return 0;

    FFT(execute)(Plan);
    FFT(destroy_plan)(Plan);
    
    /* Precompute the denominator that will be used in the u-subproblem. */    
    for(y = i = 0; y < PadHeight; y++)
        for(x = 0; x < TransWidth; x++, i++)
            DenomTrans[i] = 
                (num)(PadNumPixels*(Alpha*(KernelTrans[i][0]*KernelTrans[i][0]
                + KernelTrans[i][1]*KernelTrans[i][1])
                + 2*(2 - cos(x*M_2PI/PadWidth) - cos(y*M_2PI/PadHeight))));
    
    /* Plan Fourier transforms */
    PadSize[1] = PadWidth;
    PadSize[0] = PadHeight;

    if(!(*TransformA = FFT(plan_many_dft_r2c)(2, PadSize, NumChannels, 
        A, NULL, 1, PadNumPixels, ATrans, NULL, 1, TransWidth*PadHeight,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT))
        || !(*InvTransformA = FFT(plan_many_dft_c2r)(2, PadSize, NumChannels, 
        ATrans, NULL, 1, TransWidth*PadHeight, A, NULL, 1, PadNumPixels,         
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT))
        || !(*TransformB = FFT(plan_many_dft_r2c)(2, PadSize, NumChannels, 
        B, NULL, 1, PadNumPixels, BTrans, NULL, 1, TransWidth*PadHeight,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT))
        || !(*InvTransformB = FFT(plan_many_dft_c2r)(2, PadSize, NumChannels, 
        BTrans, NULL, 1, TransWidth*PadHeight, B, NULL, 1, PadNumPixels,         
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT)))
        return 0;
    
    return 1;
}


/** @brief Compute ATrans = Alpha . conj(KernelTrans) . DFT[ztilde] */
static void AdjBlurFourier(numcomplex *ATrans, num *A, FFT(plan) TransformA,
    const numcomplex *KernelTrans, const num *ztilde,
    int Width, int Height, int NumChannels, num Alpha)
{
    const int PadWidth = 2*Width;
    const int PadHeight = 2*Height;    
    const int TransWidth = PadWidth/2 + 1;
    const int TransNumPixels = TransWidth*PadHeight;
    num Temp;
    int i, iUpper, iLower, up, x, xr, y, k;
    
    
    /* Pad ztilde with even half-sample symmetry */
    for(k = 0; k < NumChannels; k++)
        for(y = 0; y < Height; y++, ztilde += Width)
        {
            up = (y == 0) ? 0 : -Width;
            iUpper = PadWidth*(y + PadHeight*k);
            iLower = PadWidth*((PadHeight - 1 - y) + PadHeight*k);   
            xr = PadWidth - 1;
            
            A[iUpper] = 
                A[iLower] = 
                A[xr + iUpper] = 
                A[xr + iLower] = ztilde[0];
            
            for(x = 1, xr--; x < Width; x++, xr--)
                A[x + iUpper] = 
                    A[x + iLower] = 
                    A[xr + iUpper] = 
                    A[xr + iLower] = ztilde[x];
        }
    
    /* Compute ATrans = DFT[A] */
    FFT(execute)(TransformA);
    
    /* Compute ATrans = Alpha . conj(KernelTrans) . ATrans */
    for(k = 0; k < NumChannels; k++, ATrans += TransNumPixels)
        for(i = 0; i < TransNumPixels; i++)
        {
            Temp = Alpha*(KernelTrans[i][0] * ATrans[i][1] 
                - KernelTrans[i][1] * ATrans[i][0]);
            ATrans[i][0] = Alpha*(KernelTrans[i][0] * ATrans[i][0] 
                + KernelTrans[i][1] * ATrans[i][1]);
            ATrans[i][1] = Temp;
        }
}


/** @brief Compute BTrans = ( ATrans - DFT[div(dtilde)] ) / DenomTrans */
static void UTransSolveFourier(numcomplex *BTrans, num *B, FFT(plan) TransformB,
    numcomplex *ATrans, const numvec2 *dtilde, 
    const num *DenomTrans, int Width, int Height, int NumChannels)
{
    const int PadWidth = 2*Width;
    const int PadHeight = 2*Height;    
    const int TransWidth = PadWidth/2 + 1;
    const int TransNumPixels = TransWidth*PadHeight;
    int iUpper, iLower, up, x, xr, y, k;
    
    /* Compute B = div(dtilde) and pad with even half-sample symmetry */
    for(k = 0; k < NumChannels; k++)
        for(y = 0; y < Height; y++, dtilde += Width)
        {
            up = (y == 0) ? 0 : -Width;
            iUpper = PadWidth*(y + PadHeight*k);
            iLower = PadWidth*((PadHeight - 1 - y) + PadHeight*k);   
            xr = PadWidth - 1;
            
            B[iUpper] = 
            B[iLower] = 
            B[xr + iUpper] = 
            B[xr + iLower] = dtilde[0].y - dtilde[up].y;     
            
            for(x = 1, xr--; x < Width; x++, xr--)
                B[x + iUpper] = 
                B[x + iLower] = 
                B[xr + iUpper] = 
                B[xr + iLower] = dtilde[x].x - dtilde[x - 1].x
                    + dtilde[x].y - dtilde[x + up].y;
        }
    
    /* Compute BTrans = DFT[B] */
    FFT(execute)(TransformB);

    /* Compute BTrans = ( ATrans - BTrans ) / DenomTrans */
    for(k = 0; k < NumChannels; k++, ATrans += TransNumPixels, BTrans += TransNumPixels)
        for(x = 0; x < TransNumPixels; x++)
        {
            BTrans[x][0] = (ATrans[x][0] - BTrans[x][0]) / DenomTrans[x];
            BTrans[x][1] = (ATrans[x][1] - BTrans[x][1]) / DenomTrans[x];
        }
}


/** @brief Compute u = IDFT[ ( ATrans - DFT[div(dtilde)] ) / DenomTrans ] */
static num UDeconvFourier(    
    num *u, num *B, numcomplex *BTrans, 
    FFT(plan) TransformB, FFT(plan) InvTransformB,
    numcomplex *ATrans, const numvec2 *dtilde, 
    const num *DenomTrans, int Width, int Height, int NumChannels)
{
    /* BTrans = ( ATrans - DFT[div(dtilde)] ) / DenomTrans */
    UTransSolveFourier(BTrans, B, TransformB, ATrans, dtilde,
        DenomTrans, Width, Height, NumChannels);
    /* B = IDFT[BTrans] */
    FFT(execute)(InvTransformB);
    /* Trim padding, compute ||B - u||, and assign u = B */
    return UTrimUpdate(u, B, Width, Height, 2*Width, 2*Height, NumChannels);
}


/** 
 * @brief Extension of UDeconvFourier for when UseZ = 1
 *
 * This extended version of UDeconvFourier is used when performing Fourier-
 * based deconvolution with the three-auxiliary variable algorithm (UseZ = 1),
 * that is, in a deconvolution problem with a non-symmetric kernel and non-
 * Gaussian noise model.
 */
static num UDeconvFourierZ(
    num *u, num *A, numcomplex *ATrans, num *B, numcomplex *BTrans, 
    FFT(plan) TransformA, FFT(plan) InvTransformA,
    FFT(plan) TransformB, FFT(plan) InvTransformB,
    const numvec2 *dtilde, const num *ztilde, const numcomplex *KernelTrans,
    const num *DenomTrans, int Width, int Height, int NumChannels, num Alpha)
{
    const int PadWidth = 2*Width;
    const int PadHeight = 2*Height;    
    const int TransWidth = PadWidth/2 + 1;
    const int TransNumPixels = TransWidth*PadHeight;
    int x, k;
    
    
    /* Compute ATrans = Alpha . conj(KernelTrans) . DFT[ztilde] */
    AdjBlurFourier(ATrans, A, TransformA, KernelTrans, ztilde,
        Width, Height, NumChannels, Alpha);
    /* BTrans = ( ATrans - DFT[div(dtilde)] ) / DenomTrans */
    UTransSolveFourier(BTrans, B, TransformB, ATrans, dtilde,
        DenomTrans, Width, Height, NumChannels);
    
    /* Compute ATrans = KernelTrans . BTrans */
    for(k = 0; k < NumChannels; k++, 
        ATrans += TransNumPixels, BTrans += TransNumPixels)
        for(x = 0; x < TransNumPixels; x++)
        {
            ATrans[x][0] = KernelTrans[x][0] * BTrans[x][0]
                - KernelTrans[x][1] * BTrans[x][1];
            ATrans[x][1] = KernelTrans[x][0] * BTrans[x][1]
                + KernelTrans[x][1] * BTrans[x][0];
        }
    
    /* A = IDFT[ATrans] = new Ku */
    FFT(execute)(InvTransformA);
    /* B = IDFT[BTrans] = new u */
    FFT(execute)(InvTransformB);
    
    /* Trim padding, compute ||B - u||, and assign u = B */
    return UTrimUpdate(u, B, Width, Height, 2*Width, 2*Height, NumChannels);
}


/** 
 * @brief Solve the z-subproblem with a Gaussian (L2) noise model 
 *
 * This routine is not needed since the three-variable auxiliary model is only
 * applied for non-Gaussian noise models.  This routine is here for testing
 * purposes.
 */
static void ZSolveL2(num *z, num *ztilde, 
    const num *Ku, const num *f, const num *VaryingLambda, 
    int Width, int Height, int PadWidth, int PadHeight, int NumChannels,
    num Lambda, num Gamma2)
{
    const int PadJump = PadWidth*(PadHeight - Height);
    num znew;
    int x, y, k;
    
    
    if(!VaryingLambda) /* Constant fidelity weight */
    {
        const num Beta = Lambda/Gamma2;
        const num Denom = 1 + Beta;
        
        for(k = 0; k < NumChannels; k++, Ku += PadJump)
            for(y = 0; y < Height; y++)
            {
                for(x = 0; x < Width; x++)
                {
                    znew = (Ku[x] + z[x] - ztilde[x] + Beta*f[x]) / Denom;
                    ztilde[x] += 2*znew - z[x] - Ku[x];
                    z[x] = znew;
                }
                
                z += Width;
                ztilde += Width;
                f += Width;
                Ku += PadWidth;
            }
    }
    else    /* Spatially varying fidelity weight */
    {
        const num *LambdaPtr;
        num Beta;
    
        for(k = 0; k < NumChannels; k++, Ku += PadJump)
            for(y = 0, LambdaPtr = VaryingLambda; y < Height; y++)
            {
                for(x = 0; x < Width; x++)
                {
                    Beta = LambdaPtr[x]/Gamma2;
                    znew = (Ku[x] + z[x] - ztilde[x] + Beta*f[x]) / (1 + Beta);
                    
                    ztilde[x] += 2*znew - z[x] - Ku[x];
                    z[x] = znew;
                }
                
                z += Width;
                ztilde += Width;
                f += Width;
                LambdaPtr += Width;
                Ku += PadWidth;
            }
    }
}


/** @brief Solve the z-subproblem with a Laplace (L1) noise model */
static void ZSolveL1(num *z, num *ztilde, 
    const num *Ku, const num *f, const num *VaryingLambda, 
    int Width, int Height, int PadWidth, int PadHeight, int NumChannels,
    num Lambda, num Gamma2)
{
    const int PadJump = PadWidth*(PadHeight - Height);
    num znew;
    int x, y, k;
    

    if(!VaryingLambda) /* Constant fidelity weight */
    {
        const num Beta = Lambda/Gamma2;
        
        for(k = 0; k < NumChannels; k++, Ku += PadJump)
            for(y = 0; y < Height; y++)
            {
                for(x = 0; x < Width; x++)
                {
                    znew = Ku[x] - f[x] + z[x] - ztilde[x];
            
                    if(znew > Beta)
                        znew += f[x] - Beta;
                    else if(znew < -Beta)
                        znew += f[x] + Beta;
                    else
                        znew = f[x];
                    
                    ztilde[x] += 2*znew - z[x] - Ku[x];
                    z[x] = znew;
                }
                
                z += Width;
                ztilde += Width;
                f += Width;
                Ku += PadWidth;
            }
    }
    else    /* Spatially varying fidelity weight */
    {
        const num *LambdaPtr;
        num Beta;
        
        for(k = 0; k < NumChannels; k++, Ku += PadJump)
            for(y = 0, LambdaPtr = VaryingLambda; y < Height; y++)
            {
                for(x = 0; x < Width; x++)
                {
                    Beta = LambdaPtr[x]/Gamma2;
                    znew = Ku[x] - f[x] + z[x] - ztilde[x];
                
                    if(znew > Beta)
                        znew += f[x] - Beta;
                    else if(znew < -Beta)
                        znew += f[x] + Beta;
                    else
                        znew = f[x];
                    
                    ztilde[x] += 2*znew - z[x] - Ku[x];
                    z[x] = znew;
                }
                
                z += Width;
                ztilde += Width;
                f += Width;
                LambdaPtr += Width;
                Ku += PadWidth;
            }
    }
}


/** @brief Solve the z-subproblem with a Poisson noise model */
static void ZSolvePoisson(num *z, num *ztilde, 
    const num *Ku, const num *f, const num *VaryingLambda, 
    int Width, int Height, int PadWidth, int PadHeight, int NumChannels,
    num Lambda, num Gamma2)
{
    const int PadJump = PadWidth*(PadHeight - Height);
    num znew;
    int x, y, k;
    
    
    if(!VaryingLambda) /* Constant fidelity weight */
    {
        const num Beta = Lambda/Gamma2;
        
        for(k = 0; k < NumChannels; k++, Ku += PadJump)
            for(y = 0; y < Height; y++)
            {
                for(x = 0; x < Width; x++)
                {
                    znew = (Ku[x] + z[x] - ztilde[x] - Beta)/2;
                    znew = znew + (num)sqrt(znew*znew + Beta*f[x]);
                            
                    ztilde[x] += 2*znew - z[x] - Ku[x];
                    z[x] = znew;
                }
                
                z += Width;
                ztilde += Width;
                f += Width;
                Ku += PadWidth;
            }
    }
    else    /* Spatially varying fidelity weight */
    {
        const num *LambdaPtr;
        num Beta;
        
        for(k = 0; k < NumChannels; k++, Ku += PadJump)
            for(y = 0, LambdaPtr = VaryingLambda; y < Height; y++)
            {
                for(x = 0; x < Width; x++)
                {
                    Beta = LambdaPtr[x]/Gamma2;
                    znew = (Ku[x] + z[x] - ztilde[x] - Beta)/2;
                    znew = znew + (num)sqrt(znew*znew + Beta*f[x]);
                            
                    ztilde[x] += 2*znew - z[x] - Ku[x];
                    z[x] = znew;
                }
                
                z += Width;
                ztilde += Width;
                f += Width;
                LambdaPtr += Width;
                Ku += PadWidth;
            }
    }
}


/** @brief Algorithm planning function */
static int ChooseAlgorithm(int *UseZ, int *DeconvFlag, int *DctFlag,
        void (**ZSolveFun)(num*, num*, const num*, const num*, const num*,
        int, int, int, int, int, num, num), const tvregopt *Opt)
{
    if(!Opt)
        return 0;        
    
    switch(Opt->NoiseModel)
    {
    case NOISEMODEL_L2:
        /* Use simplified algorithm (UseZ = 0) for the L2 model */
        *UseZ = 0;
        *ZSolveFun = ZSolveL2;
        break;
    case NOISEMODEL_L1:
        *UseZ = 1;
        *ZSolveFun = ZSolveL1;
        break;
    case NOISEMODEL_POISSON:
        *UseZ = 1;
        *ZSolveFun = ZSolvePoisson;
        break;
    default:
        return 0;
    }
    
    if(Opt->Kernel)
    {        
        if(Opt->VaryingLambda)
            *UseZ = 1;
        
        *DeconvFlag = 1;
        *DctFlag = IsSymmetric(Opt->Kernel, 
            Opt->KernelWidth, Opt->KernelHeight);
    }
    else
        *DeconvFlag = *DctFlag = 0;
    
    return 1;
}


/** @brief Test if Kernel is whole-sample symmetric */
static int IsSymmetric(const num *Kernel, int KernelWidth, int KernelHeight)
{
    int x, xr, y, yr;
        
    if(KernelWidth % 2 == 0 || KernelHeight % 2 == 0)
        return 0;
    
    for(y = 0, yr = KernelHeight - 1; y < KernelHeight; y++, yr--)
        for(x = 0, xr = KernelWidth - 1; x < KernelWidth; x++, xr--)
            if(Kernel[x + KernelWidth*y] != Kernel[xr + KernelWidth*y]
                || Kernel[x + KernelWidth*y] != Kernel[x + KernelWidth*yr])
                return 0;
    
    return 1;            
}


/* If GNU C language extensions are available, apply the "unused" attribute
   to avoid warnings.  TvRestoreSimplePlot is a plotting callback function
   for TvRestore, so the unused arguments are indeed required. */
#ifdef __GNUC__
int TvRestoreSimplePlot(int State, int Iter, num Delta,
    __attribute__((unused)) const num *u, 
    __attribute__((unused)) int Width, 
    __attribute__((unused)) int Height, 
    __attribute__((unused)) int NumChannels,
    __attribute__((unused)) void *Param)    
#else
int TvRestoreSimplePlot(int State, int Iter, num Delta,
    const num *u, 
    int Width, 
    int Height, 
    int NumChannels,
    void *Param)
#endif
{
    switch(State)
    {
    case 0: /* TvRestore is running */
        /* We print to stderr so that messages are displayed on the console
           immediately, during the TvRestore computation.  If we use stdout,
           messages might be buffered and not displayed until after TvRestore
           completes, which would defeat the point of having this real-time 
           plot callback. */
        fprintf(stderr, "   Iteration %4d     Delta %7.4f\r", Iter, Delta);
        break;
    case 1: /* Converged successfully */
        fprintf(stderr, "Converged in %d iterations.           \n", Iter);
        break;
    case 2: /* Maximum iterations exceeded */
        fprintf(stderr, "Maximum number of iterations exceeded.\n");
        break;
    }
    return 1;
}


/*
 * Functions for options handling 
 */

/** 
* @brief Create a new tvregopt options object
* 
* This routine creates a new tvregopt options object and initializes it to
* default values.  It is the caller's responsibility to call TvRegFreeOpt
* to free the tvregopt object when done.
*/
tvregopt *TvRegNewOpt()
{
    tvregopt *Opt;
        
    if((Opt = (tvregopt *)Malloc(sizeof(struct tvregstruct))))
        *Opt = DefaultOpt;
    
    if(!(Opt->AlgString = (char *)Malloc(sizeof(char)*ALGSTRING_SIZE)))
    {
        Free(Opt);
        return NULL;
    }
    
    return Opt;
}


/** @brief Free tvregopt options object */
void TvRegFreeOpt(tvregopt *Opt)
{
    if(Opt)
    {
        Free(Opt->AlgString);        
        Free(Opt);
    }
}


/** @brief Specify fidelity weight lambda */
void TvRegSetLambda(tvregopt *Opt, num Lambda)
{
    if(Opt)
        Opt->Lambda = Lambda;
}


/**
 * @brief Specify spatially varying fidelity weight
 * @param Opt tvregopt options object
 * @param VaryingLambda pointer to Lambda array
 * @param LambdaWidth, LambdaHeight dimensions of the array
 * 
 * VaryingLambda should be a contiguous array of size LambdaWidth by 
 * LambdaHeight in row-major order of nonnegative values,
 *    VaryingLambda[x + Width*y] = fidelity weight at pixel (x,y).
 * Smaller VaryingLambda at a point implies stronger denoising, and a value
 * of zero specifies that the point should be inpainted.
 * 
 * If VaryingLambda = NULL, the constant Lambda value is used.
 * 
 * For inpainting, set VaryingLambda as 
 *    VaryingLambda[x + Width*y] = 0 if pixel (x,y) is unknown,
 *    VaryingLambda[x + Width*y] = C if pixel (x,y) is known,
 * where C is a positive constant.  Unknown pixels are inpainted (interpolated).
 * Known pixels are denoised (and deconvolved, if a kernel is also set).  To 
 * keep the known pixels (approximately) unchanged, set C to a large value.
 */
void TvRegSetVaryingLambda(tvregopt *Opt,
    const num *VaryingLambda, int LambdaWidth, int LambdaHeight)
{
    if(Opt)
    {
        Opt->VaryingLambda = VaryingLambda;
        Opt->LambdaWidth = LambdaWidth;
        Opt->LambdaHeight = LambdaHeight;
    }
}


/** 
 * @brief Specify kernel for a deconvolution problem
 * @param Opt tvregopt options object
 * @param Kernel pointer to convolution kernel
 * @param KernelWidth, KernelHeight dimensions of the kernel
 * 
 * Kernel should be a contiguous array of size KernelWidth by KernelHeight
 * in row-major order,
 *    Kernel[x + KernelWidth*y] = K(x,y).
 * If Kernel = NULL, then no deconvolution is performed.
 */
void TvRegSetKernel(tvregopt *Opt, 
    const num *Kernel, int KernelWidth, int KernelHeight)
{
    if(Opt)
    {
        Opt->Kernel = Kernel;
        Opt->KernelWidth = KernelWidth;
        Opt->KernelHeight = KernelHeight;
    }
}


/** @brief Specify convergence tolerance */
void TvRegSetTol(tvregopt *Opt, num Tol)
{
    if(Opt)
        Opt->Tol = Tol;
}


/** @brief Specify d = grad u constraint weight */
void TvRegSetGamma1(tvregopt *Opt, num Gamma1)
{
    if(Opt)
        Opt->Gamma1 = Gamma1;
}


/** @brief Specify z = Ku constraint weight */
void TvRegSetGamma2(tvregopt *Opt, num Gamma2)
{
    if(Opt)
        Opt->Gamma2 = Gamma2;
}


/** @brief Specify the maximum number of iterations */
void TvRegSetMaxIter(tvregopt *Opt, int MaxIter)
{
    if(Opt)
        Opt->MaxIter = MaxIter;
}


/**  
 * @brief Specify noise model
 * @param Opt tvregopt options object
 * @param NoiseModel string
 * 
 * NoiseModel should be a string specifying one of the following:
 * 
 *   - 'Gaussian' or 'L2'   (default) Additive white Gaussian noise (AWGN),
 *                          this is the noise model used in the traditional 
 *                          Rudin-Osher-Fatemi model;
 * 
 *   - 'Laplace' or 'L1'    Laplace noise, effective for salt & pepper noise;
 * 
 *   - 'Poisson'            Each pixel is an independent Poisson random
 *                          variable with mean equal to the exact value.
 */
int TvRegSetNoiseModel(tvregopt *Opt, const char *NoiseModel)
{
    if(!Opt)
        return 0;
    
    if(!NoiseModel || !strcmp(NoiseModel, "L2") || !strcmp(NoiseModel, "l2") 
        || !strcmp(NoiseModel, "Gaussian") || !strcmp(NoiseModel, "gaussian"))
        Opt->NoiseModel = NOISEMODEL_L2;
    else if(!strcmp(NoiseModel, "L1") || !strcmp(NoiseModel, "l1") 
        || !strcmp(NoiseModel, "Laplace") || !strcmp(NoiseModel, "laplace")
        || !strcmp(NoiseModel, "Laplacian") || !strcmp(NoiseModel, "laplacian"))
        Opt->NoiseModel = NOISEMODEL_L1;
    else if(!strcmp(NoiseModel, "Poisson") || !strcmp(NoiseModel, "poisson"))
        Opt->NoiseModel = NOISEMODEL_POISSON;
    else
        return 0;
    
    return 1;
}


/**
 * @brief Specify plotting function
 * @param Opt tvregopt options object
 * @param PlotFun plotting function
 * @param PlotParam void pointer for passing addition parameters
 * 
 * Specifying the plotting function gives control over how TvRestore displays 
 * information.  Setting PlotFun = NULL disables all normal display (error 
 * messages are still displayed).
 * 
 * An example PlotFun is
@code
    int ExamplePlotFun(int State, int Iter, num Delta,
        const num *u, int Width, int Height, int NumChannels, void *PlotParam)
    {
        switch(State)
        {
        case 0: 
            fprintf(stderr, " RUNNING   Iter=%4d, Delta=%7.4f\r", Iter, Delta);
            break;
        case 1: 
            fprintf(stderr, " CONVERGED Iter=%4d, Delta=%7.4f\n", Iter, Delta);
            break;
        case 2: 
            fprintf(stderr, " Maximum number of iterations exceeded!\n");
            break;
        }
        return 1;
    }
@endcode
 * The State argument is either 0, 1, or 2, and indicates TvRestore's status.
 * Iter is the number of Bregman iterations completed, Delta is the change in
 * the solution Delta = ||u^cur - u^prev||_2 / ||f||_2.  Argument u gives a 
 * pointer to the current solution, which can be used to plot an animated  
 * display of the solution progress.  PlotParam is a void pointer that can be
 * used to pass additional information to PlotFun if needed.
 */
void TvRegSetPlotFun(tvregopt *Opt, 
    int (*PlotFun)(int, int, num, const num*, int, int, int, void*),
    void *PlotParam)
{
    if(Opt)
    {
        Opt->PlotFun = PlotFun;
        Opt->PlotParam = PlotParam;
    }
}


/** @brief Debugging function that prints the current options */
void TvRegPrintOpt(const tvregopt *Opt)
{
    if(!Opt)
        Opt = &DefaultOpt;
    
    printf("lambda    : ");

    if(!Opt->VaryingLambda)
        printf("%g\n", Opt->Lambda);
    else
        printf("[%d x %d]\n", 
            Opt->LambdaWidth, Opt->LambdaHeight);
    
    printf("K         : ");

    if(!Opt->Kernel)
        printf("(identity)\n");
    else
        printf("[%d x %d]\n", Opt->KernelWidth, Opt->KernelHeight);

    printf("tol       : %g\n", (double)Opt->Tol);
    printf("max iter  : %d\n", Opt->MaxIter);
    printf("gamma1    : %g\n", (double)Opt->Gamma1);
    printf("gamma2    : %g\n", (double)Opt->Gamma2);
    printf("noise     : ");

    switch(Opt->NoiseModel)
    {
    case NOISEMODEL_L2:
        printf("L2\n");
        break;
    case NOISEMODEL_L1:
        printf("L1\n");
        break;
    case NOISEMODEL_POISSON:
        printf("Poisson\n");
        break;
    default:
        printf("(invalid)\n");
        break;
    }

    printf("plotting  : ");    

    if(Opt->PlotFun == TvRestoreSimplePlot)
        printf("default\n");
    else if(!Opt->PlotFun)
        printf("none\n");
    else
        printf("custom\n");
    
    printf("algorithm : %s\n", TvRegGetAlgorithm(Opt));
}


const char *TvRegGetAlgorithm(const tvregopt *Opt)
{
    static const char *DefaultAlgorithm = 
        (char *)"split Bregman (d = grad u) Gauss-Seidel u-solver";
    static const char *Invalid = (char *)"(invalid)";
    void (*ZSolveFun)(num*, num*, const num*, const num*, const num*,
        int, int, int, int, int, num, num);
    int UseZ, DeconvFlag, DctFlag;

    
    if(!Opt)
        return DefaultAlgorithm;
    
    if(!ChooseAlgorithm(&UseZ, &DeconvFlag, &DctFlag, &ZSolveFun, Opt))
        return Invalid;
    
    sprintf(Opt->AlgString, "split Bregman (%s) %s u-solver",
            (UseZ) ?
                "d = grad u, z = Ku" : 
                "d = grad u",
            (!DeconvFlag) ? 
                "Gauss-Seidel" :
                ((DctFlag) ? 
                    "DCT" :
                    "Fourier"));
    return Opt->AlgString;
}


/* This code is included if the file is being compiled as a MEX function */
#ifdef MATLAB_MEX_FILE

#ifdef NUM_SINGLE
#define NUM_CLASSNAME   "single"
#else
#define NUM_CLASSNAME   "double"
#endif

#define CAST_FULL       0x01
#define CAST_SPARSE     0x02
#define CAST_REAL       0x04
#define CAST_COMPLEX    0x08

#define IS_SCALAR(P) (mxIsNumeric(P) && mxGetNumberOfElements(P) == 1)

/**
 * @brief Cast the datatype of a numeric or logical mxArray
 * @param Deep is 0 if the result is the same array, 1 if it is a new array
 * @param Src is a pointer to the source mxArray object
 * @param Flags specifies whether array may be full, sparse, real, complex
 * @param Class is a null terminated string of the target datatype
 * @return pointer to casted mxArray, or NULL if array is not numeric
 * 
 * This routine should be used for example as
@code
    int XDeep;
    if(!(X = ArrayNumericCast(&XDeep, X, CAST_FULL | CAST_REAL, "single")))
        mexErrMsgTxt("X must be a numeric type.");
    
    ... use X ...
    
    if(XDeep)
        mxDestroyArray(X);
@endcode
 */
mxArray *ArrayNumericCast(int *Deep, const mxArray *Array, 
    int Flags, const char *Class)
{
    mxArray *Cast;
    
    
    *Deep = 0;
        
    if(!(Flags & (CAST_FULL | CAST_SPARSE))
        || !(Flags & (CAST_REAL | CAST_COMPLEX))
        || !mxIsNumeric(Array) && !mxIsLogical(Array))
            return NULL;

    /* Sparsity */
    if(mxIsSparse(Array))
    {
        if(!(Flags & CAST_SPARSE))  /* Convert sparse to full */
        {
            if(mexCallMATLAB(1, &Cast, 1, (mxArray **)&Array, "full"))
                mexErrMsgTxt("Cast failed.");
            
            *Deep = 1;
            Array = Cast;
        }
    }
    else
    {
        if(!(Flags & CAST_FULL))    /* Convert full to sparse */
        {
            if(mexCallMATLAB(1, &Cast, 1, (mxArray **)&Array, "sparse"))
                mexErrMsgTxt("Cast failed.");
            
            *Deep = 1;
            Array = Cast;
        }
    }
    
    /* Complexity */
    if(mxIsComplex(Array))
    {
        if(!(Flags & CAST_COMPLEX)) /* Explicit conversion is not necessary */
            mexWarnMsgTxt("Ignoring imaginary data.");
    }
    else
    {
        if(!(Flags & CAST_REAL))    /* Convert real to complex */
        {
            if(mexCallMATLAB(1, &Cast, 1, (mxArray **)&Array, "complex"))
                mexErrMsgTxt("Cast failed.");
            
            if(*Deep)
                mxDestroyArray((mxArray *)Array);
                
            *Deep = 1;
            Array = Cast;
        }
    }
            
    /* Datatype */
    if(strcmp(Class, mxGetClassName(Array)))  /* Convert datatype to Class */
    {   
        if(mexCallMATLAB(1, &Cast, 1, (mxArray **)&Array, Class))
            mexErrMsgTxt("Cast failed.");
        
        if(*Deep)
            mxDestroyArray((mxArray *)Array);
        
        *Deep = 1;
        Array = Cast;
    }
    
    return (mxArray *)Array;
}


int CallMatlabPlotFun(int State, int Iter, num Delta,
    const num *u, int Width, int Height, int NumChannels, void *Param)
{
    mxArray **PlotParam = (mxArray **)Param;
    mxArray *Lhs = NULL, *FevalRhs[5] = {PlotParam[0],
        NULL, NULL, NULL, PlotParam[1]};
    int Result;

    if(PlotParam[0] && PlotParam[1])
    {
        FevalRhs[1] = mxCreateDoubleScalar((double)State);
        FevalRhs[2] = mxCreateDoubleScalar((double)Iter);
        FevalRhs[3] = mxCreateDoubleScalar((double)Delta);

        if(mexCallMATLAB(0, &Lhs, 5, (mxArray **)FevalRhs, "feval"))
            mexErrMsgTxt("Error in plot function.");

        mexCallMATLAB(0, NULL, 0, NULL, "drawnow");

        if(Lhs && !mxIsEmpty(Lhs) && mxGetNumberOfElements(Lhs) == 1)
        {
            if(mxIsNumeric(Lhs))
                Result = (int)mxGetScalar(Lhs);
            else if(mxIsLogical(Lhs))
                Result = (int)(*mxGetLogicals(Lhs));
        }

        mxDestroyArray(Lhs);
        mxDestroyArray(FevalRhs[3]);
        mxDestroyArray(FevalRhs[2]);
        mxDestroyArray(FevalRhs[1]);
    }
    else
        mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
    
#ifdef MATLAB_CTRL_C
    /* Attempt to detect Ctrl+C */
    if(Result && utIsInterruptPending())
        Result = 0;
#endif    
    
    return Result;
}


void WarnUnknownOpt(const char *OptName)
{
    const char *IgnoringMessage = (char *)"Ignoring unknown option '";
    char *MessageBuf;
    int MessageLen;
    
    MessageLen = strlen(IgnoringMessage) + strlen(OptName) + 2;
    
    if(!(MessageBuf = (char *)mxMalloc(sizeof(char)*(MessageLen + 1))))
        mexWarnMsgTxt("Ignoring unknown option.");
    else
    {
        sprintf(MessageBuf, "%s%s'.", IgnoringMessage, OptName);
        mexWarnMsgTxt(MessageBuf);
        mxFree(MessageBuf);
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Input and ouput arguments */
    #define F_IN        prhs[0]
    #define LAMBDA_IN   prhs[1]
    #define OPT_IN      prhs[2]
    #define U0_IN       prhs[3]
    #define U_OUT       plhs[0]


    mxArray *Field, *F, *Lambda, *U0, *K, *PlotParam[2] = {NULL, NULL};
    const int *Size;
    const char *FieldName;
    char *Buf;
    tvregopt *Opt;
    int FDeep = 0, LambdaDeep = 0, U0Deep = 0, KDeep = 0;
    int k, BufLen, NumChannels, Verbose = 0;


    if(nrhs < 2)
        mexErrMsgTxt("At least two input arguments required.");
    else if(nrhs > 4)
        mexErrMsgTxt("Too many input arguments.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

    if(mxGetNumberOfDimensions(F_IN) > 3
        || !(F = ArrayNumericCast(&FDeep, F_IN,
        CAST_FULL | CAST_REAL, NUM_CLASSNAME)))
        mexErrMsgTxt("First argument must be a 2D or 3D numeric array.");
    
    Size = mxGetDimensions(F);
    NumChannels = (mxGetNumberOfDimensions(F) == 2) ? 1 : Size[2];
    
    if(!(Opt = TvRegNewOpt()))
        mexErrMsgTxt("Memory allocation failed.");
    
    if(IS_SCALAR(LAMBDA_IN))
        TvRegSetLambda(Opt, mxGetScalar(LAMBDA_IN));
    else if(mxGetNumberOfDimensions(LAMBDA_IN) != 2
        || !(Lambda = ArrayNumericCast(&LambdaDeep, LAMBDA_IN,
        CAST_FULL | CAST_REAL, NUM_CLASSNAME)))
        mexErrMsgTxt("Second argument must be a numeric scalar or 2D array.");
    else
    {
        if(mxGetM(Lambda) == Size[0] && mxGetN(Lambda) == Size[1])
            TvRegSetVaryingLambda(Opt, (num *)mxGetData(Lambda), 
                mxGetM(Lambda), mxGetN(Lambda));
        else
            mexErrMsgTxt("lambda must have the same number of rows and columns as the image.");
    }

    if(nrhs < 4)
        U_OUT = mxDuplicateArray(F);
    else if(mxGetNumberOfDimensions(U0_IN) != mxGetNumberOfDimensions(F)
        || mxGetNumberOfElements(U0_IN) != mxGetNumberOfElements(F_IN))
        mexErrMsgTxt("First and fourth arguments must have the same size.");
    else if(!(U0 = ArrayNumericCast(&U0Deep, U0_IN,
        CAST_FULL | CAST_REAL, NUM_CLASSNAME)))
        mexErrMsgTxt("Fourth argument must be a numeric array.");
    else
    {
        if(U0Deep)
        {
            U_OUT = U0;
            U0Deep = 0;
        }
        else
            U_OUT = mxDuplicateArray(U0);
    }
    
    if(nrhs >= 3 && !mxIsEmpty(OPT_IN))
    {
        if(!mxIsStruct(OPT_IN))
            mexErrMsgTxt("Third argument must be a struct.");
            
        /* Parse MATLAB options struct, convert to tvregopt */
        for(k = 0; k < mxGetNumberOfFields(OPT_IN); k++)
        {
            if(mxIsEmpty(Field = mxGetFieldByNumber(OPT_IN, 0, k)))
                continue;
            
            FieldName = mxGetFieldNameByNumber(OPT_IN, k);
            
            if(!strcmp(FieldName, "K"))
            {
                if(mxGetNumberOfDimensions(Field) != 2
                    || !(K = ArrayNumericCast(&KDeep, Field,
                    CAST_FULL | CAST_REAL, NUM_CLASSNAME)))
                    mexErrMsgTxt("Kernel must be a 2D numeric array.");
                else
                    TvRegSetKernel(Opt, (num *)mxGetData(K),
                        mxGetM(K), mxGetN(K));
            }
            else if(!strcmp(FieldName, "gamma1"))
            {
                if(!IS_SCALAR(Field))
                    mexErrMsgTxt("gamma1 must be a numeric scalar.");
                else
                    TvRegSetGamma1(Opt, (num)mxGetScalar(Field));
            }
            else if(!strcmp(FieldName, "gamma2"))
            {
                if(!IS_SCALAR(Field))
                    mexErrMsgTxt("gamma2 must be a numeric scalar.");
                else
                    TvRegSetGamma2(Opt, (num)mxGetScalar(Field));
            }
            else if(!strcmp(FieldName, "tol"))
            {
                if(!IS_SCALAR(Field))
                    mexErrMsgTxt("tol must be a numeric scalar.");
                else
                    TvRegSetTol(Opt, (num)mxGetScalar(Field));
            }
            else if(!strcmp(FieldName, "maxiter"))
            {
                if(!IS_SCALAR(Field))
                    mexErrMsgTxt("maxiter must be a numeric scalar.");
                else
                    TvRegSetMaxIter(Opt, (int)mxGetScalar(Field));
            }
            else if(!strcmp(FieldName, "noise"))
            {
                if(!mxIsChar(Field))
                    mexErrMsgTxt("noise must be a string.");
                else
                {
                    BufLen = mxGetNumberOfElements(Field)*sizeof(mxChar) + 1;
                    Buf = (char *)mxMalloc(BufLen);
                    mxGetString(Field, Buf, BufLen);
                    
                    if(!TvRegSetNoiseModel(Opt, Buf))
                        mexErrMsgTxt("Unknown noise model.");

                    mxFree(Buf);
                }
            }
            else if(!strcmp(FieldName, "plotfun"))
            {
                if(!mxIsChar(Field) && !mxIsFunctionHandle(Field))
                    mexErrMsgTxt("plotfun must be the name of a function or a function handle.");
                
                PlotParam[0] = Field;
                PlotParam[1] = U_OUT;
            }
            else if(!strcmp(FieldName, "verbose"))
            {
                if(IS_SCALAR(Field))
					Verbose = (int)mxGetScalar(Field);
            }
            else
                WarnUnknownOpt(FieldName);
        }
    }

#ifndef MATLAB_CTRL_C
    if(PlotParam[0] && PlotParam[1])
#endif
    TvRegSetPlotFun(Opt, CallMatlabPlotFun, (void *)PlotParam);

	if(Verbose > 0)
	{
		mexPrintf("f         : [%d x %d x %d]\n", Size[0], Size[1], NumChannels);
		TvRegPrintOpt(Opt);
		mexPrintf("implement : MEX, "
#ifdef NUM_SINGLE
		"single datatype"
#else
    	"double datatype"
#endif
#ifdef MATLAB_CTRL_C
		", libut Ctrl-C support\n");
#else
		"\n");
#endif
	}

    if(!TvRestore((num *)mxGetData(U_OUT), (num *)mxGetData(F),
        Size[0], Size[1], NumChannels, Opt))
        mexErrMsgTxt("Error in TvRestore.");

    if(KDeep)
        mxDestroyArray(K);
    if(U0Deep)
        mxDestroyArray(U0);
    if(LambdaDeep)
        mxDestroyArray(Lambda);
    if(FDeep)
        mxDestroyArray(F);
    TvRegFreeOpt(Opt);
}
#endif /* MATLAB_MEX_FILE */
