/**
* @file chanvese.c
* @brief Chan-Vese active contours without edges image segmentation
* @author Pascal Getreuer <getreuer@gmail.com>
*
* This file implements Chan-Vese active contours without edges two-phase
* image segmentation.  This file can be used either as part of a C program
* or compiled on its own as a MATLAB MEX function.
*
* To compile as a MEX function, run the command
*   mex chanvese.c
* in the MATLAB console.  The calling syntax is
*   phi = chanvese(f,phi0,Tol,MaxIter,mu,nu,lambda1,lambda2,dt,PlotFun)
* All arguments except f are optional.  Passing the empty matrix [] specifies
* the default for that parameter.  See chanvese.m for details.
*
* Pascal Getreuer 2007-2010
*
*
* License (BSD)
* 
* Copyright (c) 2010, Pascal Getreuer
* All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without 
* modification, are permitted provided that the following conditions are met:
* 
* - Redistributions of source code must retain the above copyright 
*   notice, this list of conditions and the following disclaimer.
* - Redistributions in binary form must reproduce the above copyright 
*   notice, this list of conditions and the following disclaimer in 
*   the documentation and/or other materials provided with the distribution.
*       
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "chanvese.h"

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

#define DIVIDE_EPS       ((num)1e-6)

#ifndef M_PI
/** @brief The constant pi */
#define M_PI        3.14159265358979323846264338327950288
#endif


/** @brief Options handling for ChanVese */
struct chanvesestruct
{
    num Tol;
    int MaxIter;
    num Mu;
    num Nu;
    num Lambda1;
    num Lambda2;
    num dt;
    int (*PlotFun)(int, int, num, const num*, const num*, const num*, 
        int, int, int, void*);
    void *PlotParam;
};

#ifdef __GNUC__
int ChanVeseSimplePlot(int State, int Iter, num Delta,
    const num *c1, const num *c2,
    __attribute__((unused)) const num *Phi,                          
    __attribute__((unused)) int Width, 
    __attribute__((unused)) int Height, 
    int NumChannels,
    __attribute__((unused)) void *Param);
#else
int ChanVeseSimplePlot(int State, int Iter, num Delta,
    const num *c1, const num *c2,
    const num *Phi,                          
    int Width, 
    int Height, 
    int NumChannels,
    void *Param);
#endif

/** @brief Default options struct */
static struct chanvesestruct DefaultChanVeseOpt =
        {(num)1e-4, 500, (num)0.25, 0, 1, 1, (num)0.5, 
        ChanVeseSimplePlot, NULL};
        

/**
* @brief Chan-Vese two-phase image segmentation
*
* @param Phi pointer to array to hold the resulting segmentation
* @param f the input image
* @param Width, Height, NumChannels the size of f
* @param Tol convergence tolerance
* @param MaxIter maximum number of iterations
* @param Mu length penalty
* @param Nu area penalty (positive penalizes area inside the curve)
* @param Lambda1 fit penalty inside the curve
* @param Lambda2 fit penalty outside the curve
* @param dt timestep
* @param PlotFun function for outputting intermediate results
*
* The input f can be a grayscale image or an image with any number of
* channels, i.e., three channels for a color image, or possibly many more in a
* hyperspectral image.  If f is a multichannel image, the segmentation is done
* using the Chan, Sandberg, Vese vector extension of the Chan-Vese model.
*
* The data for f should be stored as a contiguous block of data of 
* Width*Height*NumChannels elements, where the elements are ordered so that
*   f[x + Width*(y + Height*k)] = kth component of the pixel at (x,y)
*
* The array Phi is a contiguous array of size Width by Height with the same
* order as f.  Phi is a "level set function" of the segmentation, meaning the
* segmentation is indicated by its sign:
*    Phi[x + Width*y] >= 0 means (x,y) is inside the segmentation curve,
*    Phi[x + Width*y] <  0 means (x,y) is outside.
* Before calling this routine, Phi should be initialized either by calling 
* InitPhi or by setting it to a level set function of an initial guess of the
* segmentation.  After this routine, the final segmentation is obtained from
* the sign of Phi.
*  
* The routine runs at most MaxIter number of iterations and stops when the
* change between successive iterations is less than Tol.  Set Tol=0 to force
* the routine to run exactly MaxIter iterations.
*/
int ChanVese(num *Phi, const num *f, 
    int Width, int Height, int NumChannels, const chanveseopt *Opt)
{
    int (*PlotFun)(int, int, num, const num*, const num*, const num*, 
        int, int, int, void*);
    const int ChannelStep = Width*Height;
    const int NumEl = Width*Height*NumChannels;
    const num *fPtr, *fPtr2;
    num *PhiPtr, *c1 = 0, *c2 = 0;
    num c1Scalar, c2Scalar, Mu, Nu, Lambda1, Lambda2, dt;
    num PhiLast, Delta, PhiX, PhiY, IDivU, IDivD, IDivL, IDivR;
    num Temp1, Temp2, Dist1, Dist2;
    int Iter, i, j, Channel, PhiTol, PhiDiff, MaxIter, Success = 2;
    
    
    if(!Phi || !f || Width <= 0 || Height <= 0 || NumChannels <= 0)
        return 0;
    
    if(!Opt)
        Opt = &DefaultChanVeseOpt;
    
    Mu = Opt->Mu;
    Nu = Opt->Nu;
    Lambda1 = Opt->Lambda1;
    Lambda2 = Opt->Lambda2;
    dt = Opt->dt;
    MaxIter = Opt->MaxIter;
    PlotFun = Opt->PlotFun;
    
    PhiTol = (int)(Opt->Tol*NumEl + 0.5);
    PhiDiff = (PhiTol > 0) ? PhiTol*1000 : 1000;
    
    if(NumChannels > 1)
    {
        if(!(c1 = Malloc(sizeof(num)*NumChannels)) 
            || !(c2 = Malloc(sizeof(num)*NumChannels)))
        {
            fprintf(stderr, "Out of memory\n.");
            return 0;
        }
    }
    else
    {
        c1 = &c1Scalar;
        c2 = &c2Scalar;
    }
    
    RegionAverages(c1, c2, Phi, f, Width, Height, NumChannels);
    
    if(PlotFun)
        if(!PlotFun(0, 0, ((num)PhiDiff)/NumEl, c1, c2, Phi,
                Width, Height, NumChannels, Opt->PlotParam))
            goto Done;
    
    for(Iter = 1; Iter <= MaxIter; Iter++)
    {
        PhiPtr = Phi + 1 + Width;
        fPtr = f + 1 + Width;
        PhiDiff = 0;
        
        for(j = 1; j < Height-1; j++, PhiPtr += 2, fPtr += 2)
        {
            for(i = 1; i < Width-1; i++, PhiPtr++, fPtr++)
            {
                Delta = dt/(1 + PhiPtr[0]*PhiPtr[0]);
                PhiX = PhiPtr[1] - PhiPtr[0];
                PhiY = ( PhiPtr[Width] - PhiPtr[-Width])/2;
                IDivR = (num)(1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY));
                PhiX = PhiPtr[0] - PhiPtr[-1];
                IDivL = (num)(1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY));
                PhiX = (PhiPtr[1] - PhiPtr[-1])/2;
                PhiY =  PhiPtr[Width] - PhiPtr[0];
                IDivD = (num)(1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY));
                PhiY = PhiPtr[0] - PhiPtr[-Width];
                IDivU = (num)(1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY));
                
                if(NumChannels == 1)
                {
                    Dist1 = fPtr[0] - c1Scalar;
                    Dist2 = fPtr[0] - c2Scalar;
                    Dist1 *= Dist1;
                    Dist2 *= Dist2;
                }
                else    
                {
                    Dist1 = Dist2 = 0.0;
                    
                    for(Channel = 0, fPtr2 = fPtr; 
                        Channel < NumChannels; Channel++, fPtr2 += ChannelStep)
                    {
                        Temp1 = fPtr2[0] - c1[Channel];
                        Temp2 = fPtr2[0] - c2[Channel];
                        Dist1 += Temp1*Temp1;
                        Dist2 += Temp2*Temp2;
                    }
                }
                
                /* Semi-implicit update of PHI_CENTER */
                PhiLast = PhiPtr[0];
                PhiPtr[0] = (PhiPtr[0] + Delta*(
                        Mu*(PhiPtr[1]*IDivR + PhiPtr[-1]*IDivL
                            + PhiPtr[Width]*IDivD + PhiPtr[-Width]*IDivU)
                        - Nu - Lambda1*Dist1 + Lambda2*Dist2) ) /
                    (1 + Delta*Mu*(IDivR + IDivL + IDivD + IDivU));

                if(PhiLast*PhiPtr[0] < 0)
                    PhiDiff++;
            }
        }
        
        /* Neumann boundary conditions */
        for(i = 1; i < Width-1; i++)
            Phi[i] = Phi[i+Width];
        for(i = 1; i < Width-1; i++)
            Phi[i+Width*(Height-1)] = Phi[i+Width*(Height-2)];
        for(j = 0; j < Height; j++)
            Phi[j*Width] = Phi[1+j*Width];
        for(j = 0; j < Height; j++)
            Phi[(Width-1) + j*Width] = Phi[(Width-2)+j*Width];
        
        RegionAverages(c1, c2, Phi, f, Width, Height, NumChannels);
        
        if(Iter >= 2 && PhiDiff <= PhiTol)
            break;
        
        if(PlotFun)
            if(!PlotFun(0, Iter, ((num)PhiDiff)/NumEl, c1, c2, Phi,
                    Width, Height, NumChannels, Opt->PlotParam))
                goto Done;
    }

    Success = (Iter <= MaxIter) ? 1:2;

    if(PlotFun)
        PlotFun(Success, (Iter <= MaxIter) ? Iter:MaxIter, 
                ((num)PhiDiff)/NumEl, c1, c2, Phi, 
                Width, Height, NumChannels, Opt->PlotParam);
    
Done:        
    if(NumChannels > 1)
    {
        Free(c2);
        Free(c1);
    }
    
    return Success;    
}


/** @brief Default initialization for Phi */
void ChanVeseInitPhi(num *Phi, int Width, int Height)
{
    int i, j;
    
    for(j = 0; j < Height; j++)
        for(i = 0; i < Width; i++)
            *(Phi++) = (num)(sin(i*M_PI/5.0)*sin(j*M_PI/5.0));
}


/** @brief Compute averages inside and outside of the segmentation contour */
void RegionAverages(num *c1, num *c2, const num *Phi, const num *f,
    int Width, int Height, int NumChannels)
{
    const int NumPixels = Width * Height;
    num Sum1 = 0, Sum2 = 0;    
    int n, Channel, Count1 = 0, Count2 = 0;
    
    for(Channel = 0; Channel < NumChannels; Channel++, f += NumPixels)
    {
        for(n = 0; n < NumPixels; n++)
        {
            if(Phi[n] >= 0)
            {
                Count1++;
                Sum1 += f[n];
            }
            else
            {
                Count2++;
                Sum2 += f[n];
            }
        }  
        
        c1[Channel] = (Count1) ? (Sum1/Count1) : 0;
        c2[Channel] = (Count2) ? (Sum2/Count2) : 0;
    }
}



/* If GNU C language extensions are available, apply the "unused" attribute
   to avoid warnings.  TvRestoreSimplePlot is a plotting callback function
   for TvRestore, so the unused arguments are indeed required. */
#ifdef __GNUC__
int ChanVeseSimplePlot(int State, int Iter, num Delta,
    const num *c1, const num *c2,
    __attribute__((unused)) const num *Phi,                          
    __attribute__((unused)) int Width, 
    __attribute__((unused)) int Height, 
    int NumChannels,
    __attribute__((unused)) void *Param)
#else
int ChanVeseSimplePlot(int State, int Iter, num Delta,
    const num *c1, const num *c2,
    const num *Phi,                          
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
        if(NumChannels == 1)
            fprintf(stderr, "   Iteration %4d     Delta %7.4f     c1 = %6.4f     c2 = %6.4f\r", 
                Iter, Delta, *c1, *c2);
        else
            fprintf(stderr, "   Iteration %4d     Delta %7.4f\r", Iter, Delta);
        break;
    case 1: /* Converged successfully */
        fprintf(stderr, "Converged in %d iterations.                                            \n", 
            Iter);
        break;
    case 2: /* Maximum iterations exceeded */
        fprintf(stderr, "Maximum number of iterations exceeded.                                 \n");
        break;
    }
    return 1;
}


/*
 * Functions for options handling 
 */

/** 
* @brief Create a new chanveseopt options object
* 
* This routine creates a new chanveseopt options object and initializes it to
* default values.  It is the caller's responsibility to call ChanVeseFreeOpt
* to free the chanveseopt object when done.
*/
chanveseopt *ChanVeseNewOpt()
{
    chanveseopt *Opt;
        
    if((Opt = (chanveseopt *)Malloc(sizeof(struct chanvesestruct))))
        *Opt = DefaultChanVeseOpt;
    
    return Opt;
}


/** @brief Free chanveseopt options object */
void ChanVeseFreeOpt(chanveseopt *Opt)
{
    if(Opt)
        Free(Opt);
}


/** @brief Specify mu, the edge length penalty */
void ChanVeseSetMu(chanveseopt *Opt, num Mu)
{
    if(Opt)
        Opt->Mu = Mu;
}


/** @brief Specify nu, the area penalty (may be positive or negative) */
void ChanVeseSetNu(chanveseopt *Opt, num Nu)
{
    if(Opt)
        Opt->Nu = Nu;
}


/** @brief Specify lambda1, the fit weight inside the curve */
void ChanVeseSetLambda1(chanveseopt *Opt, num Lambda1)
{
    if(Opt)
        Opt->Lambda1 = Lambda1;
}


/** @brief Specify lambda2, the fit weight outside the curve */
void ChanVeseSetLambda2(chanveseopt *Opt, num Lambda2)
{
    if(Opt)
        Opt->Lambda2 = Lambda2;
}


/** @brief Specify the convergence tolerance */
void ChanVeseSetTol(chanveseopt *Opt, num Tol)
{
    if(Opt)
        Opt->Tol = Tol;
}


/** @brief Specify the timestep */
void ChanVeseSetDt(chanveseopt *Opt, num dt)
{
    if(Opt)
        Opt->dt = dt;
}


/** @brief Specify the maximum number of iterations */
void ChanVeseSetMaxIter(chanveseopt *Opt, int MaxIter)
{
    if(Opt)
        Opt->MaxIter = MaxIter;
}


/**
 * @brief Specify plotting function
 * @param Opt chanveseopt options object
 * @param PlotFun plotting function
 * @param PlotParam void pointer for passing addition parameters
 * 
 * Specifying the plotting function gives control over how ChanVese displays 
 * information.  Setting PlotFun = NULL disables all normal display (error 
 * messages are still displayed).
 * 
 * An example PlotFun is
@code
    int ExamplePlotFun(int State, int Iter, num Delta, 
        const num *c1, const num *c2,   const num *Phi,                          
        int Width, int Height, int NumChannels, void *Param)
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
 * The State argument is either 0, 1, or 2, and indicates ChanVese's status.
 * Iter is the number of Bregman iterations completed, Delta is the change in
 * the solution Delta = ||u^cur - u^prev||_2 / ||f||_2.  Argument u gives a 
 * pointer to the current solution, which can be used to plot an animated  
 * display of the solution progress.  PlotParam is a void pointer that can be
 * used to pass additional information to PlotFun if needed.
 */
void ChanVeseSetPlotFun(chanveseopt *Opt, 
    int (*PlotFun)(int, int, num, const num*, const num*, const num*, 
        int, int, int, void*), void *PlotParam)
{
    if(Opt)
    {
        Opt->PlotFun = PlotFun;
        Opt->PlotParam = PlotParam;
    }
}


void ChanVesePrintOpt(const chanveseopt *Opt)
{
    if(!Opt)
        Opt = &DefaultChanVeseOpt;
    
    printf("tol       : %g\n", Opt->Tol);
    printf("max iter  : %d\n", Opt->MaxIter);
    printf("mu        : %g\n", Opt->Mu);
    printf("nu        : %g\n", Opt->Nu);
    printf("lambda1   : %g\n", Opt->Lambda1);
    printf("lambda2   : %g\n", Opt->Lambda2);
    printf("dt        : %g\n", Opt->dt);
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
    const num *c1, const num *c2, const num *Phi, 
    int Width, int Height, int NumChannels, void *Param)
{
    mxArray **PlotParam = (mxArray **)Param;
    mxArray *Lhs = NULL, *FevalRhs[5] = {PlotParam[0],
        NULL, NULL, NULL, PlotParam[1]};
    int Result = 1;
    
    if(PlotParam[0] && PlotParam[1])
    {
        FevalRhs[1] = mxCreateDoubleScalar((double)State);
        FevalRhs[2] = mxCreateDoubleScalar((double)Iter);
        FevalRhs[3] = mxCreateDoubleScalar((double)Delta);        
        
        if(mexCallMATLAB(0, &Lhs, 5, FevalRhs, "feval"))
            mexErrMsgTxt("Error evaluating plot function.");
        
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


void ParseParam(chanveseopt *Opt, int *Verbose, 
	const char *FieldName, const mxArray *Field,
    mxArray *Phi_out, const mxArray *PlotParam[])
{
    if(!strcmp(FieldName, "mu"))
    {
        if(!IS_SCALAR(Field))
            mexErrMsgTxt("mu must be a numeric scalar.");
        else
            ChanVeseSetMu(Opt, (num)mxGetScalar(Field));
    }
    else if(!strcmp(FieldName, "nu"))
    {
        if(!IS_SCALAR(Field))
            mexErrMsgTxt("nu must be a numeric scalar.");
        else
            ChanVeseSetNu(Opt, (num)mxGetScalar(Field));
    }
    else if(!strcmp(FieldName, "lambda1"))
    {
        if(!IS_SCALAR(Field))
            mexErrMsgTxt("lambda1 must be a numeric scalar.");
        else
            ChanVeseSetLambda1(Opt, (num)mxGetScalar(Field));
    }
    else if(!strcmp(FieldName, "lambda2"))
    {
        if(!IS_SCALAR(Field))
            mexErrMsgTxt("lambda2 must be a numeric scalar.");
        else
            ChanVeseSetLambda2(Opt, (num)mxGetScalar(Field));
    }
    else if(!strcmp(FieldName, "tol"))
    {
        if(!IS_SCALAR(Field))
            mexErrMsgTxt("tol must be a numeric scalar.");
        else
            ChanVeseSetTol(Opt, (num)mxGetScalar(Field));
    }
    else if(!strcmp(FieldName, "maxiter"))
    {
        if(!IS_SCALAR(Field))
            mexErrMsgTxt("maxiter must be a numeric scalar.");
        else
            ChanVeseSetMaxIter(Opt, (int)mxGetScalar(Field));
    }
    else if(!strcmp(FieldName, "dt"))
    {
        if(!IS_SCALAR(Field))
            mexErrMsgTxt("dt must be a numeric scalar.");
        else
            ChanVeseSetDt(Opt, (num)mxGetScalar(Field));
    }
    else if(!strcmp(FieldName, "plotfun"))
    {
        if(!mxIsChar(Field) && !mxIsFunctionHandle(Field))
            mexErrMsgTxt("plotfun must be the name of a function or a function handle.");
        
        PlotParam[0] = Field;
        PlotParam[1] = Phi_out;
    }
    else if(!strcmp(FieldName, "verbose"))
    {
    	if(IS_SCALAR(Field))
			*Verbose = (int)mxGetScalar(Field);
    }
    else
        WarnUnknownOpt(FieldName);
}


/** @brief MEX gateway */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    #define	F_IN	     prhs[0]
    #define PHI0_IN      prhs[1]
    #define PHI_OUT      plhs[0]
    
    /* Calling syntax with options struct */    
    #define OPT_IN       prhs[2]
        
    /* Second calling syntax */
    #define TOL_IN       prhs[2]
    #define MAXITER_IN   prhs[3]
    #define MU_IN        prhs[4]
    #define NU_IN        prhs[5]
    #define LAMBDA1_IN   prhs[6]    
    #define LAMBDA2_IN   prhs[7]    
    #define DT_IN        prhs[8]
    #define PLOTFUN_IN   prhs[9]
    
    const mxArray *Phi0_in = NULL, *PlotParam[2] = {NULL, NULL};
    mxArray *Field, *F, *Phi0;
    chanveseopt *Opt;
    const int *Size;
    char *Buf;
    int k, BufLen, FDeep = 0, Phi0Deep = 0, MaxIter, NumChannels, Verbose = 0;
    
    /* Parse the input arguments */
    if(nrhs < 1)
        mexErrMsgTxt("At least one input argument required.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");
    
    if(mxGetNumberOfDimensions(F_IN) > 3
        || !(F = ArrayNumericCast(&FDeep, F_IN,
        CAST_FULL | CAST_REAL, NUM_CLASSNAME)))
        mexErrMsgTxt("First argument must be a 2D or 3D numeric array.");
    
    Size = mxGetDimensions(F);
    NumChannels = (mxGetNumberOfDimensions(F) == 2) ? 1 : Size[2];   
    
    if(!(Opt = ChanVeseNewOpt()))
        mexErrMsgTxt("Memory allocation failed.");
    
    /* Initialize Phi */
    if(nrhs < 2 || mxIsEmpty(PHI0_IN))
    {
        /* Create output matrix Phi */
#ifdef NUM_SINGLE
        PHI_OUT = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxREAL);
#else
        PHI_OUT = mxCreateDoubleMatrix(0, 0, mxREAL); 
#endif
        mxSetDimensions(PHI_OUT, (const mwSize *)Size, 2);
        mxSetData(PHI_OUT, mxMalloc(Size[0]*Size[1]*sizeof(num)));
        ChanVeseInitPhi((num *)mxGetData(PHI_OUT), Size[0], Size[1]);
    }
    else if(mxGetNumberOfDimensions(PHI0_IN) != 2
        || !(Phi0 = ArrayNumericCast(&Phi0Deep, PHI0_IN,
        CAST_FULL | CAST_REAL, NUM_CLASSNAME)))
        mexErrMsgTxt("phi0 must be a 2D numeric array.");
    else if(mxGetM(Phi0) != Size[0] || mxGetN(Phi0) != Size[1])
        mexErrMsgTxt("phi0 must have the same number of rows and columns as the image.");
    else
    {
        if(Phi0Deep)
        {
            PHI_OUT = Phi0;
            Phi0Deep = 0;
        }
        else
            PHI_OUT = mxDuplicateArray(Phi0);
    }
    
    if(nrhs >= 3 && mxIsStruct(OPT_IN))
    {
        if(nrhs > 3)
            mexErrMsgTxt("Too many input arguments.");
        
        /* Parse MATLAB options struct, convert to chanveseopt */
        for(k = 0; k < mxGetNumberOfFields(OPT_IN); k++)
        {
            if(mxIsEmpty(Field = mxGetFieldByNumber(OPT_IN, 0, k)))
                continue;
            
            ParseParam(Opt, &Verbose, mxGetFieldNameByNumber(OPT_IN, k),
                Field, PHI_OUT, PlotParam);
        }
    }
    else
    {
        if(nrhs > 10)
            mexErrMsgTxt("Too many input arguments.");
        
        if(nrhs >= 3)
            ParseParam(Opt, &Verbose, "tol", TOL_IN, PHI_OUT, PlotParam);
        if(nrhs >= 4)
            ParseParam(Opt, &Verbose, "maxiter", MAXITER_IN, PHI_OUT, PlotParam);
        if(nrhs >= 5)
            ParseParam(Opt, &Verbose, "mu", MU_IN, PHI_OUT, PlotParam);
        if(nrhs >= 6)
            ParseParam(Opt, &Verbose, "nu", NU_IN, PHI_OUT, PlotParam);
        if(nrhs >= 7)
            ParseParam(Opt, &Verbose, "lambda1", LAMBDA1_IN, PHI_OUT, PlotParam);
        if(nrhs >= 8)
            ParseParam(Opt, &Verbose, "lambda2", LAMBDA2_IN, PHI_OUT, PlotParam);
        if(nrhs >= 9)
            ParseParam(Opt, &Verbose, "dt", DT_IN, PHI_OUT, PlotParam);
        if(nrhs >= 10)
            ParseParam(Opt, &Verbose, "plotfun", PLOTFUN_IN, PHI_OUT, PlotParam);
    }
    
#ifndef MATLAB_CTRL_C
    if(PlotParam[0] && PlotParam[1])
#endif    
    ChanVeseSetPlotFun(Opt, CallMatlabPlotFun, (void *)PlotParam);

    if(Verbose > 0)
    {
		mexPrintf("f         : [%d x %d x %d]\n", Size[0], Size[1], NumChannels);
    	ChanVesePrintOpt(Opt);
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

    /* Call the main solver routine */
    if(!ChanVese((num *)mxGetData(PHI_OUT), (const num *)mxGetData(F), 
        Size[0], Size[1], NumChannels, Opt))
        mexErrMsgTxt("Error in ChanVese.");
    
    if(Phi0Deep)
        mxDestroyArray(Phi0);
    if(FDeep)
        mxDestroyArray(F);
    ChanVeseFreeOpt(Opt);
    return;
}
#endif  /* MATLAB_MEX_FILE */
