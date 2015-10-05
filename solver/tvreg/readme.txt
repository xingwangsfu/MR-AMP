

   tvreg v2: Variational Imaging Methods for Denoising, Deconvolution,
          Inpainting, and Segmentation

       Pascal Getreuer, December 2010


The tvreg package applies total variation (TV) regularization to
perform image denoising, deconvolution, and inpainting. Three different
noise models are supported: Gaussian (L2), Laplace (L1), and Poisson.
The implementation solves the general TV restoration problem

                  /
   Min  TV(u)  +  | lambda F(K*u,f) dx,
    u             /

to perform denoising, deconvolution, and inpainting as special cases.
It is efficiently solved using the recent split Bregman method.  Also
included is an efficient implementation of Chan-Vese two-phase 
segmentation.  All functions support grayscale, color, and arbitrary
multichannel images.

Please see the included documentation file tvreg.pdf for details.


   ---------------------
    Get Started Quickly
   ---------------------

1. Install the FFTW library (http://www.fftw.org).  Windows users can
   download precompiled .dll files from 
       http://www.fftw.org/install/windows.html.

2. Compile the programs with GCC using "make -f makefile.gcc" or
   Microsoft Visual C++ using "nmake -f makefile.vc".  See section 7 of
   the documentation for help.

3. Try the demos

   tvdenoise_demo     Total variation denoising demo
   tvdeconv_demo      Total variation deconvolution demo
   tvinpaint_demo     Total variation inpainting demo
   chanvese_demo      Chan-Vese segmentation demo


   ------------------------------- 
    Get Started Quickly in MATLAB
   -------------------------------

Compiling is not required to use tvreg in MATLAB.  Try the demos

   tvdenoise_demo     Total variation denoising demo
   tvdeconv_demo      Total variation deconvolution demo
   tvinpaint_demo     Total variation inpainting demo
   chanvese_demo      Chan-Vese segmentation demo

For improved performance, run the included script "complex_mex.m" to
compile the main computation routines as MEX functions.  This requires
that FFTW is installed, please see section 7.3 of the documentation.

This material is based upon work supported by the National Science 
Foundation under Award No. DMS-1004694.  Any opinions, findings, and 
conclusions or recommendations expressed in this material are those of 
the author(s) and do not necessarily reflect the views of the National
Science Foundation.
