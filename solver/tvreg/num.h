/**
 * @file num.h
 * @brief num typedef
 * @author Pascal Getreuer <getreuer@gmail.com>
 * 
 * This file defines type "num", which by default is a typedef for double.
 * If NUM_SINGLE is defined, then num is a typedef for float.
 * 
 * The file also defines a name mangling macro FFT(S) for using num with
 * the FFTW library.  The macro is defined such that FFT(functionname) 
 * expands to 
 *    fftwf_functionname  if num is single,
 * or 
 *    fftw_functionname   if num is double.
 */
#ifndef _NUM_H_
#define _NUM_H_

/** 
 * @brief  Token concatenation macro 
 * 
 * This extra level of indirection is needed so that macros expand before 
 * token pasting.
 */
#define _NUM_H_CONCAT(A,B)    A ## B


#ifdef NUM_SINGLE

/* Definitions for single precision */
typedef float num;
#define FFT(S)      _NUM_H_CONCAT(fftwf_,S)

#else

/* Definitions for double precision */
typedef double num;
#define FFT(S)      _NUM_H_CONCAT(fftw_,S)

#endif

#endif /* _NUM_H_ */
