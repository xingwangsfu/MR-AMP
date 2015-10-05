/**
 * @file tvreg.h
 * @brief TV-regularized image restoration
 * @author Pascal Getreuer <getreuer@gmail.com>
 */
#ifndef _TVREG_H_
#define _TVREG_H_

#include "num.h"


typedef struct tvregstruct tvregopt;


tvregopt *TvRegNewOpt();
void TvRegFreeOpt(tvregopt *Opt);
void TvRegSetLambda(tvregopt *Opt, num Lambda);
void TvRegSetVaryingLambda(tvregopt *Opt,
    const num *VaryingLambda, int LambdaWidth, int LambdaHeight);
void TvRegSetKernel(tvregopt *Opt, 
    const num *Kernel, int KernelWidth, int KernelHeight);
void TvRegSetTol(tvregopt *Opt, num Tol);
void TvRegSetGamma1(tvregopt *Opt, num Gamma1);
void TvRegSetGamma2(tvregopt *Opt, num Gamma2);
void TvRegSetMaxIter(tvregopt *Opt, int MaxIter);
int TvRegSetNoiseModel(tvregopt *Opt, const char *NoiseModel);
void TvRegSetPlotFun(tvregopt *Opt, 
    int (*PlotFun)(int, int, num, const num*, int, int, int, void*),
    void *PlotParam);

void TvRegPrintOpt(const tvregopt *Opt);
const char *TvRegGetAlgorithm(const tvregopt *Opt);

int TvRestore(num *u, const num *f, int Width, int Height, int NumChannels,
     const tvregopt *Opt);

#endif /* _TVREG_H_ */
