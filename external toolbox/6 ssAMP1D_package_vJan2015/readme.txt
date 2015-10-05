README for ssAMP-1D solver written by Jaewook Kang in Gwangju institude of science and technology (GIST).
(jwkkang@gist.ac.kr, https://sites.google.com/site/jwkang10/ ) 
=============================================================================================================
    
  This file contains methods for performing
  "compressed sensing recontructions of one dimensional piecewise-constant signal via ssAMP-1D solver."
  
  한글명: 구분적 상수인 1차원 신호 복구를 위한 근사 메시지 전달 (AMP) 알고리즘 관련 알고리즘 비교 프로그램
  
  *Related publication: 
   Jaewook Kang, HyoYoung Jung, Heung-No Lee, and Kiseon Kim,  
   - "One-dimensional  Piecewise-Constant Signal Recovery via  Spike-and-Slab Approximate Message-Passing,"  
    proc. of the 48th Asilomar Conference (Pacific Grove, CA), Nov. 2014.
   - Jaewook Kang, Hyoyoung Jung, Heung-No Lee, Kiseon Kim,  "Spike-and-Slab Approximate Message-Passing  Recovery for 1-D Piecewise-Constant Signal,"  
    submitted to IEEE Transactions on Signal Processing Feb. 2015.
  =========================================================================================================
  
  The file "demo_1DFDrev.m" in this folder,  provides an exemplary signal 
  recovery comparison of ssAMP-1D algorithm to the recent solvers:
  
  [1] EFLA : the corresonding paper - J. Liu, L. Yuan, and J. Ye. 
%    “An efficient algorithm for a class of fused lasso problems,” 
%     proc of ACM SIGKDD Conference on Knowledge Discovery and Data Mining, 
%     2010. (obtained from http://www.public.asu.edu/~jye02/Software/SLEP/)

% [2] TV-AMP: the corresonding paper - D. L. Donoho, I. Johnstone, and A. Montanari, 
%     “Accurate prediction of phase transitions in compressed sensing via 
%      a connection to minimax denoising, ” IEEE Trans. Inform. Theory, 
%      vol. 59, no. 6, pp. 3396-3433, June 2013.

% [3] Chambolle-Pock TV: the corresonding paper - A. Chambolle, T. Pock, 
%    “A first-order primal-dual algorithm for convex
%      problems with applications to imaging,” J. Math. Imag. Vis., vol. 40, pp.
%      120-145, May 2011.

% [4] GrAMPA: the corresonding paper - M. Borgerding and P. Schiniter, 
%     “Generalized approximate message passing for the cosparse analysis
%      model,” avabilable at ArXiv:1312.3968v1 [cs.IT], Dec. 2013 
      (obtained from http://www2.ece.ohio-state.edu/~schniter/GrAMPA/download.html )

   The comparison provides 
   
   1) Normalized mean-squared-error of the signal recovery, and 
   2) CPU runtime spending for the recovery.

All the copyrights for this software are at Gwangju institude of science and technology 
(Registration No. C-2014-032048 by KOREA COPYRIGHT COMMISSION (http://www.copyright.or.kr/main.do) )

Final update 2015 Jan. 
**Special thank to prof. Schniter Phil (Ohio state University) for help to construct the experiement configuration
========================================================================================================================

The total size of the file is 1.08 MB 

Player Information:
===================
"	MATLAB Version: 8.2.0.701 (R2013b)
"	Operating System: Microsoft Windows 7 Version 6.1 (Build 7601: Service Pack 1)
"	Java Version: Java 1.7.0_11-b21 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

Contact Information:
===================
"	Lab. phone: +82-62-715-2264
"	E-mail        :  jwkkang@gist.ac.kr, jwkang10@gmail.com
      

      
