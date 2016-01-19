# MR-AMP
MR-AMP algorithms 

This toolbox is to implement the multi-resolution compressed sensing reconstruction via approximate message passing (MR-AMP)
proposed in the paper "Multi-resolution compressed sensing reconstruction via approximate message passing" (http://arxiv.org/abs/1508.02454).

Two different multi-resolution approaches are proposed, i.e., transform domain and spatial domain.

For transform domain approach, please find the demo file "DCT_WT_Demo.m" in the "demo" folder.
For spatial domain approach, please find the file "TVAMP_Demo.m" in the "demo" folder.
For the comparison between TVAL3 and AMP-TV-2D, please refer to "tval3_vs_tvamp.m" in the "demo" folder.

You can download TVAL3 toolbox here http://www.caam.rice.edu/~optimization/L1/TVAL3/ ,
and need to replace the TVAL3.m in the “solver” folder with TVAL3.m in our "demo" folder,
since we aim to find the optimal slack parameters for TVAL3.
