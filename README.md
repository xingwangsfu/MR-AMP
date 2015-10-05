# MR-AMP
MR-AMP algorithms 

This toolbox is to implement the multi-resolution compressed sensing reconstruction via approximate message passing (MR-AMP)
proposed in the paper "Multi-resolution compressed sensing reconstruction via approximate message passing" (http://arxiv.org/abs/1508.02454).

Two different multi-resolution approaches are proposed, i.e., transform domain and spatial domain.

For transform domain approach, please find the demo file "DCT_WT_Demo.m" in the "demo" folder.
For spatial domain approach, please find the file "TVAMP_Demo.m" in the "demo" folder.

Please note that to run DCT_WT_Demo.m, you need to download the SI_AMP toolbox first, since the denoiser used in DCT_WT_Demo.m is 
the soft-thresholding function implemented in SI_AMP toolbox. You can download the SI_AMP toolbox here (https://github.com/zhanzhangyu/SI-AMP).
After downloading, please add SI_AMP to the path.
