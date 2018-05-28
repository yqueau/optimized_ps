# optimized_ps

This Matlab demo code illustrates the use of state-of-the-art nonconvex optimization techniques for refining a baseline photometric stereo. The approach is described in the papers: 

[1] "Optimisation of classic photometric stereo by non-convex variational minimisation", Georg Radow, Laurent Hoeltgen, Yvain Quéau and Michael Breuss. To appear in Journal of Mathematical Imaging and Vision (2018)

[2] "Optimised photometric stereo via non-convex variational minimisation", Laurent Hoeltgen, Yvain Quéau, Michael Breuss and Georg Radow. In: Proceedings of BMVC (2016)

## Description

The demo code consists in two parts:

1. Baseline photometric stereo. Normals times albedo are estimated using pseudo-inverse (i.e., least-squares estimation). An estimate of the depth gradient is then obtained, which is eventually integrated in a least-squares sense into a depth map. 

2. Optimised photometric stereo. Depth and albedo are then refined in a least-squares fashion as well, but considering the original non-convex variational model. Optimisation is carried out using the state-of-the-art iPiano algorithm. 

Both approaches are based on least-squares and should thus provide similar resuts. However, the demo shows that the baseline approach is biased due to the non-integrability of the estimated normal field. The refined result gives a "better" explanation of the input images in terms of the quadratic reprojection error which is the metric used in both methods. 
