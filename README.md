# Spherical Density-Equalizing Map

<img src = "https://github.com/garyptchoi/spherical-density-equalizing-map/blob/main/cover.jpg" height="400" />

* Spherical density-equalizing map (SDEM): Compute a spherical density-equalizing map of a genus-0 closed surface using the method in [1].

* Landmark-aligned spherical density-equalizing map (LSDEM): Compute a landmark-aligned spherical density-equalizing map of a genus-0 closed surface using the method in [1].

Any comments and suggestions are welcome. 

If you use this code in your work, please cite the following paper:

[1] Z. Y. Lyu,  L. M. Lui, and G. P. T. Choi,
    "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
    Preprint, 2024.

Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi

===============================================================

Usage:
* `map = SDEM(v,f,population,dt,epsilon,max_iter)`
* `map = LSDEM(v,f,population,S,landmark,target,alpha,beta,gamma,dt,epsilon,max_iter)`

Input:
* `v`: nv x 3 vertex coordinates of a genus-0 triangle mesh
* `f`: nf x 3 triangulations of a genus-0 triangle mesh
* `population`: nf x 1 positive quantity
* `S`: nv x 3 vertex coordinates of the initial spherical conformal parameterization
* `dt`: step size
* `epsilon`: stopping parameter
* `max_iter`: maximum number of iterations
* `landmark`: k x 1 vertex indices of the landmarks
* `target`: k x 3 target positions of the landmarks on the unit sphere
* `alpha`: nonnegative weighting parameter for the density-equalizing term
* `beta`: nonnegative weighting parameter for the harmonic term
* `gamma`: nonnegative weighting parameter for the landmark mismatch term

Output:
* `map`: nv x 3 vertex coordinates of the spherical parameterization

