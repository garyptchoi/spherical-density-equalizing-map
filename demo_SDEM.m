% SDEM: Computing spherical density-equalizing maps of genus-0 closed 
% surfaces using the proposed SDEM method in [1].
%
% Main program:
% map = SDEM(v,f,population,S,dt,epsilon,max_iter)
% 
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed triangle mesh
% f: nf x 3 triangulations of a genus-0 closed triangle mesh
% population: nf x 1 positive quantity
% S: nv x 3 vertex coordinates of the initial spherical conformal parameterization 
%   (set S=[] if you want the algorithm to compute it automatically)
% dt: step size (optional, default = 0.1)
% epsilon: stopping parameter (optional, default = 1e-3)
% max_iter: maximum number of iterations (optional, default = 200)
%
% Output:
% map: nv x 3 vertex coordinates of the spherical density-equalizing map
% 
% If you use this code in your own work, please cite the following paper:
% [1] Z. Lyu, L. M. Lui, and G. P. T. Choi,
%     "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
%     Preprint, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi

addpath('code');
addpath('data');

%% Example 1: Mapping a spherical surface

load('sphere.mat')

% For an input spherical surface, we can directly run the SDEM algorithm 
% with the initial spherical parameterization being itself:
S = v;
map = SDEM(v,f,population,S,0.1,1e-3,200);

density = population./face_area(f,S);

plot_mesh_with_density(v,f,density);
title('Input surface');

plot_mesh_with_density(map,f,density);
title('SDEM result');

%% Example 2: Mapping a general genus-0 closed surface

load('david1.mat'); % population: enlarging the nose region
% load('david2.mat'); % population: shrinking the nose region

% For a general input surface, we first compute an initial spherical
% conformal parameterization using the FLASH method in [Choi et al., SIAM J. Imaging Sci. 2015] 
% together with the Mobius area correction scheme in [Choi et al., SIAM J. Imaging Sci. 2020]
S1 = spherical_conformal_map(v,f);
S = mobius_area_correction_spherical(v,f,S1);

% We can then run the SDEM algorithm with the initial spherical conformal parameterization
map = SDEM(v,f,population,S,0.1,1e-3,200);

density = population./face_area(f,S);

plot_mesh_with_density(v,f,density);
title('Input surface');
view([-120 -5])

plot_mesh_with_density(S,f,density);
title('Initial spherical conformal parameterization');
view([-150 90])

plot_mesh_with_density(map,f,density);
title('SDEM result');
view([-150 90])

%% Example 3: Area-preserving parameterization of genus-0 closed surfaces

load('twistedball.mat');
% load('maxplanck.mat');

% First compute an initial spherical conformal parameterization
S1 = spherical_conformal_map(v,f);
S = mobius_area_correction_spherical(v,f,S1);

% Set the population as the face area for achieving area-preserving map
population = face_area(f,v);

% Run the SDEM algorithm
map = SDEM(v,f,population,S,0.1,1e-3,200);

plot_mesh(v,f);
title('Input surface');
view([-180 -80])

plot_mesh(S,f);
title('Initial spherical conformal parameterization');
view([90 0])

plot_mesh(map,f);
title('SDEM result');
view([90 0])


