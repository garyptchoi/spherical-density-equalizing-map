% LSDEM: Computing landmark-aligned spherical density-equalizing maps of
% genus-0 closed surfaces using the proposed LSDEM method in [1].
% 
% Main program:
% map = LSDEM(v,f,population,S,landmark,target,alpha,beta,gamma,dt,epsilon,max_iter)
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed triangle mesh
% f: nf x 3 triangulations of a genus-0 closed triangle mesh
% population: nf x 1 positive quantity
% S: nv x 3 vertex coordinates of the initial spherical conformal parameterization 
%   (set S=[] if you want the algorithm to compute it automatically)
% landmark: k x 1 vertex indices of the landmarks
% target: k x 3 target positions of the landmarks on the unit sphere
% alpha: nonnegative weighting parameter for the density-equalizing term
% beta: nonnegative weighting parameter for the harmonic term
% gamma: nonnegative weighting parameter for the landmark mismatch term
% dt: step size (optional, default = 0.01)
% epsilon: stopping parameter (optional, default = 1e-3)
% max_iter: maximum number of iterations (optional, default = 200)
%
% Output:
% map: nv x 3 vertex coordinates of the landmark-aligned spherical density-equalizing map
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

load('sphere_with_landmark.mat'); 

% For an input spherical surface, we can directly run the LSDEM algorithm 
% with the initial spherical parameterization being itself:
S = v;
tic;
map = LSDEM(v,f,population,S,landmark,target,1,2,5,0.01,1e-3,200);
toc;

density = population./face_area(f,S);

plot_mesh_with_density(S,f,density);
hold on;
plot3(S(landmark,1),S(landmark,2),S(landmark,3),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10);
plot3(target(:,1),target(:,2),target(:,3),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10);
title('Input surface');

plot_mesh_with_density(map,f,density);
hold on;
plot3(map(landmark,1),map(landmark,2),map(landmark,3),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10);
plot3(target(:,1),target(:,2),target(:,3),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10);
title('LSDEM result');

%% Example 2: Mapping a general genus-0 closed surface

load('brain_with_landmark.mat'); 

% For a general input surface, we first compute an initial spherical
% conformal parameterization using the FLASH method in [Choi et al., SIAM J. Imaging Sci. 2015] 
% together with the Mobius area correction scheme in [Choi et al., SIAM J. Imaging Sci. 2020]
S1 = spherical_conformal_map(v,f);
S = mobius_area_correction_spherical(v,f,S1);

% Run the LSDEM algorithm
map = LSDEM(v,f,population,S,landmark,target,1,2,5,0.01,1e-3,200);

plot_mesh(v,f);
hold on;
plot3(v(landmark,1),v(landmark,2),v(landmark,3),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10);
title('Input surface');
view([90 25])

plot_mesh(S,f);
hold on;
plot3(S(landmark,1),S(landmark,2),S(landmark,3),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10);
plot3(target(:,1),target(:,2),target(:,3),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10);
title('Initial spherical conformal parameterization');
view([-5 5])

plot_mesh(map,f);
hold on;
plot3(map(landmark,1),map(landmark,2),map(landmark,3),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10);
plot3(target(:,1),target(:,2),target(:,3),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10);
title('LSDEM result');
view([-5 5])
