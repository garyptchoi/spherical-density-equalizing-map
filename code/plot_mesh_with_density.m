function plot_mesh_with_density(v,f,density)
% Plot a mesh with density.
% 
% Input: 
% v: nv x 3 vertex coordinates
% f: nf x 3 triangulations
% density: nf x 1 quantity defined on faces
% 
% If you use this code in your work, please cite the following paper:
% [1] Z. Lyu, L. M. Lui, and G. P. T. Choi,
%     "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
%     Preprint, arXiv:2401.11795, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/spherical-density-equalizing-map

figure; 
patch('Faces',f,'Vertices',v,'FaceColor','flat','FaceVertexCData',density);
colormap spring; colorbar
axis equal tight off
ax = gca; ax.Clipping = 'off';