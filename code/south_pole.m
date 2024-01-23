function south_f = south_pole(f,v,bigtri)
% Find the south pole with respect to a given north pole triangle.
% 
% If you use this code in your work, please cite the following paper:
% [1] Z. Lyu, L. M. Lui, and G. P. T. Choi,
%     "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
%     Preprint, arXiv:2401.11795, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/spherical-density-equalizing-map

% get face centers
f_center = (v(f(:,1),:)+v(f(:,2),:)+v(f(:,3),:))/3;

% project onto the sphere
radius = sqrt(sum(f_center.^2,2));
f_center = f_center./[radius,radius,radius];

% find the most distant one
[~,south_f] = max(sum((f_center-f_center(bigtri)).^2,2));





