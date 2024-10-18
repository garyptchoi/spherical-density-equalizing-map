function [v_temp,M,Minv] = rotate_sphere(f,v,north_f)
% Rotate the sphere such that the triangle north_f becomes the north pole.
% 
% If you use this code in your work, please cite the following paper:
% [1] Z. Lyu, L. M. Lui, and G. P. T. Choi,
%     "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
%     SIAM Journal on Imaging Sciences, 17(4), 2110-2141, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/spherical-density-equalizing-map


center_v = (v(f(north_f,1),:) + v(f(north_f,2),:) + v(f(north_f,3),:))/3;
center_v = center_v/sqrt(sum(center_v.^2));

sin_z = -center_v(2)/sqrt(sum(center_v(1)^2 + center_v(2)^2));
cos_z = center_v(1)/sqrt(sum(center_v(1)^2 + center_v(2)^2));
rot_z = [cos_z,-sin_z,0;sin_z,cos_z,0;0,0,1];
B = rot_z*[center_v(1),center_v(2),center_v(3)]';

b1 = B(1); b2 = B(2); b3 = B(3);

sin_y = -b1/sqrt(sum(b1^2 + b3^2));
cos_y = b3/sqrt(sum(b1^2 + b3^2));
rot_y = [cos_y,0,sin_y;0,1,0;-sin_y,0,cos_y];

M = rot_y*rot_z;

Minv = [cos_z,sin_z,0;-sin_z,cos_z,0;0,0,1]*[cos_y,0,-sin_y;0,1,0;sin_y,0,cos_y];

v_temp = (M*v')';