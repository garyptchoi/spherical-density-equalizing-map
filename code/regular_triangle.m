function regular_triangle = regular_triangle(f,v)
% Find the triangle with the most regular 1-ring neighborhood.
%
% If you use this code in your own work, please cite the following paper:
% [1] Z. Lyu, L. M. Lui, and G. P. T. Choi,
%     "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
%     Preprint, arXiv:2401.11795, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/spherical-density-equalizing-map


nv = length(v);
nf = length(f);

% face regularity
temp = v(reshape(f',1,length(f)*3),1:3);
e1 = sqrt(sum((temp(2:3:end,1:3) - temp(3:3:end,1:3))'.^2))';
e2 = sqrt(sum((temp(1:3:end,1:3) - temp(3:3:end,1:3))'.^2))';
e3 = sqrt(sum((temp(1:3:end,1:3) - temp(2:3:end,1:3))'.^2))';
R_f = abs(e1./(e1+e2+e3)-1/3)+...
    abs(e2./(e1+e2+e3)-1/3)+abs(e3./(e1+e2+e3)-1/3);

% create face-to-vertex matrix
row = [f(:,1); f(:,2); f(:,3)];
col = [1:length(f), 1:length(f), 1:length(f)]';
val = 1/3*ones(3*length(f),1);
H = sparse(row,col,val,nv,nf);

% vertex regularity
R_v = H*R_f;

% average vertex regularity for each face
R_average = 1/3*(R_v(f(:,1))+R_v(f(:,2))+R_v(f(:,3)));
[~,regular_triangle] = min(R_average);

end

