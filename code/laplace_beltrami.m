function L = laplace_beltrami(v,f)
% Compute the cotangent Laplacian.
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi and C. H. Rycroft, 
%     "Density-equalizing maps for simply connected open surfaces."
%     SIAM Journal on Imaging Sciences, 11(2), pp. 1134-1178, 2018.
% 
% Copyright (c) 2017-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

nv = length(v);

f1 = f(:,1); 
f2 = f(:,2); 
f3 = f(:,3);

% edge length
l = [sqrt(sum((v(f2,:) - v(f3,:)).^2,2)),...
    sqrt(sum((v(f3,:) - v(f1,:)).^2,2)),...
    sqrt(sum((v(f1,:) - v(f2,:)).^2,2))];
l1 = l(:,1); 
l2 = l(:,2); 
l3 = l(:,3);

% Heron's formula
s = (l1 + l2 + l3)*0.5;
area = sqrt( s.*(s-l1).*(s-l2).*(s-l3));
 
% cotangent weight
cot12 = (l1.^2 + l2.^2 - l3.^2)./area/4;
cot23 = (l2.^2 + l3.^2 - l1.^2)./area/4; 
cot31 = (l1.^2 + l3.^2 - l2.^2)./area/4; 

% construct matrix
II = [f1; f2; f2; f3; f3; f1; f1; f2; f3];
JJ = [f2; f1; f3; f2; f1; f3; f1; f2; f3];
V = [-cot12; -cot12; -cot23; -cot23; -cot31; -cot31; ...
    cot12+cot31; cot12+cot23; cot31+cot23]/2;
L = sparse(II,JJ,V,nv,nv);

end