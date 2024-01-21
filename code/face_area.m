function fa = face_area(f,v)
% Compute the area of all faces.
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi and C. H. Rycroft, 
%     "Density-equalizing maps for simply connected open surfaces."
%     SIAM Journal on Imaging Sciences, 11(2), pp. 1134-1178, 2018.
% 
% Copyright (c) 2017-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

v12 = v(f(:,2),:) - v(f(:,1),:);
v23 = v(f(:,3),:) - v(f(:,2),:);
v31 = v(f(:,1),:) - v(f(:,3),:);

a = sqrt(dot(v12,v12,2));
b = sqrt(dot(v23,v23,2));
c = sqrt(dot(v31,v31,2));

s = (a+b+c)/2.0;
fa = sqrt(s.*(s-a).*(s-b).*(s-c)); 
