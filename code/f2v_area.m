function M = f2v_area(v,f)
% Face to vertex interpolation with area weighting.
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi and C. H. Rycroft, 
%     "Density-equalizing maps for simply connected open surfaces."
%     SIAM Journal on Imaging Sciences, 11(2), pp. 1134-1178, 2018.
% 
% Copyright (c) 2017-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

nv = length(v);
nf = length(f);

if size(v,2) == 2
    v = [v,zeros(nv,1)];
end

% find area
area = face_area(f,v);

% create matrix
row = [f(:,3); f(:,1); f(:,2)];
col = [1:length(f), 1:length(f), 1:length(f)]';
val = [area; area; area];
M = sparse(row,col,val,nv,nf);
 
% normalize
vertex_area_sum = sum(M,2);
[Mrow,Mcol,Mval] = find(M);
M = sparse(Mrow,Mcol,Mval./vertex_area_sum(Mrow),nv,nf);