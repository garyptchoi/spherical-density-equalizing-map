function grad = compute_gradient_3D(v,f,g)
% compute the 3D gradient.
%
% If you use this code in your own work, please cite the following paper:
% [1] Z. Lyu, L. M. Lui, and G. P. T. Choi,
%     "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
%     Preprint, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi

if size(v,2) ~= 3
    v = [v, v(:,1).*0];
end

% compute edges
e1 = v(f(:,3),:) - v(f(:,2),:);
e2 = v(f(:,1),:) - v(f(:,3),:);
e3 = v(f(:,2),:) - v(f(:,1),:);

% compute area
cross12 = cross(e1,e2);
area = abs(1/2*(cross12(:,1).^2+cross12(:,2).^2+cross12(:,3).^2).^(1/2));
N = [(1./(2*area)).*cross12(:,1), (1./(2*area)).*cross12(:,2), (1./(2*area)).*cross12(:,3)];

% compute gradient
temp = [g(f(:,1)).*e1(:,1),g(f(:,1)).*e1(:,2),g(f(:,1)).*e1(:,3)] + ...
    [g(f(:,2)).*e2(:,1),g(f(:,2)).*e2(:,2),g(f(:,2)).*e2(:,3)] + ...
    [g(f(:,3)).*e3(:,1),g(f(:,3)).*e3(:,2),g(f(:,3)).*e3(:,3)];
grad = cross(N,temp);
grad = [1./(2*area).*grad(:,1),1./(2*area).*grad(:,2),1./(2*area).*grad(:,3)];

