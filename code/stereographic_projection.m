function p = stereographic_projection(v)
% Stereographic projection and inverse stereographic projection.
% 
% If you use this code in your own work, please cite the following paper:
% [1] Z. Lyu, L. M. Lui, and G. P. T. Choi,
%     "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
%     Preprint, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi

    if size(v,2) == 3
        % stereographic projection
        p = [v(:,1)./(1-v(:,3)), v(:,2)./(1-v(:,3))];
        p(isnan(p)) = Inf;
    else
        if size(v,2) == 1
        % turn the complex vector into xy coordinates
        v = [real(v), imag(v)];
        end
        z = 1+v(:,1).^2+v(:,2).^2;
        % inverse stereographic projection
        p = [2*v(:,1)./z, 2*v(:,2)./z, (-1+v(:,1).^2+v(:,2).^2)./z];
        p(isnan(z)|(~isfinite(z)),1) = 0;
        p(isnan(z)|(~isfinite(z)),2) = 0;
        p(isnan(z)|(~isfinite(z)),3) = 1;
    end
end