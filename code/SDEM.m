function map = SDEM(v,f,population,S,dt,epsilon,max_iter)
% Computing spherical density-equalizing maps of genus-0 closed surfaces 
% using the proposed SDEM method in [1].
%
% Input:
% v: nv x 3 vertex coordinates of a spherical surface mesh
% f: nf x 3 triangulations of a spherical surface mesh
% population: nf x 1 positive quantity
% S: nv x 3 vertex coordinates of the initial spherical conformal parameterization 
%   (set S=[] if you want the algorithm to compute it automatically)
% dt: step size (optional, default = 0.1)
% epsilon: stopping parameter (optional, default = 1e-3)
% max_iter: maximum number of iterations (optional, default = 200)
%
% Output:
% map: nv x 3 vertex coordinates of the spherical density-equalizing map
% 
% If you use this code in your work, please cite the following paper:
% [1] Z. Lyu, L. M. Lui, and G. P. T. Choi,
%     "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
%     SIAM Journal on Imaging Sciences, 17(4), 2110-2141, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/spherical-density-equalizing-map


if isempty(S)
    % Compute an initial spherical conformal parameterization 
    
    % Apply the FLASH method in [Choi et al., SIAM J. Imaging Sci. 2015] 
    S1 = spherical_conformal_map(v,f);
    
    % Reduce the area distortion of the parameterization using the Mobius 
    % area correction scheme in [Choi et al., SIAM J. Imaging Sci. 2020]
    S = mobius_area_correction_spherical(v,f,S1);
end

if nargin < 5
    dt = 0.1;
end

if nargin < 6
    epsilon = 1e-3;
end

if nargin < 7
    max_iter = 200;
end

% normalize the input spherical parameterization
r = [S(:,1)./sqrt(sum(S.^2,2)),...
     S(:,2)./sqrt(sum(S.^2,2)),...
     S(:,3)./sqrt(sum(S.^2,2))];
 
bigtri = regular_triangle(f,r);

% compute density
rho_f = population./face_area(f,r);
rho_v = f2v_area(r,f)*rho_f;

step = 0;
rho_v_error = std(rho_v)/mean(rho_v);
disp('Step     std(rho)/mean(rho)');
disp([num2str(step), '        ',num2str(rho_v_error)]);
    
while rho_v_error >= epsilon && step < max_iter
    
    % update rho 
    L = laplace_beltrami(r,f);
    A = lumped_mass_matrix(r,f);
    rho_v_temp = (A+dt*L)\(A*rho_v);

    % update density gradient
    grad_rho_temp_f = compute_gradient_3D(r,f,rho_v_temp);
    grad_rho_temp_v = f2v_area(r,f)*grad_rho_temp_f;
    
    % update displacement
    dr = -[grad_rho_temp_v(:,1)./rho_v_temp,...
           grad_rho_temp_v(:,2)./rho_v_temp,...
           grad_rho_temp_v(:,3)./rho_v_temp];
       
    dr_proj = dr-[sum(dr.*r,2),sum(dr.*r,2),sum(dr.*r,2)].*r;
    
    r = update_and_correct_overlap(f,S,r,bigtri,dr_proj,dt);
    
    step = step + 1;
    rho_v_error = std(rho_v_temp)/mean(rho_v_temp);
    disp([num2str(step), '        ',num2str(rho_v_error)]);

    % re-coupling scheme
    rho_f = population./face_area(f,r);
    rho_v = f2v_area(r,f)*rho_f;
    
end

map = r;

end

