function map = LSDEM(v,f,population,S,landmark,target,alpha,beta,gamma,dt,epsilon,max_iter)
% Computing landmark-aligned spherical density-equalizing maps of genus-0 
% closed surfaces using the proposed LSDEM method in [1].
% 
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed triangle mesh
% f: nf x 3 triangulations of a genus-0 closed triangle mesh
% population: nf x 1 positive quantity
% S: nv x 3 vertex coordinates of the initial spherical conformal parameterization 
%   (set S=[] if you want the algorithm to compute it automatically)
% landmark: k x 1 vertex indices of the landmarks
% target: k x 3 target positions of the landmarks on the unit sphere
% alpha: nonnegative weighting parameter for the density-equalizing term
% beta: nonnegative weighting parameter for the harmonic term
% gamma: nonnegative weighting parameter for the landmark mistmatch term
% dt: step size (optional, default = 0.01)
% epsilon: stopping parameter (optional, default = 1e-3)
% max_iter: maximum number of iterations (optional, default = 200)
%
% Output:
% map: nv x 3 vertex coordinates of the landmark-aligned spherical density-equalizing map
%
% If you use this code in your own work, please cite the following paper:
% [1] Z. Lyu, L. M. Lui, and G. P. T. Choi,
%     "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
%     Preprint, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi


if isempty(S)
    % Compute an initial spherical conformal parameterization 
    
    % Apply the FLASH method in [Choi et al., SIAM J. Imaging Sci. 2015] 
    S1 = spherical_conformal_map(v,f);
    
    % Reduce the area distortion of the parameterization using the Mobius 
    % area correction scheme in [Choi et al., SIAM J. Imaging Sci. 2020]
    S = mobius_area_correction_spherical(v,f,S1);
end

if nargin < 10
    dt = 0.01;
end

if nargin < 11
    epsilon = 1e-3;
end

if nargin < 12
    max_iter = 200;
end

if alpha < 0 || beta < 0 || gamma < 0
    error('The weighting parameters should be nonnegative.');
end

if max(abs(sqrt(sum(target.^2,2))-1)) > 1e-4
    error('The target positions should lie on the unit sphere.');
end

% normalize the input spherical parameterization
r = [S(:,1)./sqrt(sum(S.^2,2)), ...
     S(:,2)./sqrt(sum(S.^2,2)), ...
     S(:,3)./sqrt(sum(S.^2,2))];

bigtri = regular_triangle(f,r); 

% compute density
population = population/sum(population); 
rho_f = population./face_area(f,r);
rho_v = f2v_area(r,f)*rho_f;

step = 0;
f_diff = Inf;
disp('Step     ||f_n-f_{n-1}||');

r_old = r;
while f_diff >= epsilon && step < max_iter
    
    % optimal rotation
    r = optimal_rotation(r,landmark,target);

    % update rho 
    L = laplace_beltrami(r,f);
    A = lumped_mass_matrix(r,f);
    rho_v_temp = (A+dt*L)\(A*rho_v);
    
    % update density gradient
    grad_rho_temp_f = compute_gradient_3D(r,f,rho_v_temp);
    grad_rho_temp_v = f2v_area(r,f)*grad_rho_temp_f;
    
    % update displacement using gradient descent
    dr = -[grad_rho_temp_v(:,1)./rho_v_temp,...
           grad_rho_temp_v(:,2)./rho_v_temp,...
           grad_rho_temp_v(:,3)./rho_v_temp];
    dE1 = alpha*(dr-[sum(dr.*r,2),sum(dr.*r,2),sum(dr.*r,2)].*r); % density term
    Lr = -A\(L*r);
    dE2 = -beta*(Lr-[sum(Lr.*r,2),sum(Lr.*r,2),sum(Lr.*r,2)].*r); % harmonic term
    dE3 = zeros(size(r));
    dE3(landmark,:) = -gamma*(r(landmark,:)-target); % landmark mismatch term
    dE = dE1 + dE2 + dE3;
    
    % overlap correction
    r = update_and_correct_overlap(f,S,r,bigtri,dE,dt);
    
    step = step + 1;
    f_diff = norm(r-r_old,Inf);
    disp([num2str(step), '        ',num2str(f_diff)]);
    
    % re-coupling scheme
    rho_f = population./face_area(f,r);
    rho_v = f2v_area(r,f)*rho_f;
    
    r_old = r;
end

map = r;

end