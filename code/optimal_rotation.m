function v = optimal_rotation(v,landmark,target)
% Optimal rotation to reduce landmark mismatch.
% 
% If you use this code in your work, please cite the following paper:
% [1] Z. Lyu, L. M. Lui, and G. P. T. Choi,
%     "Spherical Density-Equalizing Map for Genus-0 Closed Surfaces."
%     Preprint, arXiv:2401.11795, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/spherical-density-equalizing-map

S = v;
S_landmark = S(landmark,:);

% rotation matrices
R_x = @(t) [1,0,0;0,cos(t),-sin(t);0,sin(t),cos(t)];
R_y = @(g) [cos(g),0,sin(g);0,1,0;-sin(g),0,cos(g)];
R_z = @(h) [cos(h),-sin(h),0;sin(h),cos(h),0;0,0,1];

% derivatives of rotation matrices
dR_x = @(t) [0,0,0;0,-sin(t),-cos(t);0,cos(t),-sin(t)];
dR_y = @(g) [-sin(g),0,cos(g);0,0,0;-cos(g),0,-sin(g)];
dR_z = @(h) [-sin(h),-cos(h),0;cos(h),-sin(h),0;0,0,0];

% Landmark mismatch error
L = @(w) sum(sum((R_x(w(1))*R_y(w(2))*R_z(w(3))*S_landmark')' - target).^2,2);

% initialization
E = 1;
para_t = 0;
para_g = 0;
para_h = 0;
dt = 0.001;
step = 0;
L_old = L([para_t,para_g,para_h]);

while E > 1e-6 && step < 1000
    % update parameters
    para_t = para_t - dt*...
        sum(sum(2*((R_x(para_t)*R_y(para_g)*R_z(para_h)*S_landmark')' - target).*...
        (dR_x(para_t)*R_y(para_g)*R_z(para_h)*S_landmark')'));
    
    para_g = para_g - dt*...
        sum(sum(2*((R_x(para_t)*R_y(para_g)*R_z(para_h)*S_landmark')' - target).*...
        (R_x(para_t)*dR_y(para_g)*R_z(para_h)*S_landmark')'));
    
    para_h = para_h - dt*...
        sum(sum(2*((R_x(para_t)*R_y(para_g)*R_z(para_h)*S_landmark')' - target).*...
        (R_x(para_t)*R_y(para_g)*dR_z(para_h)*S_landmark')'));
    
    % update landmark mismatch error
    L_temp = L([para_t,para_g,para_h]);
    
    E = abs(L_temp - L_old);
    
    L_old = L_temp;
    
    step = step + 1;
end

% final rotation result
v =  (R_x(para_t)*R_y(para_g)*R_z(para_h)*S')';

end


