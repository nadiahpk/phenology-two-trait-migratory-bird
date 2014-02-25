function [x_cV,yzV,nV] = carlo(u_q);

% -- [x_cV,yzV,nV] = carlo(u_q);
%
%
% The purpose of this function is to create vectors ready
% for a python slider that Carlo wants to create.
%
% INPUTS
%
% u_q: A parameter for the slider.
%
%
% OUTPUTS
%
% x_cV: A vector of x_c values in the range specified. This
% is the x-axis of Figure 1 a, b, and c.
%
% yzV: A matrix of [y*,z*] values in the range specified,
% corresponding to x_cV for easy plotting. This is used to
% create the various curves and filled areas in the bottom
% pane of Figure 1a,b, and c, the y-axis marked "bird 
% phenology".
%
% nV: A vector of population sizes corresponding to x_cV. It
% is the y-axis of the top pane of Figure 1a, b, and c.
%

params; % Get dictionary of parameter values p
change_p = 'u_q';
change_p_val = u_q;
x_c_hi = 175; x_c_lo = 135;
y0 = 130; 
z0 = x_c_hi-y0-p.z_n; 
print_flag = 0;
plot_flag = 0;

[x_cV,yzV,nV] = solvedbo(p,change_p,change_p_val,x_c_hi,x_c_lo,y0,z0,print_flag,plot_flag);
