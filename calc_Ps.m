function res = calc_Ps(p,y,z);

% -- res = calc_Ps(p,y,z);
%
% The purpose of this function is to calculate the 
% early-season adult survival probability given 
% parameters, arrival time, prelaying period, and optimal
% hatching time. This function will accept either a single
% point or vector of points.
%
%
% INPUTS
%
% p: A dictionary of parameter values. Use the file params.m
% for an example of how to specify one of these.
% 
% y: Arrival time, either single value or vector.
% 
% z: Prelaying period, either single value or vector.
% 
% x_c: Optimal hatching time, either single value or vector.
%
%
% OUTPUTS
%
% res: Early-season adult survival probability
%

L_m = p.L_m;
b_s = p.b_s;
z_n = p.z_n;
z_f = p.z_f;

res = exp( (L_m/b_s)*( exp(-b_s*(y+z+z_n+z_f))-exp(-b_s*y) ));
