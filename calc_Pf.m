function res = calc_Pf(p,y,z,x_c);

% -- res = calc_Pf(p,y,z,x_c);
%
% The purpose of this function is to calculate the 
% fledging rate x clutch size given parameters, arrival
% time, prelaying period, and optimal hatching time. This
% function will accept either a single point or vector of
% points.
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
% res: Probability of fledging x clutch size
%

% Note that the clutch size parameter $a$ has been 
% included here
res = p.a*exp(-((y+z+p.z_n)-x_c).^2/(2*p.sigma^2));
