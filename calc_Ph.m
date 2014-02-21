function res = calc_Ph(p,y,z);

% -- res = calc_Ph(p,y,z);
%
% The purpose of this function is to calculate the 
% hatching success given parameters, arrival
% time, and prelaying period. This function will accept
% either a single point or vector of points.
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
%
% OUTPUTS
%
% res: Hatching success
%

u_q = p.u_q;
b_q = p.b_q;
Q_0 = p.Q_0;

intRyi = y-u_q + log(exp(-b_q*y+b_q*u_q)+1)/b_q;
C = Q_0./(exp(intRyi)*(1-Q_0));
intRt = (y+z)-u_q + log(exp(-b_q*(y+z)+b_q*u_q)+1)/b_q;

res = C.*exp(intRt)./(1+C.*exp(intRt));
