function res = calc_Pe(p,Y,y)

% -- res = calc_Pe(p,Y,y)
%
%
% The purpose of this function is to calculate the 
% territory acquisition probability given parameters, arrival
% time in population, arrival time of variant. This
% function will accept either a single point or vector of
% points.
%
%
% INPUTS
%
% p: A dictionary of parameter values. Use the file params.m
% for an example of how to specify one of these.
% 
% Y: Arrival time of population ("resident strategy"),
% either single value or vector.
% 
% Y: Arrival time of variant ("mutant strategy"),
% either single value or vector.
% 
%
% OUTPUTS
%
% res: Territory acquisition probability of variant strategy
%

K = p.K;
b_e = p.b_e;
n = p.n;
if n < K
    res = 1; 
else
    if nargin < 3;
        % No variant, just want resident strategy value
        res = K/n; % It's equiv to below when Y=y but shorter
    else
        % Want variant in population
        t_h = Y + log(-K/((K/n-1)*n))/b_e;
        res = e.^(b_e*t_h)./(e.^(b_e*t_h)+e.^(b_e*y));
        % FYI, this is also equivalent to
        % res = 1/(1+exp(b_e*(yd-t_h)));
        % So Pe(y'-y,n(y,z)).
    end
end
