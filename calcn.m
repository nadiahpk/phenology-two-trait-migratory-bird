function nV = calcn(p,yzV,x_cV);

% -- nV = calcn(p,yzV,x_cV);
%
%
% The purpose of this function is to calculate the steady
% state population value for given parameters, arrival time,
% prelaying period, and optimal hatching time. This function
% will accept either a single point or a vector of points.
%
% 
% INPUTS
%
% p: A dictionary of parameter values. Use the file params.m
% for an example of how to specify one of these.
% 
% yzV: [Arrival time,Prelaying period], either single row or
% matrix.
% 
% x_cV: Optimal hatching time, either single value or vector.
%
%
% OUTPUTS
%
% nV: Population steady state, either single value or
% vector.
%


nV = [];

for ind = 1:length(x_cV)
    y_i = yzV(ind,1);
    z_i = yzV(ind,2);
    x_c = x_cV(ind);

    % Calculate steady-state values 
    Pf = calc_Pf(p,y_i,z_i,x_c);
    Ps = calc_Ps(p,y_i,z_i);
    Ph = calc_Ph(p,y_i,z_i);
    n = p.K*Ps*Pf*Ph/((1/p.M)-Ps);

    nV = [nV;n];
end
