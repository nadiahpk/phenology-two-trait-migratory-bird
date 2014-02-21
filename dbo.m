function res = dbo(p,yz,x_c,YZ)

% -- res = dbo(p,yz,x_c,YZ)
% -- res = dbo(p,yz,x_c)
%
%
% The purpose of this function is to calculate the fitness
% gradient at the given parameters, arrival time, prelaying
% period, and optimal hatching time. It's used in solvedbo.m
% to find the evolutionarily singular strategy, and in
% numcheck_ess.m and numcheck_conv.m to estimate the
% ESS-stability and convergence stability numerically.
%
%
% INPUTS
%
% p: A dictionary of parameter values. Use the file params.m
% for an example of how to specify one of these.
% 
% [yz]: The variant strategy arrival time and prelaying
% period.
%
% [YZ]: The population strategy arrival time and prelaying
% period. If this is not specified, it is taken as equal to
% the variant strategy.
% 
% x_c: Optimal hatching time.
% 
%
% OUTPUTS
%
% res: [dy;dz] The gradient of the invasion fitness with
% respect to arrival time and prelaying period. In vector 
% form so it's easy to use with a solver.
%
%
% See also: solvedbo

% Pull out parameters, just makes code easier to read
K = p.K; b_e = p.b_e; a = p.a; sigma = p.sigma; Q_0 = p.Q_0;
b_q = p.b_q; u_q = p.u_q; L_m = p.L_m; b_s = p.b_s; M = p.M; z_f = p.z_f; z_n = p.z_n;

if nargin == 3
    % In a situation don't need to distinguish,
    % like when finding the singular strategy in solvedbo.m
    YZ = yz; 
end
% In other situations, like when numerically testing for ESS stability (see numcheck_ess.m) we do need to make a distinction between the variant and population traits.

y_i = yz(1); z_i = yz(2);
Y = YZ(1); Z = YZ(2);

% Evaluate steady state n
Ps = calc_Ps(p,Y,Z,x_c);
Pf = calc_Pf(p,Y,Z,x_c);
Ph = calc_Ph(p,Y,Z);
n = p.K*Ps*Pf*Ph/((1/p.M)-Ps);
p.n = n;

% Evaluate each component at the mutant trait
Ps = calc_Ps(p,y_i,z_i,x_c);
Pf = calc_Pf(p,y_i,z_i,x_c);
Ph = calc_Ph(p,y_i,z_i);
Pe = calc_Pe(p,Y,y_i);

% The expressions for the derivatives below were verified
% with Sage

dPsdy = -(b_s*e^(-(y_i + z_n + z_f + z_i)*b_s) - b_s*e^(-b_s*y_i))*L_m*e^((e^(-(y_i + z_n + z_f + z_i)*b_s) - e^(-b_s*y_i))*L_m/b_s)/b_s;
dPsdz = -L_m*e^(-(y_i + z_n + z_f + z_i)*b_s + (e^(-(y_i + z_n + z_f + z_i)*b_s) - e^(-b_s*y_i))*L_m/b_s);

dPfdy = (x_c - y_i - z_i - z_n)*a*e^(-1/2*(x_c - y_i - z_i - z_n)^2/sigma^2)/sigma^2;
dPfdz = dPfdy;

dPhdy = -(e^(-(y_i + z_i)*b_q + b_q*u_q)/(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1) - e^(b_q*u_q - b_q*y_i)/(e^(b_q*u_q - b_q*y_i) + 1))*Q_0*e^(z_i + log(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y_i) + 1)/b_q)/((Q_0 - 1)*(Q_0*e^(z_i + log(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y_i) + 1)/b_q)/(Q_0 - 1) - 1)) + (e^(-(y_i + z_i)*b_q + b_q*u_q)/(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1) - e^(b_q*u_q - b_q*y_i)/(e^(b_q*u_q - b_q*y_i) + 1))*Q_0^2*e^(2*z_i + 2*log(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1)/b_q - 2*log(e^(b_q*u_q - b_q*y_i) + 1)/b_q)/((Q_0 - 1)^2*(Q_0*e^(z_i + log(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y_i) + 1)/b_q)/(Q_0 - 1) - 1)^2);

dPhdz = -(e^(-(y_i + z_i)*b_q + b_q*u_q)/(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1) - 1)*Q_0*e^(z_i + log(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y_i) + 1)/b_q)/((Q_0 - 1)*(Q_0*e^(z_i + log(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y_i) + 1)/b_q)/(Q_0 - 1) - 1)) + (e^(-(y_i + z_i)*b_q + b_q*u_q)/(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1) - 1)*Q_0^2*e^(2*z_i + 2*log(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1)/b_q - 2*log(e^(b_q*u_q - b_q*y_i) + 1)/b_q)/((Q_0 - 1)^2*(Q_0*e^(z_i + log(e^(-(y_i + z_i)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y_i) + 1)/b_q)/(Q_0 - 1) - 1)^2);

% NOTE Density dependence here
if Pe == 1
    dPedy = 0;
else
    t_h = Y + (1/b_e)*log((K/n)/(1-(K/n)));
    dPedy = -b_e*exp(b_e*t_h+b_e*y_i)/((exp(b_e*t_h)+exp(b_e*y_i))^2);
end

% Put the components together to find the fitness gradient
dy = M*(dPsdy + dPedy*Ps*Ph*Pf + Pe*dPsdy*Ph*Pf + Pe*Ps*dPhdy*Pf + Pe*Ps*Ph*dPfdy);
dz = M*(dPsdz + Pe*( dPhdz*Pf*Ps + Ph*dPfdz*Ps + Ph*Pf*dPsdz));

res = [dy; dz];
