function [resV,Hess] = check_ess(p,yzV,x_cV)

% -- [resV,Hess] = check_ess(p,yzV,x_cV)
%
% The purpose of this function is to check the ESS
% stability of each evolutionarily singular strategy defined
% by the inputs. The check is done by finding the largest
% eigenvalue of the Hessian matrix and returning it to the
% user, who can verify that they are zero.
% 
% This code uses a semi-analytic approach in
% that the expressions for each of the derivatives and
% second derivatives were found using Sage. This code also
% has a counterpart numcheck_ess.m which finds the
% eigenvalues numerically. Agreement between the two methods
% was used a check of the correctness of the solutions.
%
%
% INPUTS
%
% p: A dictionary of parameter values. Use the file params.m
% for an example of how to specify one of these.
%
% yzV: A matrix of [y*,z*] values corresponding to x_cV for 
% easy plotting.
%
% x_cV: A vector of x_c values at which [y*,z*] were
% calculated.
%
%
% OUTPUTS
%
% resV: A vector of largest eigenvalues of the Hessian
% matrix corresponding to the evolutionarily singular
% strategies specified by the inputs.
%
% Hess: The Jacobian matrix of the last evolutionarily
% singular strategy specified. It is useful if you want to
% take a closer look at just one singular strategy.
%
%
% see also: numcheck_ess, check_conv
%

% Pull out parameters, just makes code easier to read
K = p.K; b_e = p.b_e; a = p.a; sigma = p.sigma; Q_0 = p.Q_0; b_q = p.b_q; u_q = p.u_q; L_m = p.L_m; b_s = p.b_s; M = p.M; z_f = p.z_f; z_n = p.z_n;

resV = []; % Storage vector for output

% This loop steps through each of our x_c values
for ind = 1:length(x_cV);

    yz = yzV(ind,:);
    x_c = x_cV(ind);
    y = yz(1); z = yz(2);

    % Calculate each component's value
    Ps = calc_Ps(p,y,z);
    Pf = calc_Pf(p,y,z,x_c);
    Ph = calc_Ph(p,y,z);
    n = calcn(p,yz,x_c);
    Pe = K/n;

    % Calculate first derivatives 
    % These expressions were found using Sage. 

    dPsdy = -(b_s*e^(-(y + z_n + z_f + z)*b_s) - b_s*e^(-b_s*y))*L_m*e^((e^(-(y + z_n + z_f + z)*b_s) - e^(-b_s*y))*L_m/b_s)/b_s;
    dPsdz = -L_m*e^(-(y + z_n + z_f + z)*b_s + (e^(-(y + z_n + z_f + z)*b_s) - e^(-b_s*y))*L_m/b_s);

    dPfdy = (x_c - y - z - z_n)*a*e^(-1/2*(x_c - y - z - z_n)^2/sigma^2)/sigma^2;
    dPfdz = dPfdy;

    dPhdy = -(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - e^(b_q*u_q - b_q*y)/(e^(b_q*u_q - b_q*y) + 1))*Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)) + (e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - e^(b_q*u_q - b_q*y)/(e^(b_q*u_q - b_q*y) + 1))*Q_0^2*e^(2*z + 2*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 2*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^2*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^2);

    dPhdz = -(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - 1)*Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)) + (e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - 1)*Q_0^2*e^(2*z + 2*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 2*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^2*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^2);

    t_h = y + log(-K/((K/n-1)*n))/b_e;
    dPedy = -b_e*e.^((y + log(-Pe/(Pe-1))/b_e)*b_e + b_e*y)./(e^((y + log(-Pe/(Pe-1))/b_e)*b_e) + e.^(b_e*y))^2;
    dPedz = 0;

    % Calculate second derivatives at yz steady states


    dPsdydy = (b_s*e^(-(y + z + z_f + z_n)*b_s) - b_s*e^(-b_s*y))^2*L_m^2*e^((e^(-(y + z + z_f + z_n)*b_s) - e^(-b_s*y))*L_m/b_s)/b_s^2 + (b_s^2*e^(-(y + z + z_f + z_n)*b_s) - b_s^2*e^(-b_s*y))*L_m*e^((e^(-(y + z + z_f + z_n)*b_s) - e^(-b_s*y))*L_m/b_s)/b_s;

    dPsdydz = (b_s*e^(-(y + z + z_f + z_n)*b_s) - b_s*e^(-b_s*y))*L_m^2*e^(-(y + z + z_f + z_n)*b_s + (e^(-(y + z + z_f + z_n)*b_s) - e^(-b_s*y))*L_m/b_s)/b_s + L_m*b_s*e^(-(y + z + z_f + z_n)*b_s + (e^(-(y + z + z_f + z_n)*b_s) - e^(-b_s*y))*L_m/b_s);

    dPsdzdz = (L_m*e^(-(y + z + z_f + z_n)*b_s) + b_s)*L_m*e^(-(y + z + z_f + z_n)*b_s + (e^(-(y + z + z_f + z_n)*b_s) - e^(-b_s*y))*L_m/b_s);


    dPfdydy = (x_c - y - z - z_n)^2*a*e^(-1/2*(x_c - y - z - z_n)^2/sigma^2)/sigma^4 - a*e^(-1/2*(x_c - y - z - z_n)^2/sigma^2)/sigma^2;
    dPfdydz = dPfdydy;
    dPfdzdz = dPfdydy;


    dPhdydy = (e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - e^(b_q*u_q - b_q*y)/(e^(b_q*u_q - b_q*y) + 1))^2*Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)) + (b_q*e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - b_q*e^(b_q*u_q - b_q*y)/(e^(b_q*u_q - b_q*y) + 1) - b_q*e^(-2*(y + z)*b_q + 2*b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1)^2 + b_q*e^(2*b_q*u_q - 2*b_q*y)/(e^(b_q*u_q - b_q*y) + 1)^2)*Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)) - 3*(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - e^(b_q*u_q - b_q*y)/(e^(b_q*u_q - b_q*y) + 1))^2*Q_0^2*e^(2*z + 2*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 2*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^2*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^2) - (b_q*e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - b_q*e^(b_q*u_q - b_q*y)/(e^(b_q*u_q - b_q*y) + 1) - b_q*e^(-2*(y + z)*b_q + 2*b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1)^2 + b_q*e^(2*b_q*u_q - 2*b_q*y)/(e^(b_q*u_q - b_q*y) + 1)^2)*Q_0^2*e^(2*z + 2*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 2*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^2*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^2) + 2*(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - e^(b_q*u_q - b_q*y)/(e^(b_q*u_q - b_q*y) + 1))^2*Q_0^3*e^(3*z + 3*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 3*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^3*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^3);

    dPhdydz = (e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - 1)*(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - e^(b_q*u_q - b_q*y)/(e^(b_q*u_q - b_q*y) + 1))*Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)) - 3*(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - 1)*(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - e^(b_q*u_q - b_q*y)/(e^(b_q*u_q - b_q*y) + 1))*Q_0^2*e^(2*z + 2*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 2*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^2*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^2) + (b_q*e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - b_q*e^(-2*(y + z)*b_q + 2*b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1)^2)*Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)) + 2*(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - 1)*(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - e^(b_q*u_q - b_q*y)/(e^(b_q*u_q - b_q*y) + 1))*Q_0^3*e^(3*z + 3*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 3*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^3*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^3) - (b_q*e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - b_q*e^(-2*(y + z)*b_q + 2*b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1)^2)*Q_0^2*e^(2*z + 2*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 2*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^2*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^2);

    dPhdzdz = (e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - 1)^2*Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)) - 3*(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - 1)^2*Q_0^2*e^(2*z + 2*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 2*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^2*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^2) + (b_q*e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - b_q*e^(-2*(y + z)*b_q + 2*b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1)^2)*Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)) + 2*(e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - 1)^2*Q_0^3*e^(3*z + 3*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 3*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^3*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^3) - (b_q*e^(-(y + z)*b_q + b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1) - b_q*e^(-2*(y + z)*b_q + 2*b_q*u_q)/(e^(-(y + z)*b_q + b_q*u_q) + 1)^2)*Q_0^2*e^(2*z + 2*log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - 2*log(e^(b_q*u_q - b_q*y) + 1)/b_q)/((Q_0 - 1)^2*(Q_0*e^(z + log(e^(-(y + z)*b_q + b_q*u_q) + 1)/b_q - log(e^(b_q*u_q - b_q*y) + 1)/b_q)/(Q_0 - 1) - 1)^2);

    dPedydy = -b_e^2*e^(b_e*t_h + b_e*y)/(e^(b_e*t_h) + e^(b_e*y))^2 + 2*b_e^2*e^(b_e*t_h + 2*b_e*y)/(e^(b_e*t_h) + e^(b_e*y))^3;
    dPedydz = 0;
    dPedzdz = 0;

    % Put Hessian Matrix together 
    Hess = zeros(2,2); % order: [dydy dydz; dydz dzdz]

    Hess(1,1) = ...
        dPsdydy*Pe*Ph*Pf + Ps*dPedydy*Ph*Pf + Ps*Pe*dPhdydy*Pf + Ps*Pe*Ph*dPfdydy + ...
        2*dPsdy*dPedy*Ph*Pf + 2*dPsdy*Pe*dPhdy*Pf + 2*dPsdy*Pe*Ph*dPfdy + ...
                             2*Ps*dPedy*dPhdy*Pf + 2*Ps*dPedy*Ph*dPfdy + ...
                                                  2*Ps*Pe*dPhdy*dPfdy + ...
        dPsdydy;

    Hess(1,2) = ...
        dPsdydz*Pe*Ph*Pf + Ps*dPedydz*Ph*Pf + Ps*Pe*dPhdydz*Pf + Ps*Pe*Ph*dPfdydz + ...
        dPsdy*dPedz*Ph*Pf + dPsdz*dPedy*Ph*Pf + dPsdy*Pe*dPhdz*Pf + dPsdz*Pe*dPhdy*Pf + dPsdy*Pe*Ph*dPfdz + dPsdz*Pe*Ph*dPfdy + Ps*dPedy*dPhdz*Pf + Ps*dPedz*dPhdy*Pf + Ps*dPedy*Ph*dPfdz + Ps*dPedz*Ph*dPfdy + Ps*Pe*dPhdy*dPfdz + Ps*Pe*dPhdz*dPfdy + ...
        dPsdydz;

    Hess(2,1) = Hess(1,2);

    Hess(2,2) = ...
        dPsdzdz*Pe*Ph*Pf + Ps*dPedzdz*Ph*Pf + Ps*Pe*dPhdzdz*Pf + Ps*Pe*Ph*dPfdzdz + ...
        2*dPsdz*dPedz*Ph*Pf + 2*dPsdz*Pe*dPhdz*Pf + 2*dPsdz*Pe*Ph*dPfdz + ...
                             2*Ps*dPedz*dPhdz*Pf + 2*Ps*dPedz*Ph*dPfdz + ...
                                                  2*Ps*Pe*dPhdz*dPfdz + ...
        dPsdzdz;


    Hess = M*Hess; % Rest of year survival
    eigH = eig(Hess);
    eigH = real(eigH);
    res = max(eigH);

    resV = [resV;res];
end

% Plot it within Octave
if length(x_cV) > 1
    plot(x_cV,resV);
    xlabel('Optimal hatching time x_c')
    ylabel('Largest eigenvalue of Hessian')
    title('Singular strategy is ESS-stable where largest eig < 0')
end
