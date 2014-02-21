function [resV,Hess] = numcheck_ess(p,yzV,x_cV)

% -- [resV,Hess] = numcheck_ess(p,yzV,x_cV)
%
% The purpose of this function is to check the ESS
% stability of each evolutionarily singular strategy defined
% by the inputs. The check is done by finding the largest
% eigenvalue of the Hessian matrix and returning it to the
% user, who can verify that they are zero.
% 
% This code uses a numerical approach in
% that the expressions for each of the derivatives and
% second derivatives are found by estimating the derivatives
% of the fitness gradient returned by dbo.m. 
% This code also has a counterpart check_conv.m which
% finds the eigenvalues numerically. Agreement between the
% two methods was used a check of the correctness of the
% solutions.
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
% see also: check_ess, numcheck_conv
%

resV = []; % Storage vector for output
del = 1e-4; % Small number for estimating derivative

% This loop steps through each of our x_c values
for ind = 1:length(x_cV);

    YZ = yzV(ind,:);
    x_c = x_cV(ind);

    Y = YZ(1); Z = YZ(2);
    y = Y; z = Z;

    % Little step over which to estimate derivative
    yhi = y+del/2; ylo = y-del/2;
    zhi = z+del/2; zlo = z-del/2;

    % Estimate Hessian numerically by estimating the
    % derivatives of the fitness gradient returned by dbo
    g_yhi = dbo(p,[yhi,z],x_c,YZ);
    g_ylo = dbo(p,[ylo,z],x_c,YZ);
    g_dy = (g_yhi-g_ylo)/del;

    g_zhi = dbo(p,[y,zhi],x_c,YZ);
    g_zlo = dbo(p,[y,zlo],x_c,YZ);
    g_dz = (g_zhi-g_zlo)/del;

    Hess = [g_dy,g_dz];
    eigH = max(real(eig(Hess)));

    resV = [resV;eigH];
end

% Plot it within Octave
if length(x_cV) > 1
    plot(x_cV,resV);
    xlabel('Optimal hatching time x_c')
    ylabel('Largest eigenvalue of Hessian')
    title('Singular strategy is convergence stable where largest eig < 0')
end
