function [resV,Jac] = numcheck_conv(p,yzV,x_cV)

% -- [resV,Jac] = numcheck_conv(p,yzV,x_cV)
%
% The purpose of this function is to check the convergence
% stability of each evolutionarily singular strategy defined
% by the inputs. The check is done by finding the largest
% eigenvalue of the Jacobian matrix and returning it to the
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
% resV: A vector of largest eigenvalues of the Jacobian
% matrix corresponding to the evolutionarily singular
% strategies specified by the inputs.
%
% Jac: The Jacobian matrix of the last evolutionarily
% singular strategy specified. It is useful if you want to
% take a closer look at just one singular strategy.
%
%
% see also: check_conv, numcheck_ess
%

resV = []; % Storage vector for output
del = 1e-4; % Small number for estimating derivative

% This loop steps through each of our x_c values
for ind = 1:length(x_cV);

    YZ = yzV(ind,:);
    x_c = x_cV(ind);

    Y = YZ(1); Z = YZ(2);

    % Little step over which to estimate derivative
    Yhi = Y+del/2; Ylo = Y-del/2;
    Zhi = Z+del/2; Zlo = Z-del/2;


    % Estimate Jacobian numerically by estimating the
    % derivatives of the fitness gradient returned by dbo
    g_Yhi = dbo(p,[Yhi,Z],x_c,[Yhi,Z]);
    g_Ylo = dbo(p,[Ylo,Z],x_c,[Ylo,Z]);
    g_dY = (g_Yhi-g_Ylo)/del;

    g_Zhi = dbo(p,[Y,Zhi],x_c,[Y,Zhi]);
    g_Zlo = dbo(p,[Y,Zlo],x_c,[Y,Zlo]);
    g_dZ = (g_Zhi-g_Zlo)/del;

    Jac = [g_dY,g_dZ];
    eigJ = max(real(eig(Jac)));

    resV = [resV;eigJ];
end

% Plot it within Octave
if length(x_cV) > 1
    plot(x_cV,resV);
    xlabel('Optimal hatching time x_c')
    ylabel('Largest eigenvalue of Jacobian')
    title('Singular strategy is convergence stable where largest eig < 0')
end
