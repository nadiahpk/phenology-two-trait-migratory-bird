function [x_cV,yzV,nV] = solvedbo(p,change_p,change_p_val,x_c_hi,x_c_lo,y0,z0,print_flag);

% -- [x_cV,yzV,nV] = solvedbo()
% -- [x_cV,yzV,nV] = solvedbo(p,change_p,change_p_val,x_c_hi,x_c_lo,y0,z0,print_flag);
%
%
% The purpose of this function is to plot evolutionarily
% singular bird phenology with respect to optimal hatching
% time. When run without any parameters, it produces the
% left-most pane of Figure 1. Otherwise parameters can be
% specified to explore how Figure 1 changes with changes in
% a parameter value.
%
% Note that the components were all taken as f(y,z), where 
% z is the prelaying period. In the general model in the
% paper we use f(x,y). The choice can be made on 
% convenience, so long as one is consistent throughout a
% given analysis.
%
%
% INPUTS
%
% p: A dictionary of parameter values. Use the file params.m
% for an example of how to specify one of these.
%
% change_p: A string specifying the key of one of the
% parameter values defined in p.
%
% change_p_val: The value to which change_p should be set.
%
% x_c_hi: The upper value of the optimal hatching time range
% over which the figure should be generated.
%
% x_c_lo: The low value of the optimal hatching time range
% over which the figure should be generated.
%
% y0: The starting guess for what the evolutionarily
% singular arrival time is at x_c = x_c_hi.
% Note that if the solver is failing, it's probably because
% either this value or z0 have not been chosen close enough
% to the true values.
%
% z0: The starting guess for what the evolutionarily
% singular prelaying period is at x_c = x_c_hi. If x_c_hi is
% in the 'flat' region for y, it's a good idea to set this
% such that x is close to the optimal hatching time. 
% Note that if the solver is failing, it's probably because
% either this value or y0 have not been chosen close enough
% to the true values.
%
% print_flag: Determines if the return values will be
% written to a .dat file or not. Set this to 1 if you'd like
% that or any other value if not.
%
%
% OUTPUTS
% 
% x_cV: A vector of x_c values in the range specified. 
%
% yzV: A matrix of [y*,z*] values in the range specified,
% corresponding to x_cV for easy plotting.
%
% nV: A vector of population sizes corresponding to x_cV.
%

if nargin == 0;
    % User would like default, which I've set to the
    % left-most pane of Figure 1

    params; % This gets our parameter-value dictionary p

    % Change parameter values to see what effect it has on
    % the way in which the ESS responds to advancing food peak
    p.u_q = 125
    change_p = 'u_q';

    % This is the advancing food peak. We start at a high value (late
    % peak) and move it forward in time to simulate late-spring warming
    % due to climate change. At each x_c value we solve the trait ESS, and
    % then plot it like in Figure 2
    x_cV = linspace(175,135,50); % Should be enough to get the whole range at reasonable resolution

    % Initial guess for trait values
    y0 = 120; % Arrival time
    z0 = x_cV(1)-y0-p.z_n; % Prelaying period. Our initial guess is the optimum, which works well when we start on the rhs of Figure 2

    print_flag = 1; % Print the output to the .dat file so it can be plotted with Gnuplot

else

    % Modify parameter value 
    % For example, if change_p is the string 'u_q' and
    % change_p_val = 125, this evaluates 
    % > p.u_q = 125;
    str = ['p.',change_p,' = ',num2str(change_p_val),';'];
    eval(str);

    x_cV = linspace(x_c_hi,x_c_lo,50); % Should be enough to get the whole range at reasonable resolution
    % Note we go backwards through the x_c values because
    % that helps the solver
end

% Prepare storage vectors
yzV = []; 
nV = [];

% Initialise counters and flags
count = 1;
viable_flag = 1; % True while we have a viable steady state population

% This loop steps through each x_c value and finds the trait ESS
while (viable_flag == 1) && (count <= length(x_cV));

    x_c = x_cV(count);
    yz0 = [y0,z0];

    % Solve for the evolutionarily singular strategy

    % The functiond dbo defines our fitness gradient
    f = @(v) dbo(p,v,x_c); 
    % An evolutionarily singular strategy is when that
    % gradient is zero
    [yz,fval] = fsolve(f,yz0); 

    % At this strategy, we'll need a viable population n >= K
    n = calcn(p,yz,x_c);
    
    % If we're returning reasonable trait values and our population is
    % still viable
    if (min(yz) > 0) && (n >= p.K)

        % Store our result
        yzV = [yzV;yz];
        nV = [nV;n/p.K];

        % Update - use this step's value as the next step's first
        % guess, as this makes it much more likely to find the solution
        z0 = yz(2);
        y0 = yz(1);

        count += 1;
    else

        % Don't bother storing this result.
        x_cV(count:end) = []; % Omit this x_c value onwards so we can plot easily.
        viable_flag = 0; % Update flag to end the while-loop

    end
end

% Turn this on to print results to a file
if print_flag == 1
    pval = eval(['p.',change_p]);
    str = [change_p,num2str(pval),'.dat'];
    fid = fopen(str,'w');
    str = ['# This is a run of ',change_p,' = ',num2str(pval),'\n'];
    fprintf(fid,str);
    fprintf(fid,'# It was generated by solvedbo.m \n');
    fprintf(fid,'# Here are the parameter values that were used \n');
    names = fieldnames(p);
    for ind = 1:length(names)
        entry = names(ind);
        fprintf(fid,'# ');
        fprintf(fid,cell2mat(entry));
        fprintf(fid,' = ');
        fprintf(fid,num2str(p.(cell2mat(entry))));
        fprintf(fid,'; ');
        str = ['p.',cell2mat(entry),' = ',cell2mat(entry),';\n'];
        fprintf(fid,str);
    end

    fprintf(fid,'# Columns are:\n');
    fprintf(fid,'# x_c \t \t \t y \t \t \t z \t \t \t n \n');
    fprintf(fid,'%1.4f \t %1.4f \t %1.4f \t %1.4f \n',[x_cV',yzV,nV]');
    fclose(fid);
end

% Plot the results within Octave
plot(x_cV,yzV(:,1),'k;Arrival date;')
hold on
plot(x_cV,yzV(:,1)+yzV(:,2),'r;Laying date;')
plot(x_cV,x_cV,'k.') 
plot(x_cV,yzV(:,1)+yzV(:,2)+p.z_n,'b;Hatching date;')
hold off
xlabel('Optimal hatching time x_c')
ylabel('Bird phenology')
axis equal

