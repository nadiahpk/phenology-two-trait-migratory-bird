% Pe
p.K = 100; % Doesn't matter, scales n
p.b_e = 0.03; % Range of 0 to 0.125 in Johansson (2012; Ecology Letters), chosen so arrival date ~ 130 when Pq ~ 1

% Pf
p.a=3; % Half of the average clutch size, taken from Askenmo (1981; Ardea);
p.sigma = 5; % Jonzen et al. (2007)

% Ph
p.Q_0 = 0.05; % Initial egg-resource? Only thing we know is that it's an income breeder, so this should be 'low'
p.b_q = 0.5; % Slope of resource increase? Depends on the resource.  Jonzen used a time of 10 days for the minimum prelaying  period, this combination will get to 0.9 in about 5 days.
p.u_q = 134; % Flat arrival time, choose as midpoint for graphs
%p.u_q = 125; % Good looking graph, little upkick at end, looks like flycatchers.

% Ps
p.L_m = 170; % Jonzen et al. (2007)
p.b_s = 0.075; % 0.05 > bs > 0.1 was the range investigated by Jonzen. Choose the middle point.

p.M = 0.5; % Sanz (2001,The Auk) puts females at 45 - 52%

% Key time
p.z_f = 40; % Jonzen et al.'s time to end of breeding minus time to laying
p.z_n = 14; % Incubation time

