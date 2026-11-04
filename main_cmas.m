%% main_cmas.m  -  Energy-efficient rendezvous for heterogeneous MAS (C-MAS)
%
%  Driver script: configure the scenario, generate a heterogeneous rendezvous
%  instance, inspect the communication graph, visualize the instance, and run
%  the centralized (C-MAS) controller. Later stages (energy and communication-
%  overhead accounting, performance plots) are added as separate functions.
%
%  Usage:
%      >> main_cmas
%
%  Repository: <GitHub URL>
%  Reference : Energy-Efficient Rendezvous Algorithms for Heterogeneous
%              Multi-Agent Systems, IEEE T-SMC: Systems.

clear; clc;

% 1) Configuration -------------------------------------------------------
params = defaultParams();

% Example overrides (uncomment to change the scenario):
% params.sim.N                   = 50;     % more agents
% params.arena.width             = 28;     % smaller arena -> denser graph
params.sim.seed                = 6;      % Sec. VI-B demonstration instance
% params.sim.forceDeniedCentroid = true;   % guarantee H-CRM-D activation

% 2) Instance generation -------------------------------------------------
inst = generateInstance(params);
printInstanceSummary(params, inst);

% 3) Communication graph and connectivity --------------------------------
commGraph = buildCommGraph(params, inst);
metrics   = graphConnectivity(commGraph);
printGraphMetrics(metrics, params.sim.N);

% 4) Visualize the instance ----------------------------------------------
% h = plotInstance(params, inst);

% 5) Run the C-MAS controller --------------------------------------------
results = runCMAS(params, inst);

% 6) Performance metrics (energy, overhead, messages, duty) ---------------
perf = computeMetrics(params, results);

% 7) Performance plots (pose convergence + cumulative energy) -------------
hPerf = plotPerformance(params, results); 
% plotConnectivity(params, results);

% 8) Trajectory visualization and animation ------------------------------
% hTraj = plotTrajectories(params, results);
% animateRendezvous(params, results);                       % interactive
% animateRendezvous(params, results, 'Record', true, 'Stride', 5);  % save video
% animateRendezvous(params, results, 'Record', true);
