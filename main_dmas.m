%% main_dmas.m  -  Energy-efficient rendezvous for heterogeneous MAS (D-MAS)
%
%  Driver script for the distributed controller: configure the scenario,
%  generate a heterogeneous rendezvous instance, inspect the communication
%  graph, visualize the instance, and run the distributed (D-MAS) controller
%  with DRMi-D consensus rendezvous and H-DRMi-D denied-region handling.
%
%  The instance, communication graph, aggregation dynamics, consensus analysis
%  and all plots are shared with the centralized study (see MAIN_CMAS); only the
%  rendezvous computation (SOLVERENDEZVOUSDISTRIBUTED, called inside RUNDMAS) and
%  the communication-cost metrics (COMPUTEMETRICSDMAS) differ.
%
%  Usage:
%      >> main_dmas
%
%  Repository: <GitHub URL>
%  Reference : Energy-Efficient Rendezvous Algorithms for Heterogeneous
%              Multi-Agent Systems, IEEE T-SMC: Systems.

clear; clc;

% 1) Configuration -------------------------------------------------------
params = defaultParams();

% Example overrides (uncomment to change the scenario):
% params.sim.N                   = 50;     % more agents
params.sim.seed                = 6;      % Sec. VI-B demonstration instance
% params.sim.forceDeniedCentroid = true;   % guarantee H-DRMi-D activation

% 2) Instance generation -------------------------------------------------
inst = generateInstance(params);
printInstanceSummary(params, inst);

% 3) Communication graph and connectivity --------------------------------
commGraph = buildCommGraph(params, inst);
metrics   = graphConnectivity(commGraph);
printGraphMetrics(metrics, params.sim.N);

% 4) Visualize the instance ----------------------------------------------
% h = plotInstance(params, inst);

% 5) Run the D-MAS controller --------------------------------------------
results = runDMAS(params, inst);

% 6) Performance metrics (energy, D-MAS overhead, messages, consensus) ----
perf = computeMetricsDMAS(params, results);

% 7) Performance plots (pose convergence + cumulative energy) -------------
% hPerf = plotPerformance(params, results);
% plotConnectivity(params, results);

% 8) Trajectory visualization and animation ------------------------------
% hTraj = plotTrajectories(params, results);
% animateRendezvous(params, results);                       % interactive
% animateRendezvous(params, results, 'Record', true, 'Stride', 5);  % save video
% animateRendezvous(params, results, 'Record', true);
