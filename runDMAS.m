function results = runDMAS(params, inst)
%RUNDMAS  Run the distributed (D-MAS) rendezvous controller.
%
%   results = RUNDMAS(params, inst) executes the multi-round D-MAS scenario.
%   It mirrors RUNCMAS but replaces the centralized rendezvous solve with the
%   distributed DRMi-D consensus: for each outer round it (optionally) applies
%   an inter-round random event, rebuilds the communication graph, computes the
%   rendezvous point by local DRMi-OPT + power-weighted consensus (+ H-DRMi-D),
%   drives the team there with the shared aggregation controller, and records
%   convergence and communication statistics.
%
%   Returns a struct:
%     results.rendezvous  nRounds x dim        rendezvous point per round
%     results.traj        1 x nRounds cell     trajectories (traj.P, traj.t, ...)
%     results.consensus   1 x nRounds struct   AnalyzeConsensus output per round
%     results.info        1 x nRounds struct   DRMi-D info (consIters, gap, ...)
%     results.comm        1 x nRounds struct   graph stats (sumDeg, avgDeg,
%                                               connected, lambda2)
%     results.conn        1 x nRounds struct   joint-connectivity check per
%                                               round (minLambda2, allConnected,
%                                               window, windowLambda2)
%     results.instFinal   final instance
%
%   The per-round sumDeg = sum_i |N_i| is what the D-MAS communication-overhead
%   and message metrics integrate over the convergence time (see COMPUTEMETRICSDMAS).
%
%   See also SOLVERENDEZVOUSDISTRIBUTED, BUILDCOMMGRAPH, SIMULATEROUND,
%   ANALYZECONSENSUS, COMPUTEMETRICSDMAS, RUNCMAS.

arguments
    params (1,1) struct
    inst   (1,1) struct
end

nRounds = round(params.controller.simTime * params.controller.freqOuter);
dim     = params.sim.dim;

results.rendezvous = zeros(nRounds, dim);
results.traj       = cell(1, nRounds);
results.consensus  = struct('step', {}, 'time', {}, 'duty', {});
results.info       = struct('proposals', {}, 'consIters', {}, 'gap', {}, ...
    'deniedHit', {}, 'offsets', {}, 'agreed', {});
results.comm       = struct('sumDeg', {}, 'avgDeg', {}, 'connected', {}, 'lambda2', {});
results.conn       = struct('minLambda2', {}, 'allConnected', {}, 'window', {}, 'windowLambda2', {});

for r = 1:nRounds
    % inter-round events: add an agent at round 2, remove one at round 3
    if r == 2
        inst = applyRandomEvent(params, inst, 'add');
    elseif r == 3
        inst = applyRandomEvent(params, inst, 'remove');
    end

    commGraph           = buildCommGraph(params, inst);
    gm                  = graphConnectivity(commGraph);
    [rStar, info]       = solveRendezvousDistributed(params, inst, commGraph);
    [traj, inst]        = simulateRound(params, inst, rStar);
    cons                = analyzeConsensus(params, traj, rStar);
    jc                  = jointConnectivity(params, traj);

    traj.id   = inst.id;             % persistent ids present this round
    traj.type = inst.agentType;      % their types

    results.rendezvous(r, :) = rStar;
    results.traj{r}          = traj;
    results.consensus(r)     = cons;
    results.info(r)          = info;
    results.comm(r).sumDeg    = sum(commGraph.deg);
    results.comm(r).avgDeg    = mean(commGraph.deg);
    results.comm(r).connected = gm.connected;
    results.comm(r).lambda2   = gm.lambda2;
    results.conn(r).minLambda2    = jc.minLambda2;
    results.conn(r).allConnected  = jc.allConnected;
    results.conn(r).window        = jc.window;
    results.conn(r).windowLambda2 = jc.windowLambda2;

    fprintf(['Round %d: N=%d  rStar=(%s)  cons-iters=%d  gap=%.3f m  ' ...
        'denied-hit=%d  joint-conn=%d (min lambda2=%.3f)  Tc=%.1f s\n'], ...
        r, numel(inst.agentType), num2str(rStar, '%.2f '), ...
        info.consIters, info.gap, info.deniedHit, jc.allConnected, ...
        jc.minLambda2, cons.time);
end

results.instFinal = inst;
end
