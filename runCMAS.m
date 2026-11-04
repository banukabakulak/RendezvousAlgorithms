function results = runCMAS(params, inst)
%RUNCMAS  Run the centralized (C-MAS) rendezvous controller.
%
%   results = RUNCMAS(params, inst) executes the multi-round C-MAS scenario:
%   for each outer round it (optionally) applies an inter-round random event,
%   computes the rendezvous point with CRM-OPT / H-CRM-D, drives the team there
%   with the aggregation controller, and records convergence statistics.
%
%   Returns a struct:
%     results.rendezvous  nRounds x dim   rendezvous point per round
%     results.traj        1 x nRounds cell of trajectories (traj.P, traj.t)
%     results.consensus   1 x nRounds struct array (step, time, duty)
%     results.info        1 x nRounds struct array (CRM-OPT / H-CRM-D info)
%     results.conn        1 x nRounds struct array joint-connectivity check
%                         (minLambda2, allConnected, window, windowLambda2)
%     results.instFinal   final instance (positions, velocities, agent count)
%
%   See also SOLVERENDEZVOUSCENTRAL, SIMULATEROUND, ANALYZECONSENSUS,
%   JOINTCONNECTIVITY, APPLYRANDOMEVENT.

arguments
    params (1,1) struct
    inst   (1,1) struct
end

nRounds = round(params.controller.simTime * params.controller.freqOuter);
dim     = params.sim.dim;

results.rendezvous = zeros(nRounds, dim);
results.traj       = cell(1, nRounds);
results.consensus  = struct('step', {}, 'time', {}, 'duty', {});
results.info       = struct('lambda', {}, 'iters', {}, 'deniedHit', {}, 'offsets', {});
results.conn       = struct('minLambda2', {}, 'allConnected', {}, 'window', {}, 'windowLambda2', {});

for r = 1:nRounds
    % inter-round events: add an agent at round 2, remove one at round 3
    if r == 2
        inst = applyRandomEvent(params, inst, 'add');
    elseif r == 3
        inst = applyRandomEvent(params, inst, 'remove');
    end

    [rStar, info]        = solveRendezvousCentral(params, inst);
    [traj, inst]         = simulateRound(params, inst, rStar);
    cons                 = analyzeConsensus(params, traj, rStar);
    jc                   = jointConnectivity(params, traj);

    traj.id   = inst.id;             % persistent ids present this round
    traj.type = inst.agentType;      % their types

    results.rendezvous(r, :) = rStar;
    results.traj{r}          = traj;
    results.consensus(r)     = cons;
    results.info(r)          = info;
    results.conn(r).minLambda2    = jc.minLambda2;
    results.conn(r).allConnected  = jc.allConnected;
    results.conn(r).window        = jc.window;
    results.conn(r).windowLambda2 = jc.windowLambda2;

    fprintf(['Round %d: N=%d  rStar=(%s)  denied-hit=%d (%d offsets)  ' ...
        'joint-conn=%d (min lambda2=%.3f)  Tc=%.1f s\n'], ...
        r, numel(inst.agentType), num2str(rStar, '%.2f '), ...
        info.deniedHit, info.offsets, jc.allConnected, jc.minLambda2, cons.time);
end

results.instFinal = inst;
end
