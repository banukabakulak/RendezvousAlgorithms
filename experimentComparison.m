function out = experimentComparison(nInst, baseParams)
%EXPERIMENTCOMPARISON  Compare C-MAS/D-MAS against three rendezvous baselines.
%
%   out = EXPERIMENTCOMPARISON(nInst, baseParams) runs the two-round scenario of
%   Sec.~\ref{sec:SoTAComparison} with N = 30 agents and reports energy E (mAh),
%   convergence time Tc (s) and communication overhead O (Kbits) for five methods:
%
%     C-MAS  centralized energy-aware rendezvous (CRM-OPT / aggregation)
%     D-MAS  distributed energy-aware rendezvous (DRMi-OPT + consensus)
%     CA     Consensus Averaging: the unweighted centroid of agent positions,
%            reached by neighbor averaging (no energy weighting)
%     RC     RoI Center: the geometric center of the region of interest
%     CC     Circumcenter: the minimax center (smallest enclosing circle)
%
%   Each method differs only in the rendezvous POINT it selects; that point is then
%   driven to by the SAME aggregation dynamics (SIMULATEROUND), so E and Tc isolate
%   the effect of the point-selection objective. Overhead uses each method's
%   communication model and packet size kappa (C-MAS 0.5, D-MAS 1.0, baselines 0.8):
%     C-MAS, RC, CC : O = kappa * step * N        (coordinator broadcast)
%     D-MAS, CA     : O = kappa * step * sum_i|N_i| (neighbor exchange)
%   with step the convergence time in control steps.
%
%   The denied region is disabled so every method's target is feasible and the
%   comparison reflects the rendezvous objective rather than denied-region handling;
%   the inter-round agent addition of the main scenario is retained.
%
%   Results are averaged over nInst instances; each method runs on the same per-seed
%   instance (GENERATEINSTANCE resets rng) for a fair comparison.
%
%   Defaults: nInst = 20, baseParams = defaultParams().
%
%   OUT has fields method, E, Tc, O, kappa (5x1 each) and OUT.dEpct (energy % over
%   C-MAS). A readable table and paste-ready LaTeX rows are printed.
%
%   See also RUNCMAS, RUNDMAS, SIMULATEROUND, COMPUTEMETRICS, COMPUTEMETRICSDMAS.

    arguments
        nInst      (1,1) double = 20
        baseParams (1,1) struct = defaultParams()
    end

    p = baseParams;
    p.sim.N                   = 30;
    p.sim.denied              = false;     % all targets feasible; isolate objective
    p.sim.forceDeniedCentroid = false;
    p.controller.logGraph     = false;     % joint-connectivity not needed here

    method = {'C-MAS','D-MAS','CA','RC','CC'};
    kappa  = [0.5; 1.0; 0.8; 0.8; 0.8];
    nM = numel(method);
    rawE = nan(nInst,nM);  rawT = nan(nInst,nM);  rawO = nan(nInst,nM);

    for s = 1:nInst
        seed = 9000 + s;

        % ---- C-MAS (reuse validated pipeline) ----
        p.sim.seed = seed;  inst = generateInstance(p);
        evalc('resC = runCMAS(p, inst); mC = computeMetrics(p, resC);');
        rawE(s,1) = mC.totalEnergy;  rawT(s,1) = mC.avgConsensusTime;  rawO(s,1) = mC.totalOverhead;

        % ---- D-MAS ----
        p.sim.seed = seed;  inst = generateInstance(p);
        evalc('resD = runDMAS(p, inst); mD = computeMetricsDMAS(p, resD);');
        rawE(s,2) = mD.totalEnergy;  rawT(s,2) = mD.avgConsensusTime;  rawO(s,2) = mD.totalOverhead;

        % ---- baselines CA / RC / CC ----
        for mi = 3:nM
            p.sim.seed = seed;  inst = generateInstance(p);
            rb = runBaseline(p, inst, method{mi});
            rawE(s,mi) = rb.E;  rawT(s,mi) = rb.T;  rawO(s,mi) = rb.O;
        end
    end

    out.method = method(:);
    out.E     = mean(rawE,1).';
    out.Tc    = mean(rawT,1).';
    out.O     = mean(rawO,1).';
    out.kappa = kappa;
    out.dEpct = 100 * (out.E - out.E(1)) / out.E(1);   % energy % over C-MAS
    out.raw   = struct('E',rawE,'Tc',rawT,'O',rawO);

    printTable(out);
    printLatex(out);
end

% =====================================================================
function rb = runBaseline(params, inst, method)
%RUNBASELINE  Two-round scenario driving the team to a baseline rendezvous point.
    nRounds = round(params.controller.simTime * params.controller.freqOuter);
    nStep   = round(params.controller.freqInner / params.controller.freqOuter);
    kap     = 0.8;
    E = 0;  O = 0;  times = nan(1,nRounds);

    for r = 1:nRounds
        if     r == 2, inst = applyRandomEvent(params, inst, 'add');
        elseif r == 3, inst = applyRandomEvent(params, inst, 'remove');
        end

        cg     = buildCommGraph(params, inst);         % round comm graph (pre-motion)
        sumDeg = sum(cg.deg);  Nr = size(inst.pos,1);

        rStar        = baselinePoint(inst, method);
        [traj, inst] = simulateRound(params, inst, rStar);
        cons         = analyzeConsensus(params, traj, rStar);

        step = cons.step;  if isnan(step), step = nStep; end
        switch method
            case 'CA',       ov = kap * step * sumDeg;  % neighbor exchange
            otherwise,       ov = kap * step * Nr;       % RC, CC coordinator broadcast
        end

        E = E + traj.roundEnergy;  O = O + ov;  times(r) = cons.time;
    end
    rb.E = E;  rb.T = mean(times,'omitnan');  rb.O = O;
end

% ---------------------------------------------------------------------
function r = baselinePoint(inst, method)
%BASELINEPOINT  Rendezvous point for the three baselines.
    P = inst.pos;  dim = size(P,2);
    switch method
        case 'CA'                                        % unweighted centroid
            r = mean(P, 1);
        case 'RC'                                        % RoI geometric center
            cxy = mean(inst.roi.vertices, 1);
            r = [cxy, meanZ(P, dim)];
            r = r(1:dim);
        case 'CC'                                        % minimax (smallest enclosing circle)
            c0 = mean(P(:,1:2), 1);
            f  = @(c) max(vecnorm(P(:,1:2) - c, 2, 2));
            c  = fminsearch(f, c0, optimset('Display','off','TolX',1e-4,'TolFun',1e-4));
            r  = [c, meanZ(P, dim)];
            r  = r(1:dim);
    end
end

function z = meanZ(P, dim)
    if dim >= 3, z = mean(P(:,3)); else, z = []; end
end

% ---------------------------------------------------------------------
function printTable(out)
    fprintf('\n=== Comparative performance (N=30, two-round, denied off) ===\n');
    fprintf('%-8s %-11s %-9s %-11s %-7s %-9s\n', 'Method','E(mAh)','Tc(s)','O(Kbits)','kappa','dE%%');
    for i = 1:numel(out.method)
        fprintf('%-8s %-11.2f %-9.2f %-11.1f %-7.1f %-9.2f\n', ...
            out.method{i}, out.E(i), out.Tc(i), out.O(i), out.kappa(i), out.dEpct(i));
    end
    fprintf('=============================================================\n');
    fprintf('Energy over C-MAS: CA %+.2f%%, RC %+.2f%%, CC %+.2f%%\n', ...
        out.dEpct(3), out.dEpct(4), out.dEpct(5));
end

% ---------------------------------------------------------------------
function printLatex(out)
    fprintf('\n%% --- LaTeX rows for tab:Comparison ---\n');
    for i = 1:numel(out.method)
        fprintf('\\texttt{%s} & $%.2f$ & $%.2f$ & $%.1f$ & %.1f \\\\\n', ...
            out.method{i}, out.E(i), out.Tc(i), out.O(i), out.kappa(i));
    end
    fprintf('%% -------------------------------------\n');
end
