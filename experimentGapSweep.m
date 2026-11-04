function out = experimentGapSweep(Lvalues, nInst, doEnergy, baseParams)
%EXPERIMENTGAPSWEEP  Centralized--distributed equivalence vs topology sparsity.
%
%   out = EXPERIMENTGAPSWEEP(Lvalues, nInst, doEnergy, baseParams) sweeps the
%   arena side L (with N and the communication ranges fixed) so the initial
%   communication graph ranges from near-complete (small L) to sparse multi-hop
%   (large L), and quantifies how far the distributed (D-MAS) rendezvous point
%   departs from the centralized (C-MAS) one as the graph sparsifies.
%
%   Only CONNECTED instances are measured: for each draw the instance is
%   regenerated (new seed) until its initial communication graph is a single
%   component, so every row satisfies the connectivity premise the section rests
%   on (disconnected graphs, on which the static-graph consensus cannot agree,
%   are excluded rather than averaged in).
%
%   For each L it keeps nInst connected instances and, on each, compares the
%   C-MAS and D-MAS rendezvous points on the SAME instance and records:
%
%     diam       initial communication-graph diameter (hops)
%     meanDeg    initial mean degree
%     ptDev      normalized point deviation, % of arena side
%                = 100 * || r_D - r_C || / L   (scale-invariant; this is the
%                centralized-distributed deviation bounded by Theorem S2)
%     ptGap      || r_D - r_C ||, m  (raw point deviation)
%     objGap     CRM-objective gap of the D-MAS point vs the C-MAS point, %
%                (secondary; can be either sign because both points are
%                feasibility heuristics, not the true constrained optimum)
%     rounds     D-MAS consensus iterations to reach agreement
%     consMsg    rendezvous-consensus messages = rounds * sum_i |N_i|
%                (the cost of AGREEING on r*, independent of the aggregation
%                phase and of T_c; distinct from the aggregation-phase message
%                count Gamma in the main performance table)
%     energyGap  (only if doEnergy) D-MAS/C-MAS single-round mission-energy gap, %
%
%   The regions (RoI, denied) are rescaled with L so the geometry is preserved
%   and only the graph density changes. Each metric is averaged over the kept
%   instances; raw per-instance values are in OUT.raw, and OUT.nUsed reports how
%   many connected instances were actually found per L.
%
%   Defaults: Lvalues = [12 16 20 24 28 32] m, nInst = 20, doEnergy = false,
%   baseParams = defaultParams().
%
%   See also SOLVERENDEZVOUSCENTRAL, SOLVERENDEZVOUSDISTRIBUTED, PLOTGAPSWEEP.

    arguments
        Lvalues    (1,:) double  = [12 16 20 24 28 32]
        nInst      (1,1) double  = 20
        doEnergy   (1,1) logical = false
        baseParams (1,1) struct  = defaultParams()
    end

    baseParams.controller.logGraph = false;                               % not needed here
    baseParams.controller.simTime  = 1 / baseParams.controller.freqOuter;  % single round

    % This section studies the CONSENSUS equivalence (CRM <-> DRMi), so the
    % denied region is disabled: it isolates the topology effect from the
    % H-CRM-D/H-DRMi-D offset, which otherwise adds a density-independent
    % deviation that dominates the gap even on dense graphs. With the denied
    % region off, r_C is the true RoI-constrained optimum, so objGap >= 0, is
    % zero on a complete graph, and grows monotonically with sparsity.
    baseParams.sim.denied             = false;
    baseParams.sim.forceDeniedCentroid = false;
    roi0 = baseParams.region.roi;  denied0 = baseParams.region.denied;  L0 = baseParams.arena.width;
    N    = baseParams.sim.N;
    maxTries = 300;                                      % attempts to draw a connected instance

    nL = numel(Lvalues);
    z  = nan(nL,1);
    out.L = Lvalues(:);  out.diam = z;  out.meanDeg = z;  out.nUsed = z;
    out.ptDev = z;  out.ptGap = z;  out.objGap = z;  out.rounds = z;  out.consMsg = z;  out.energyGap = z;
    out.raw = struct('L',{},'diam',{},'meanDeg',{},'ptDev',{},'ptGap',{}, ...
                     'objGap',{},'rounds',{},'consMsg',{},'energyGap',{});

    for li = 1:nL
        L  = Lvalues(li);
        sc = L / L0;

        d = nan(nInst,1); md = nan(nInst,1); pd = nan(nInst,1); pg = nan(nInst,1);
        og = nan(nInst,1); rr = nan(nInst,1); cm = nan(nInst,1); eg = nan(nInst,1);
        kept = 0;

        for s = 1:nInst
            % ---- draw a CONNECTED instance ----
            cg = []; g0 = []; inst = []; p = [];
            ok = false;
            for tr = 1:maxTries
                p = baseParams;
                p.arena.width  = L;
                p.region.roi    = roi0    * sc;          % preserve geometry, scale size
                p.region.denied = denied0 * sc;
                p.sim.seed = 100000*li + 1000*s + tr;    % distinct, reproducible
                inst = generateInstance(p);
                cg   = buildCommGraph(p, inst);
                g0   = graphConnectivity(cg);
                if g0.connected, ok = true; break; end
            end
            if ~ok, continue; end                        % skip if no connected draw found
            kept = kept + 1;

            sumDeg = sum(cg.deg);                         % sum_i |N_i| (undirected)
            [rD, infoD] = solveRendezvousDistributed(p, inst, cg);
            [objGap, ptGap] = objectiveGap(p, inst, rD);

            d(s)  = g0.diameter;  md(s) = g0.meanDegree;
            pg(s) = ptGap;        pd(s) = 100 * ptGap / L;     % normalized by arena side
            og(s) = objGap;
            rr(s) = infoD.consIters;
            cm(s) = infoD.consIters * sumDeg;

            if doEnergy
                s0 = rng;                                 % identical RNG for both controllers
                [~, resC] = evalc('runCMAS(p, inst)');    %#ok<ASGLU>
                rng(s0);
                [~, resD] = evalc('runDMAS(p, inst)');    %#ok<ASGLU>
                eg(s) = 100*(resD.traj{1}.roundEnergy - resC.traj{1}.roundEnergy) ...
                              / resC.traj{1}.roundEnergy;
            end

            out.raw(end+1) = struct('L',L,'diam',d(s),'meanDeg',md(s),'ptDev',pd(s), ...
                'ptGap',pg(s),'objGap',og(s),'rounds',rr(s),'consMsg',cm(s), ...
                'energyGap',eg(s)); %#ok<AGROW>
        end

        out.diam(li)     = mean(d,  'omitnan');
        out.meanDeg(li)  = mean(md, 'omitnan');
        out.nUsed(li)    = kept;
        out.ptDev(li)    = mean(pd, 'omitnan');
        out.ptGap(li)    = mean(pg, 'omitnan');
        out.objGap(li)   = mean(og, 'omitnan');
        out.rounds(li)   = mean(rr, 'omitnan');
        out.consMsg(li)  = mean(cm, 'omitnan');
        out.energyGap(li)= mean(eg, 'omitnan');

        fprintf('L=%2d m (%2d/%2d connected): deg=%.1f diam=%.1f | ptDev=%.2f%% ptGap=%.2f m objGap=%.2f%% rounds=%.0f consMsg=%.0f', ...
            L, kept, nInst, out.meanDeg(li), out.diam(li), out.ptDev(li), ...
            out.ptGap(li), out.objGap(li), out.rounds(li), out.consMsg(li));
        if doEnergy, fprintf(' Egap=%.2f%%', out.energyGap(li)); end
        fprintf('\n');
    end

    printTable(out, nInst, N, doEnergy);
end

% =====================================================================
function [objGapPct, ptGap] = objectiveGap(params, inst, rD)
%OBJECTIVEGAP  CRM-objective gap of the D-MAS point RD vs the C-MAS point.
%   f(r) = sum_i c_i ||p_i - r||^2, c_i = w_i(v_i)/||v_i||. Both points are the
%   controllers' final (post-offset) rendezvous points, so the gap is zero on a
%   complete graph and grows as the topology sparsifies.
    speed = vecnorm(inst.vel, 2, 2);
    c     = (inst.cap.powerA .* speed.^2 + inst.cap.powerB) ./ speed;

    rC = solveRendezvousCentral(params, inst);          % C-MAS point (H-CRM-D)

    f = @(r) sum(c .* sum((inst.pos - r).^2, 2));
    fC = f(rC);
    objGapPct = 100 * (f(rD) - fC) / max(fC, eps);
    ptGap     = norm(rD - rC);
end

% =====================================================================
function printTable(out, nInst, N, doEnergy)
    fprintf('\n=== D-MAS vs C-MAS: equivalence gap and consensus cost ===\n');
    fprintf('(N=%d, up to %d connected instances per row)\n', N, nInst);
    if doEnergy
        fprintf('%-6s %-6s %-8s %-7s %-9s %-9s %-9s %-8s %-9s\n', ...
            'L(m)','diam','meanDeg','nUsed','ptDev(%)','objGap(%)','Egap(%)','rounds','consMsg');
        for li = 1:numel(out.L)
            fprintf('%-6d %-6.1f %-8.1f %-7.0f %-9.2f %-9.2f %-9.2f %-8.0f %-9.0f\n', ...
                out.L(li), out.diam(li), out.meanDeg(li), out.nUsed(li), ...
                out.ptDev(li), out.objGap(li), out.energyGap(li), out.rounds(li), out.consMsg(li));
        end
    else
        fprintf('%-6s %-6s %-8s %-7s %-9s %-9s %-8s %-9s\n', ...
            'L(m)','diam','meanDeg','nUsed','ptDev(%)','objGap(%)','rounds','consMsg');
        for li = 1:numel(out.L)
            fprintf('%-6d %-6.1f %-8.1f %-7.0f %-9.2f %-9.2f %-8.0f %-9.0f\n', ...
                out.L(li), out.diam(li), out.meanDeg(li), out.nUsed(li), ...
                out.ptDev(li), out.objGap(li), out.rounds(li), out.consMsg(li));
        end
    end
    fprintf('==========================================================\n');
end
