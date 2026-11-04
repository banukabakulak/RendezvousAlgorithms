function [rStar, info] = solveRendezvousDistributed(params, inst, commGraph)
%SOLVERENDEZVOUSDISTRIBUTED  Distributed rendezvous point (DRMi-D + H-DRMi-D).
%
%   [rStar, info] = SOLVERENDEZVOUSDISTRIBUTED(params, inst, commGraph) computes
%   the rendezvous point by local computation and consensus, rather than by a
%   single centralized solve:
%
%     1. DRMi-OPT. Each agent i forms a local, feasible rendezvous estimate by
%        solving CRM-OPT over its own neighbourhood N_i (itself plus its mutual
%        communication neighbours), using SOLVERENDEZVOUSCENTRAL on the agent
%        subset. The local estimate is already offset out of the denied region
%        (H-CRM-D), so every proposal is individually feasible.
%
%     2. Consensus. Agents run the self-normalizing power-weighted average
%        consensus update (Eq. DRM_i_Consensus)
%           r_i <- r_i + sum_{j in N_i} w_j / ( w_i + sum_{k in N_i} w_k ) (r_j - r_i),
%        over the undirected (mutual) communication graph COMMGRAPH.A, until the
%        estimates agree (disagreement spread < agreeTol) or MAXITER is reached.
%        The update is row-stochastic with positive self-weight, so on a
%        connected graph it converges geometrically to a common point; on a
%        complete graph that point is the centralized optimum.
%
%     3. H-DRMi-D. Because R \ D is non-convex, the agreed point can itself fall
%        in the denied region even though every proposal was feasible; a final
%        ESCAPEDENIED offset moves it out (check--offset).
%
%   Returns rStar (1 x dim) and info with fields:
%     proposals  N x dim  per-agent local DRMi-OPT estimates (consensus seeds)
%     consIters  consensus iterations to agreement
%     gap        || agreed - centralized CRM-OPT point ||  (topology deviation)
%     deniedHit  true if the agreed point was offset out of the denied region
%     offsets    number of offset steps applied (0 or 1)
%     agreed     1 x dim  consensus point before the H-DRMi-D offset
%
%   See also SOLVERENDEZVOUSCENTRAL, BUILDCOMMGRAPH, ESCAPEDENIED, RUNDMAS.

    arguments
        params    (1,1) struct
        inst      (1,1) struct
        commGraph (1,1) struct
    end

    dim = params.sim.dim;
    N   = size(inst.pos, 1);
    A   = commGraph.A;            % undirected (mutual) adjacency
    w   = inst.power;             % consensus mixing weights w_i(v_i) = powerA v^2 + powerB

    % ---- 1) local DRMi-OPT: each agent's feasible estimate over N_i ----
    R = zeros(N, dim);
    for i = 1:N
        idx    = [i, find(A(i, :))];              % self + mutual neighbours
        sub    = subInstance(inst, idx);
        R(i,:) = solveRendezvousCentral(params, sub);
    end
    proposals = R;

    % ---- 2) power-weighted average consensus (Eq. DRM_i_Consensus) ----
    % Each agent moves toward the power-weighted mean of its neighbours'
    % estimates with a self-normalizing, fixed gain:
    %   r_i <- r_i + sum_{j in N_i}  w_j / ( w_i + sum_{k in N_i} w_k )  (r_j - r_i).
    % The induced update matrix is row-stochastic with a strictly positive
    % self-weight w_i / ( w_i + sum_k w_k ), so on a (jointly) connected graph
    % the estimates converge geometrically to a common agreement point.
    % Stopping is on the disagreement SPREAD (distance of the estimates from
    % their common mean), the true measure of agreement.
    maxIter   = 1000;
    agreeTol  = 1e-2;                             % agreement tolerance, m
    consIters = maxIter;
    for t = 1:maxIter
        Rnew = R;
        for i = 1:N
            nb = find(A(i, :));
            if isempty(nb), continue; end
            wn    = reshape(w(nb), 1, []);        % 1 x |N_i| neighbour weights
            denom = w(i) + sum(wn);               % w_i + sum_{k in N_i} w_k
            if denom <= 0, continue; end
            delta      = (wn * (R(nb, :) - R(i, :))) / denom;  % 1 x dim
            Rnew(i, :) = R(i, :) + delta;
        end
        R      = Rnew;
        spread = max(vecnorm(R - mean(R, 1), 2, 2));
        if spread < agreeTol, consIters = t; break; end
    end
    agreed = mean(R, 1);                          % common agreement point

    % ---- deviation from the centralized CRM-OPT optimum (pre-offset) ----
    instFeasible = inst;
    instFeasible.denied.vertices = [];            % compare like-for-like (no offset)
    rCentral = solveRendezvousCentral(params, instFeasible);
    gap      = norm(agreed - rCentral);

    % ---- 3) H-DRMi-D: offset the agreed point out of the denied region ----
    rStar     = agreed;
    deniedHit = false;  offsets = 0;
    if params.sim.denied && ~isempty(inst.denied.vertices)
        [rStar, applied] = escapeDenied(params, inst, rStar);
        deniedHit = applied;  offsets = double(applied);
    end

    info.proposals = proposals;
    info.consIters = consIters;
    info.gap       = gap;
    info.deniedHit = deniedHit;
    info.offsets   = offsets;
    info.agreed    = agreed;
end

% =====================================================================
function s = subInstance(inst, idx)
%SUBINSTANCE  Instance restricted to the agents in IDX (shared RoI / denied).
    s = inst;
    s.pos = inst.pos(idx, :);
    s.vel = inst.vel(idx, :);
    f = fieldnames(inst.cap);
    for k = 1:numel(f)
        s.cap.(f{k}) = inst.cap.(f{k})(idx);
    end
end
