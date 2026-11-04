function inst = generateInstance(params)
%GENERATEINSTANCE  Create a random heterogeneous rendezvous instance.
%
%   inst = GENERATEINSTANCE(params) draws a reproducible random instance
%   from the configuration PARAMS (see DEFAULTPARAMS) and returns a struct:
%
%     inst.agentType   N x 1   type index of each agent (1..K)
%     inst.pos         N x dim initial positions, m
%     inst.vel         N x dim initial velocities, m/s
%     inst.speed       N x 1   initial speed magnitudes, m/s
%     inst.cap         struct  per-agent capability vectors (N x 1 each):
%                              powerA, powerB, vMin, vMax, aMin, aMax,
%                              commRange, battery
%     inst.power       N x 1   initial power consumption w_i(v_i), W
%     inst.weight      N x 1   rendezvous weights c_i = power_i / speed_i
%     inst.centroid    1 x dim weighted centroid (unconstrained CRM optimum)
%     inst.battery     struct  voltage (N x 1), qDischarged (N x 1), soc (N x 1)
%     inst.roi         struct  region of interest: vertices, hull
%     inst.denied      struct  denied region: vertices, hull, theta, sigma
%
%   Initial positions satisfy, by construction:
%     (1) no agent's (x, y) lies in the denied no-fly column (the 2-D polygon
%         extruded over the arena height);
%     (2) every agent has at least one mutual communication link (degree >= 1),
%         so no agent is isolated at t = 0 -- required for D-MAS, where an
%         isolated agent cannot form a local rendezvous estimate or move;
%     (3) if params.sim.forceDeniedCentroid is true, the denied column is
%         placed on the weighted centroid so the unconstrained optimum is
%         infeasible and H-CRM-D always activates.
%
%   Initial speeds are drawn with magnitude in [0.1, vMax_i] so the weight
%   c_i = power_i / speed_i is always well defined.
%
%   See also DEFAULTPARAMS, INDENIEDREGION, BUILDCOMMGRAPH.

    arguments
        params (1,1) struct
    end

    rng(params.sim.seed);                 % reproducible instance

    N   = params.sim.N;
    dim = params.sim.dim;
    K   = numel(params.types);
    W   = params.arena.width;
    H   = params.arena.height;

    % ---------------------------------------------------- agent types
    if params.sim.homogeneous
        agentType = ones(N, 1);
    else
        agentType = ceil((1:N)' * K / N); % split as evenly as possible across K types
    end
    inst.agentType = agentType;
    inst.id        = (1:N).';          % persistent agent identifiers

    % --------------------------------------- per-agent capability vectors
    T = params.types;
    inst.cap.powerA    = [T(agentType).powerA].';
    inst.cap.powerB    = [T(agentType).powerB].';
    inst.cap.vMin      = [T(agentType).vMin].';
    inst.cap.vMax      = [T(agentType).vMax].';
    inst.cap.aMin      = [T(agentType).aMin].';
    inst.cap.aMax      = [T(agentType).aMax].';
    inst.cap.commRange = [T(agentType).commRange].';
    inst.cap.battery   = [T(agentType).battery].';

    % --------------------------------- velocities, power, weights
    % (independent of position; the rendezvous weight uses the speed only)
    direction = randn(N, dim);
    direction = direction ./ vecnorm(direction, 2, 2);  % unit heading
    speedMag  = 0.1 + (inst.cap.vMax - 0.1) .* rand(N, 1);
    inst.vel   = direction .* speedMag;
    inst.speed = vecnorm(inst.vel, 2, 2);
    inst.power = inst.cap.powerA .* inst.speed.^2 + inst.cap.powerB;
    weight     = inst.power ./ inst.speed;               % c_i
    inst.weight = weight;

    % ------------------------------------- region of interest (RoI)
    roiV = params.region.roi;
    inst.roi.vertices = roiV;
    inst.roi.hull     = convhull(roiV(:, 1), roiV(:, 2));

    % ------------------------------------------------ denied region
    if params.sim.denied
        dV = params.region.denied;
        inst.denied.vertices = dV;
        inst.denied.hull     = convhull(dV(:, 1), dV(:, 2));
        [theta, sigma]       = deniedRegionParams(dV);
        inst.denied.theta    = theta;
        inst.denied.sigma    = sigma;
    else
        inst.denied = struct('vertices', [], 'hull', [], 'theta', [], 'sigma', []);
    end

    % --------------------------------------------- initial positions
    % One rejection loop enforces every constraint simultaneously: not in the
    % denied column, degree >= 1, and (optional) denied pinned to the centroid.
    range  = inst.cap.commRange;
    denied = params.sim.denied;
    forced = denied && params.sim.forceDeniedCentroid;
    if forced
        baseV     = inst.denied.vertices;
        thetaBase = inst.denied.theta;
    end

    pos = W * rand(N, dim);
    if dim == 3
        pos(:, 3) = H * rand(N, 1);                 % altitude
    end

    maxIter = 2000;
    for it = 1:maxIter
        if forced
            c = sum(weight .* pos, 1) ./ sum(weight);      % weighted centroid
            offset = c(1:2) - thetaBase;                   % pin denied onto centroid
            inst.denied.vertices = baseV + offset;
            inst.denied.theta    = thetaBase + offset;
        end

        bad = mutualDegree(pos, range) < 1;                % isolated agents
        if denied
            bad = bad | inDeniedRegion(pos, inst.denied);  % agents in denied column
        end
        if ~any(bad), break; end

        m = nnz(bad);
        pos(bad, 1) = W * rand(m, 1);
        pos(bad, 2) = W * rand(m, 1);
        if it == maxIter
            warning('generateInstance:constraints', ...
                'Could not satisfy all placement constraints in %d iterations; %d agent(s) remain.', ...
                maxIter, m);
        end
    end

    inst.pos      = pos;
    inst.centroid = sum(weight .* pos, 1) ./ sum(weight);

    % ------------------------------------------------ battery state
    inst.battery.voltage     = 4.2 * ones(N, 1);   % full-charge voltage, V
    inst.battery.qDischarged = zeros(N, 1);         % charge drawn, mAh
    inst.battery.soc         = ones(N, 1);          % state of charge in [0,1]
end

% =====================================================================
function d = mutualDegree(P, range)
%MUTUALDEGREE  Degree of each agent in the undirected (mutual) comm graph.
%   A mutual link exists between i and j when each lies within the other's
%   range, i.e. dist(i,j) <= min(range_i, range_j). Matches BUILDCOMMGRAPH.
    N  = size(P, 1);
    sq = sum(P.^2, 2);
    D  = sqrt(max(sq + sq.' - 2 * (P * P.'), 0));
    D(1:N+1:end) = Inf;                 % exclude self
    Adir = D <= range;                  % i hears j (row threshold range_i)
    d    = sum(Adir & Adir.', 2);
end
