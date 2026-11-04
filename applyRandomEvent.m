function inst = applyRandomEvent(params, inst, eventType)
%APPLYRANDOMEVENT  Inter-round event: add or remove an agent, with perturbation.
%
%   inst = APPLYRANDOMEVENT(params, inst, eventType) applies one inter-round
%   event and returns the updated instance:
%     'add'    perturb every agent (up to 2 m), then add one new agent of a
%              random type near a corner with a random initial velocity.
%     'remove' remove one random agent, then perturb the remaining agents.
%     'none'   no change.
%
%   Positions are clipped to the arena. All per-agent capability, weight and
%   battery fields are kept consistent.
%
%   See also RUNCMAS, GENERATEINSTANCE.

    arguments
        params    (1,1) struct
        inst      (1,1) struct
        eventType (1,:) char {mustBeMember(eventType, {'add','remove','none'})} = 'add'
    end

    switch eventType
        case 'add'
            inst = perturbPositions(params, inst);
            inst = addAgent(params, inst);
        case 'remove'
            inst = removeAgent(params, inst);
            inst = perturbPositions(params, inst);
        case 'none'
            % no-op
    end
end

% =====================================================================
function inst = perturbPositions(params, inst)
%PERTURBPOSITIONS  Perturb all positions by up to 2 m (single random sign).
    N   = numel(inst.agentType);
    dim = params.sim.dim;
    s   = (-1)^randi([0 1]);
    inst.pos = inst.pos + s * 2 * rand(N, dim);
    inst.pos(:, 1:2) = min(max(inst.pos(:, 1:2), 0), params.arena.width);
    if dim == 3
        inst.pos(:, 3) = min(max(inst.pos(:, 3), 0), params.arena.height);
    end
end

% =====================================================================
function inst = addAgent(params, inst)
%ADDAGENT  Append one new agent of a random type near a corner.
    dim = params.sim.dim;
    K   = numel(params.types);
    k   = randi([1 K]);
    T   = params.types(k);

    newPos = 0.25 * params.arena.width * rand(1, dim);
    if dim == 3, newPos(3) = 0.25 * params.arena.height * rand; end

    dir    = randn(1, dim); dir = dir / norm(dir);
    sp     = 0.1 + (T.vMax - 0.1) * rand;
    newVel = dir * sp;

    inst.agentType(end+1, 1) = k;
    inst.id(end+1, 1)        = max(inst.id) + 1;   % fresh persistent id
    inst.pos(end+1, :)       = newPos;
    inst.vel(end+1, :)       = newVel;

    inst.cap.powerA(end+1,1)    = T.powerA;
    inst.cap.powerB(end+1,1)    = T.powerB;
    inst.cap.vMin(end+1,1)      = T.vMin;
    inst.cap.vMax(end+1,1)      = T.vMax;
    inst.cap.aMin(end+1,1)      = T.aMin;
    inst.cap.aMax(end+1,1)      = T.aMax;
    inst.cap.commRange(end+1,1) = T.commRange;
    inst.cap.battery(end+1,1)   = T.battery;

    inst.speed(end+1,1)  = sp;
    inst.power(end+1,1)  = T.powerA * sp^2 + T.powerB;
    inst.weight(end+1,1) = inst.power(end) / sp;

    inst.battery.voltage(end+1,1)     = 4.2;
    inst.battery.qDischarged(end+1,1) = 0;
    inst.battery.soc(end+1,1)         = 1;
end

% =====================================================================
function inst = removeAgent(params, inst) %#ok<INUSD>
%REMOVEAGENT  Remove one random agent (keeps at least 2 agents).
    N = numel(inst.agentType);
    if N <= 2, return; end
    keep = true(N, 1);
    keep(randi([1 N])) = false;

    inst.agentType = inst.agentType(keep);
    inst.id        = inst.id(keep);
    inst.pos       = inst.pos(keep, :);
    inst.vel       = inst.vel(keep, :);
    inst.speed     = inst.speed(keep);
    inst.power     = inst.power(keep);
    inst.weight    = inst.weight(keep);

    f = fieldnames(inst.cap);
    for q = 1:numel(f), inst.cap.(f{q}) = inst.cap.(f{q})(keep); end

    f = fieldnames(inst.battery);
    for q = 1:numel(f), inst.battery.(f{q}) = inst.battery.(f{q})(keep); end
end
