function [traj, inst] = simulateRound(params, inst, rStar)
%SIMULATEROUND  Integrate the aggregation dynamics over one rendezvous round.
%
%   [traj, inst] = SIMULATEROUND(params, inst, rStar) advances the MAS toward
%   RSTAR for one outer round (freqInner/freqOuter control steps) using an
%   explicit fixed-step integrator at the control rate freqInner, while
%   accounting for per-agent energy consumption and battery discharge.
%
%   Returns:
%     traj.P       nStep x N x dim   position trajectory
%     traj.t       nStep x 1         time stamps, s
%     traj.speed   nStep x N         speed magnitudes, m/s
%     traj.energy  nStep x N         cumulative charge drawn per agent, mAh
%     traj.roundEnergy                total charge drawn this round, mAh
%     traj.adj     nStep x N x N      per-step mutual communication graph
%                                     (logical); used for the joint-connectivity
%                                     check (JOINTCONNECTIVITY). Logged only when
%                                     params.controller.logGraph is true (default);
%                                     set it false to save memory for large N.
%     inst                           updated positions, velocities, battery
%
%   Energy per agent per step is  dE = w_i(||v_i||) / (3.6 * freqInner * V_i)
%   mAh (= 5 w /(18 nu2 V)); the battery voltage follows the piecewise SOC
%   model. Capacities and power coefficients are per-agent (heterogeneous).
%
%   See also AGGREGATECONTROL, RUNCMAS.

    arguments
        params (1,1) struct
        inst   (1,1) struct
        rStar  (1,:) double
    end

    nStep = round(params.controller.freqInner / params.controller.freqOuter);
    nu2   = params.controller.freqInner;
    N     = size(inst.pos, 1);
    dim   = size(inst.pos, 2);

    traj.P      = zeros(nStep, N, dim);
    traj.t      = (1:nStep).' / nu2;
    traj.speed  = zeros(nStep, N);
    traj.energy = zeros(nStep, N);

    % per-step mutual communication graph (for the joint-connectivity check)
    logGraph = ~isfield(params.controller, 'logGraph') || params.controller.logGraph;
    if logGraph
        traj.adj  = false(nStep, N, N);
        rangeMin  = min(inst.cap.commRange, inst.cap.commRange.');  % N x N, min(omega_i,omega_j)
        notSelf   = ~logical(eye(N));
    else
        traj.adj = [];
    end

    pos = inst.pos;
    vel = inst.vel;
    aPow = inst.cap.powerA;  bPow = inst.cap.powerB;  cap = inst.cap.battery;
    volt = inst.battery.voltage;  qDis = inst.battery.qDischarged;

    roundEnergy = 0;
    for k = 1:nStep
        [vel, pos] = aggregateControl(params, inst, pos, vel, rStar);

        speed = vecnorm(vel, 2, 2);
        power = aPow .* speed.^2 + bPow;            % W
        dE    = power ./ (3.6 * nu2 * volt);        % mAh this step (= 5 w /(18 nu2 V))
        qDis  = qDis + dE;
        soc   = 1 - qDis ./ cap;
        volt  = batteryVoltage(soc);
        roundEnergy = roundEnergy + sum(dE);

        traj.P(k, :, :)  = pos;
        traj.speed(k, :) = speed;
        traj.energy(k, :) = qDis;

        if logGraph                                  % instantaneous mutual graph
            sq = sum(pos.^2, 2);
            D  = sqrt(max(sq + sq.' - 2 * (pos * pos.'), 0));
            traj.adj(k, :, :) = (D <= rangeMin) & notSelf;
        end
    end

    inst.pos = pos;
    inst.vel = vel;
    inst.speed = vecnorm(vel, 2, 2);
    inst.battery.voltage     = volt;
    inst.battery.qDischarged = qDis;
    inst.battery.soc         = 1 - qDis ./ cap;

    traj.roundEnergy = roundEnergy;
end

% =====================================================================
function v = batteryVoltage(soc)
%BATTERYVOLTAGE  Piecewise LiPo discharge model (per-agent, vectorized).
    v = zeros(size(soc));
    hi = soc >= 0.2;
    mid = soc >= 0 & ~hi;
    v(hi)  = 4.2 - 0.625 * (1 - soc(hi));
    v(mid) = 3.0 + 3.5 * soc(mid);
    v(~(hi | mid)) = 3.0;                 % depleted: clamp to minimum
    if any(~(hi | mid))
        warning('simulateRound:batteryDepleted', 'One or more agents have depleted batteries.');
    end
end
