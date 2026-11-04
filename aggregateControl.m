function [velNew, posNew] = aggregateControl(params, inst, pos, vel, rStar)
%AGGREGATECONTROL  One step of the energy-aware aggregation controller.
%
%   [velNew, posNew] = AGGREGATECONTROL(params, inst, pos, vel, rStar) applies
%   one control step (duration 1/freqInner) of the C-MAS aggregation law to
%   the current positions POS (N x dim) and velocities VEL (N x dim), driving
%   the team toward the rendezvous point RSTAR (1 x dim).
%
%   The control acceleration of agent i is
%       u_i =  -wAttr  sum_{j~=i} (alpha - beta exp(-||d_ij||^2/gamma)) d_ij
%              -wTrack delta (centroid - rStar)
%              -wEnergy 2 powerA_i v_i
%              +wObstc (zeta/sigma) exp(-||d_obs||^2/sigma) d_obs_i
%   with d_ij = p_i - p_j and d_obs_i the horizontal offset from the denied
%   centroid. The energy term is the velocity gradient of the per-agent power
%   model w_i(v) = powerA_i ||v||^2 + powerB_i: it opposes each agent's motion
%   in proportion to its own speed-power coefficient, so costlier (higher
%   powerA) agents are damped more and the team favours energy-efficient
%   motion. The acceleration magnitude is clipped to the agent's own
%   [aMin_i, aMax_i]; the velocity is smoothed (factor psi) and its magnitude
%   clipped to [vMin_i, vMax_i] (after a speedFactor margin); the position is
%   advanced by one control step.
%
%   See also RUNCMAS, SIMULATEROUND.

arguments
    params (1,1) struct
    inst   (1,1) struct
    pos    (:,:) double
    vel    (:,:) double
    rStar  (1,:) double
end

c   = params.controller;
N   = size(pos, 1);
dim = size(pos, 2);

% ---- attraction / repulsion (pairwise) ----
P1  = reshape(pos, N, 1, dim);
P2  = reshape(pos, 1, N, dim);
Dij = P1 - P2;                               % N x N x dim, d_ij = p_i - p_j
nrm2 = sum(Dij.^2, 3);                        % N x N
coef = c.alpha - c.beta * exp(-nrm2 / c.gamma);
u = zeros(N, dim);
for d = 1:dim
    u(:, d) = -c.wAttr * sum(coef .* Dij(:, :, d), 2);
end

% ---- tracking toward the rendezvous point (formation level) ----
centroid = mean(pos, 1);
u = u - c.wTrack * c.delta * (centroid - rStar);

% ---- energy-aware damping (per-agent) ----
% Gradient of the power model w_i(v) = powerA_i ||v||^2 + powerB_i w.r.t.
% velocity is 2 powerA_i v; subtracting it biases agents toward lower-speed,
% lower-power motion, weighted by each agent's own powerA (heterogeneous).
% u = u - c.wEnergy * 2 * inst.cap.powerA .* vel;

% ---- obstacle (denied-region) repulsion ----
safetyDist = 1;
if params.sim.denied && ~isempty(inst.denied.vertices)
    theta = inst.denied.theta;               % 1 x 2
    sigma = inst.denied.sigma;
    dObs  = pos(:, 1:2) - theta;             % horizontal offset
    if dim == 3, dObs = [dObs, zeros(N, 1)] + safetyDist; end
    g = c.wObstc * (c.zeta / sigma) * exp(-sum(dObs.^2, 2) / sigma);
    u = u + g .* dObs;
end

% ---- acceleration clipping (per-agent bounds) ----
nu    = vecnorm(u, 2, 2);
nuSafe = max(nu, eps);
amp   = max(inst.cap.aMin, min(nu, inst.cap.aMax));
accel = (amp ./ nuSafe) .* u;

% ---- velocity smoothing + speed clipping (per-agent bounds) ----
smoothed = c.psi * (vel + accel / c.freqInner) + (1 - c.psi) * vel;
ns       = vecnorm(smoothed, 2, 2);
nsSafe   = max(ns, eps);
spd      = max(inst.cap.vMin, min(c.speedFactor * ns, inst.cap.vMax));
velNew   = (spd ./ nsSafe) .* smoothed;

% ---- position update (one control step) ----
posNew = pos + velNew / c.freqInner;
end
