function [rStar, applied] = escapeDenied(params, inst, rStar)
%ESCAPEDENIED  H-CRM-D offset of a denied rendezvous point (ApplyRendezvousOffset).
%
%   [rStar, applied] = ESCAPEDENIED(params, inst, rStar) offsets the CRM-OPT
%   rendezvous point RSTAR out of the denied region when it falls inside it,
%   using the energy-weighted offset direction of the H-CRM-D algorithm. APPLIED
%   is true when an offset was performed. A single offset step is applied, as in
%   the original algorithm.
%
%   Port of ApplyRendezvousOffset with the per-agent energy weight
%   c_i = w_i(v_i)/||v_i||. The step length rho is taken over all denied
%   vertices, and the offset-direction sums are computed independently per
%   coordinate.
%
%   See also SOLVERENDEZVOUSCENTRAL, INDENIEDREGION.

arguments
    params (1,1) struct
    inst   (1,1) struct
    rStar  (1,:) double
end

eta   = 1;      % offset regularization
scale = 1; %0.5;    % offset step scale

denied = inst.denied;
applied = false;
if isempty(denied.vertices) || ~inDeniedRegion(rStar, denied)
    return;                                   % feasible: nothing to do
end
applied = true;

dim   = params.sim.dim;
pos   = inst.pos;
N     = size(pos, 1);
speed = vecnorm(inst.vel, 2, 2);
c     = (inst.cap.powerA .* speed.^2 + inst.cap.powerB) ./ speed;   % weights c_i

V  = denied.vertices;                          % denied extreme points (x,y)
Md = size(V, 1);

% ---- step length rho = max distance from any denied vertex to rStar ----
% A vertex sits on the ground (z = 0), so the distance includes the
% rendezvous altitude in 3-D.
rho = 0;
for i = 1:Md
    if dim == 3
        vd = [V(i,:), 0] - rStar;
    else
        vd = V(i,:) - rStar;
    end
    rho = max(rho, norm(vd));
end

% ---- offset direction (independent per coordinate) ----
nCoord = dim - 1;                              % offset in the horizontal plane
dBar   = zeros(1, nCoord);
for l = 1:nCoord
    sum1 = 0;
    sum2 = 0;
    for i = 1:N
        sum1 = sum1 + 2 * rho   * c(i) * (pos(i,l) - rStar(l));
        sum2 = sum2 + 2 * rho^2 * c(i);
    end
    dBar(l) = (sum1 - eta) / sum2;
end

roiCentroid = mean(inst.roi.vertices(:, 1:nCoord), 1);
dBar = dBar + roiCentroid;
d    = dBar / norm(dBar);

rStar(1:nCoord) = rStar(1:nCoord) + scale * rho * d;
end
