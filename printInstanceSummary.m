function printInstanceSummary(params, inst)
%PRINTINSTANCESUMMARY  Print a concise summary of a generated instance.
%
%   PRINTINSTANCESUMMARY(params, inst) reports the agent-type breakdown,
%   arena size, the weighted centroid, and whether that centroid lies inside
%   the region of interest and the denied region (which determines whether
%   the H-CRM-D offset heuristic will be exercised).
%
%   See also GENERATEINSTANCE, DEFAULTPARAMS.

arguments
    params (1,1) struct
    inst   (1,1) struct
end

K = numel(params.types);
fprintf('--- Rendezvous instance ---------------------------------\n');
fprintf('Agents      : N = %d in %d-D\n', params.sim.N, params.sim.dim);
for k = 1:K
    fprintf('  Type %d    : %2d  %s\n', k, nnz(inst.agentType == k), params.types(k).name);
end
fprintf('Arena       : %g x %g m footprint, %g m height\n', ...
    params.arena.width, params.arena.width, params.arena.height);

c = inst.centroid;
if params.sim.dim == 3
    fprintf('Weighted centroid : (%.2f, %.2f, %.2f) m\n', c(1), c(2), c(3));
else
    fprintf('Weighted centroid : (%.2f, %.2f) m\n', c(1), c(2));
end

inRoI = inpolygon(c(1), c(2), inst.roi.vertices(:,1), inst.roi.vertices(:,2));
fprintf('Centroid in RoI   : %s\n', ternary(inRoI, 'yes', 'no'));

if params.sim.denied && ~isempty(inst.denied.vertices)
    inD = inpolygon(c(1), c(2), inst.denied.vertices(:,1), inst.denied.vertices(:,2));
    fprintf('Centroid in denied: %s  (H-CRM-D %s)\n', ...
        ternary(inD, 'yes', 'no'), ternary(inD, 'will activate', 'not triggered'));
    fprintf('Agents in denied  : %d  (should be 0)\n', nnz(inDeniedRegion(inst.pos, inst.denied)));
end
fprintf('---------------------------------------------------------\n');
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
