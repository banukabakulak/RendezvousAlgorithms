function tf = inDeniedRegion(pos, denied)
%INDENIEDREGION  Test whether points lie inside the denied region.
%
%   tf = INDENIEDREGION(pos, denied) returns an N x 1 logical vector that is
%   true for every point whose horizontal (x, y) projection falls inside the
%   denied region. The denied region is modelled as a vertical column (the
%   2-D polygon extruded over the full arena height), so altitude is ignored.
%
%   Inputs
%     pos    N x dim  positions (only the first two columns are used)
%     denied struct    with fields .vertices (Mx2) and .hull (boundary order);
%                       an empty region yields all-false.
%
%   See also GENERATEINSTANCE.

arguments
    pos    (:, :) double
    denied (1, 1) struct
end

if isempty(denied.vertices)
    tf = false(size(pos, 1), 1);
    return;
end

V = denied.vertices;
k = denied.hull;                      % boundary vertices in order
tf = inpolygon(pos(:, 1), pos(:, 2), V(k, 1), V(k, 2));
end
