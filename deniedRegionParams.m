function [theta, sigma] = deniedRegionParams(vertices)
%DENIEDREGIONPARAMS  Centroid and bounding radius of a denied-region polygon.
%
%   [theta, sigma] = DENIEDREGIONPARAMS(vertices) returns
%       theta  1 x 2   centroid of the (distinct) extreme points
%       sigma  scalar  maximum distance from theta to any extreme point
%
%   VERTICES is an Mx2 list of polygon vertices; a repeated closing vertex
%   (first == last) is ignored.
%
%   See also GENERATEINSTANCE.

arguments
    vertices (:, 2) double
end

V     = unique(vertices, 'rows', 'stable');   % drop a closing duplicate
theta = mean(V, 1);
sigma = max(vecnorm(V - theta, 2, 2));
end
