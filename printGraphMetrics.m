function printGraphMetrics(metrics, N)
%PRINTGRAPHMETRICS  Print a concise communication-graph summary.
%
%   PRINTGRAPHMETRICS(metrics)      prints the connectivity metrics returned
%   by GRAPHCONNECTIVITY.
%   PRINTGRAPHMETRICS(metrics, N)   also shows the degree relative to N-1.
%
%   See also GRAPHCONNECTIVITY, BUILDCOMMGRAPH.

arguments
    metrics (1,1) struct
    N       (1,1) double = NaN
end

fprintf('--- Communication graph ---------------------------------\n');
if metrics.connected
    fprintf('Connected   : yes\n');
else
    fprintf('Connected   : no  (%d components)\n', metrics.numComponents);
end
fprintf('lambda2     : %.4f  (algebraic connectivity)\n', metrics.lambda2);
if isfinite(metrics.diameter)
    fprintf('Diameter    : %d hops\n', metrics.diameter);
else
    fprintf('Diameter    : Inf (disconnected)\n');
end
if isnan(N)
    fprintf('Degree      : mean %.1f (min %d, max %d)\n', ...
        metrics.meanDegree, metrics.minDegree, metrics.maxDegree);
else
    fprintf('Degree      : mean %.1f of %d (min %d, max %d)\n', ...
        metrics.meanDegree, N-1, metrics.minDegree, metrics.maxDegree);
end
fprintf('---------------------------------------------------------\n');
end
