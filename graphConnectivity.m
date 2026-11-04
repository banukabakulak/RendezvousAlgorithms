function metrics = graphConnectivity(commGraph)
%GRAPHCONNECTIVITY  Connectivity metrics of the communication graph.
%
%   metrics = GRAPHCONNECTIVITY(commGraph) analyses the undirected (mutual)
%   graph in COMMGRAPH (see BUILDCOMMGRAPH) and returns a struct:
%
%     metrics.connected      true if the graph is a single component
%     metrics.numComponents  number of connected components
%     metrics.lambda2        algebraic connectivity (second-smallest Laplacian
%                            eigenvalue, the Fiedler value); > 0 iff connected
%     metrics.diameter       graph diameter in hops (Inf if disconnected)
%     metrics.meanDegree     mean undirected degree
%     metrics.minDegree      minimum undirected degree
%     metrics.maxDegree      maximum undirected degree
%
%   Uses base-MATLAB graph functions only (graph, conncomp, distances).
%
%   See also BUILDCOMMGRAPH, PRINTGRAPHMETRICS.

arguments
    commGraph (1,1) struct
end

A = double(commGraph.A);
G = graph(A ~= 0, 'omitselfloops');

comp = conncomp(G);
metrics.numComponents = max(comp);
metrics.connected     = metrics.numComponents == 1;

% algebraic connectivity (Fiedler value) of the Laplacian L = D - A
L  = diag(sum(A, 2)) - A;
ev = sort(eig(L));
metrics.lambda2 = ev(2);

% diameter in hops
hop = distances(G);                         % unweighted shortest paths
if metrics.connected
    metrics.diameter = max(hop(:));
else
    metrics.diameter = Inf;
end

deg = commGraph.deg;
metrics.meanDegree = mean(deg);
metrics.minDegree  = min(deg);
metrics.maxDegree  = max(deg);
end
