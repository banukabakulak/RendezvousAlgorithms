function commGraph = buildCommGraph(params, inst)
%BUILDCOMMGRAPH  Distance-based communication graph with heterogeneous ranges.
%
%   commGraph = BUILDCOMMGRAPH(params, inst) builds the inter-agent
%   communication graph for the current positions, using each agent's own
%   effective range omega_i (inst.cap.commRange). Returns a struct:
%
%     commGraph.N      number of agents
%     commGraph.dist   N x N  Euclidean distances (all spatial dimensions)
%     commGraph.range  N x 1  per-agent communication range omega_i, m
%     commGraph.Adir   N x N  logical directed adjacency: Adir(i,j) is true
%                              when agent i can hear agent j, i.e.
%                              dist(i,j) <= omega_i. Row i lists the
%                              neighbour set N_i used by the message/overhead
%                              metrics.
%     commGraph.A      N x N  logical undirected adjacency: a mutual link
%                              exists when both agents hear each other
%                              (dist(i,j) <= min(omega_i, omega_j)). Used for
%                              the connectivity analysis.
%     commGraph.deg    N x 1  undirected degree (sum of A over rows)
%     commGraph.indeg  N x 1  directed in-degree |N_i| (sum of Adir over rows)
%
%   The undirected (mutual) graph is the conservative choice for heterogeneous
%   ranges and is what the consensus/connectivity theory refers to; the
%   directed adjacency is retained for the communication-cost metrics.
%
%   See also GENERATEINSTANCE, GRAPHCONNECTIVITY.

arguments
    params (1,1) struct %#ok<INUSA>  (kept for interface symmetry)
    inst   (1,1) struct
end

P     = inst.pos;
N     = size(P, 1);
range = inst.cap.commRange;                 % N x 1

% pairwise Euclidean distances (toolbox-free)
sq   = sum(P.^2, 2);
dist = sqrt(max(sq + sq.' - 2 * (P * P.'), 0));
dist(1:N+1:end) = 0;                         % clean diagonal

notSelf = ~logical(eye(N));
Adir = (dist <= range) & notSelf;           % i hears j  (row threshold omega_i)
A    = Adir & Adir.';                        % mutual (undirected) link

commGraph.N     = N;
commGraph.dist  = dist;
commGraph.range = range;
commGraph.Adir  = Adir;
commGraph.A     = A;
commGraph.deg   = sum(A, 2);
commGraph.indeg = sum(Adir, 2);
end
