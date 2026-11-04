function jc = jointConnectivity(params, traj)
%JOINTCONNECTIVITY  Empirical joint-connectivity check over sliding windows.
%
%   jc = JOINTCONNECTIVITY(params, traj) tests the joint-connectivity
%   assumption on the per-step mutual communication graphs logged in TRAJ.adj
%   (nStep x N x N logical, from SIMULATEROUND). For every window of
%   params.controller.jointWindow consecutive control steps it forms the
%   \emph{union} graph (an edge is present if it appears at any step in the
%   window) and evaluates:
%     * connectivity of the union (single connected component), and
%     * its algebraic connectivity lambda_2 (second-smallest Laplacian
%       eigenvalue, the Fiedler value; > 0 iff the union is connected).
%
%   Returns a struct:
%     jc.window          window length in control steps (jointWindow)
%     jc.nWindows        number of sliding windows evaluated
%     jc.windowLambda2   nWindows x 1  lambda_2 of each window's union graph
%     jc.connectedWindow nWindows x 1  logical, union connected per window
%     jc.allConnected    true if every window's union is connected
%     jc.minLambda2      smallest lambda_2 over all windows (the binding case)
%
%   The assumption requires the union over EVERY length-jointWindow interval to
%   be connected, so all overlapping (stride-1) windows are checked. If TRAJ.adj
%   was not logged (params.controller.logGraph false), the fields are empty/NaN.
%
%   See also SIMULATEROUND, BUILDCOMMGRAPH, GRAPHCONNECTIVITY, RUNDMAS.

arguments
    params (1,1) struct
    traj   (1,1) struct
end

if ~isfield(traj, 'adj') || isempty(traj.adj)
    jc = struct('window', NaN, 'nWindows', 0, 'windowLambda2', [], ...
        'connectedWindow', [], 'allConnected', NaN, 'minLambda2', NaN);
    return;
end

A     = traj.adj;                       % nStep x N x N logical
nStep = size(A, 1);
N     = size(A, 2);
win   = min(params.controller.jointWindow, nStep);

ends  = win:nStep;                      % window [e-win+1, e], stride 1
nWin  = numel(ends);
lam   = zeros(nWin, 1);
conn  = false(nWin, 1);

for w = 1:nWin
    e = ends(w);  s = e - win + 1;
    U = squeeze(any(A(s:e, :, :), 1));  % N x N union adjacency
    U = (U | U.');                       % symmetrize (mutual graph already is)
    U(1:N+1:end) = false;                % clear diagonal

    G       = graph(U, 'omitselfloops');
    conn(w) = max(conncomp(G)) == 1;

    L  = diag(sum(U, 2)) - double(U);    % Laplacian
    ev = sort(eig(L));
    if numel(ev) >= 2, lam(w) = ev(2); else, lam(w) = 0; end
end

jc.window          = win;
jc.nWindows        = nWin;
jc.windowLambda2   = lam;
jc.connectedWindow = conn;
jc.allConnected    = all(conn);
jc.minLambda2      = min(lam);
end
