function [rStar, info] = solveRendezvousCentral(params, inst)
%SOLVERENDEZVOUSCENTRAL  Centralized rendezvous point (CRM-OPT + H-CRM-D).
%
%   [rStar, info] = SOLVERENDEZVOUSCENTRAL(params, inst) computes the
%   energy-weighted rendezvous point for the current state:
%
%     1. CRM-OPT solves  min_r  sum_i c_i || r - p_i ||^2
%        over r = sum_j lambda_j v_j  with  sum lambda_j = 1, lambda_j >= 0,
%        where v_j are the RoI extreme points and c_i = w_i(v_i)/||v_i|| is the
%        per-agent energy weight (heterogeneous). The planar (x,y) component is
%        obtained from the KKT system by an active-set method; the altitude is
%        the weighted centroid height.
%     2. If a denied region is present and the optimum lies inside it, H-CRM-D
%        offsets the point out of the denied column (see ESCAPEDENIED).
%
%   Returns rStar (1 x dim) and info with fields:
%     lambda   M x 1  convex weights on the RoI vertices
%     iters    active-set iterations
%     deniedHit, offsets  H-CRM-D activation and number of offset steps
%
%   See also ESCAPEDENIED, GENERATEINSTANCE.

    arguments
        params (1,1) struct
        inst   (1,1) struct
    end

    dim   = params.sim.dim;
    pos   = inst.pos;
    speed = vecnorm(inst.vel, 2, 2);
    c     = (inst.cap.powerA .* speed.^2 + inst.cap.powerB) ./ speed;  % weights c_i

    V = inst.roi.vertices(:, 1:2);        % RoI extreme points (planar)
    M = size(V, 1);
    N = size(pos, 1);

    scaleK = 1e3;                          % KKT conditioning scale
    epsA   = 0.01;                         % active-set tolerance

    % ---- assemble the KKT system  A x = b,  x = [lambda; pi; mu] ----
    A = zeros(2*M+1, 2*M+1);
    b = zeros(2*M+1, 1);

    A(M+1, 1:M) = 1;                       % sum(lambda) = 1  (pi coupling)
    A(1:M, M+1) = 1;
    for i = (M+2):(2*M+1)                  % nonnegativity multipliers (inactive init)
        A(i-M-1, i) = -1;
        A(i, i)     =  1;
    end

    % Hessian block  H(j,jj) = ( sum_i 2 c_i ) <v_j, v_jj> / scaleK,
    % and linear term  b(j) = sum_i 2 c_i <v_j, p_i> / scaleK.
    W2 = 2 * sum(c);                       % sum_i 2 c_i  (vertices are shared)
    A(1:M, 1:M) = (W2 / scaleK) * (V * V.');
    b(M+1) = 1;
    b(1:M) = (2 / scaleK) * (V * (pos(:,1:2).' * c));   % sum_i 2 c_i <v_j, p_i>

    % ---- active-set iterations ----
    loop  = true;
    iters = 0;
    while loop
        x    = A \ b;
        loop = false;
        iters = iters + 1;
        for i = 1:M                        % activate a violated lambda_i >= 0
            if x(i) + epsA < 0
                A(M+1+i, i)     = 1;
                A(M+1+i, M+1+i) = 0;
                loop = true; break;
            end
        end
        if ~loop
            for i = 1:M                    % release an over-constrained multiplier
                if x(M+1+i) + epsA < 0
                    A(M+1+i, i)     = 0;
                    A(M+1+i, M+1+i) = 1;
                    loop = true; break;
                end
            end
        end
    end

    lambda = x(1:M);
    rxy    = (lambda.' * V);               % 1 x 2

    % altitude = weighted-centroid height
    if dim == 3
        zc    = sum(c .* pos(:,3)) / sum(c);
        rStar = [rxy, zc];
    else
        rStar = rxy;
    end

    info.lambda    = lambda;
    info.iters     = iters;
    info.deniedHit = false;
    info.offsets   = 0;

    % ---- H-CRM-D: offset out of a denied region if needed ----
    if params.sim.denied && ~isempty(inst.denied.vertices)
        [rStar, applied] = escapeDenied(params, inst, rStar);
        info.deniedHit   = applied;
        info.offsets     = double(applied);
    end
end
