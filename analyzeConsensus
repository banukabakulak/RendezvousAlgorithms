function cons = analyzeConsensus(params, traj, rStar)
%ANALYZECONSENSUS  Convergence time and duty cycle for one round.
%
%   cons = ANALYZECONSENSUS(params, traj, rStar) scans the position
%   trajectory TRAJ.P (nStep x N x dim) and detects rendezvous: the first time
%   all agents stay within the tolerance of RSTAR for consensusConsec
%   consecutive control steps. Returns:
%
%     cons.step   control step at which consensus is first sustained
%     cons.time   convergence time, s  (= step / freqInner)
%     cons.duty   duty cycle, %         (= step * 100 / (freqInner/freqOuter))
%
%   If consensus is never sustained, cons.step is NaN.
%
%   See also RUNCMAS, SIMULATEROUND.

arguments
    params (1,1) struct
    traj   (1,1) struct
    rStar  (1,:) double
end

epsC   = params.controller.consensusEps;      % squared-distance tolerance
need   = params.controller.consensusConsec;   % consecutive steps required
P      = traj.P;
nStep  = size(P, 1);
N      = size(P, 2);

startStep = NaN;
count     = 0;
for t = 1:nStep
    allIn = true;
    for i = 1:N
        d2 = sum((squeeze(P(t, i, :)).' - rStar).^2);
        if d2 > epsC
            allIn = false; break;
        end
    end
    if allIn
        if count == 0
            startStep = t; count = 1;
        else
            count = count + 1;
        end
    else
        count = 0; startStep = NaN;
    end
    if count >= need
        break;
    end
end

cons.step = startStep;
if isnan(startStep)
    cons.time = NaN;  cons.duty = NaN;
else
    cons.time = startStep / params.controller.freqInner;
    cons.duty = startStep * 100 / (params.controller.freqInner / params.controller.freqOuter);
end
end
