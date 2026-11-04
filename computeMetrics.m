function perf = computeMetrics(params, results)
%COMPUTEMETRICS  Aggregate C-MAS performance metrics from a run.
%
%   perf = COMPUTEMETRICS(params, results) summarizes the RESULTS returned by
%   RUNCMAS into per-round and total metrics:
%
%     perf.round(r).N          agent count in round r
%     perf.round(r).step       consensus step (full round if never reached)
%     perf.round(r).time       convergence time, s
%     perf.round(r).duty       duty cycle, %
%     perf.round(r).overhead   communication overhead, Kbits  (C-MAS)
%     perf.round(r).messages   number of messages exchanged    (C-MAS)
%     perf.round(r).energy     charge drawn in the round, mAh
%
%     perf.totalEnergy     total charge drawn, mAh
%     perf.totalOverhead   total communication overhead, Kbits
%     perf.totalMessages   total messages exchanged
%     perf.avgConsensusTime, perf.avgDutyCycle
%
%   Overhead and message models are the centralized ones:
%     O_CMAS = N * kappaC * (nu2 * Tc),   M_CMAS = N * (nu2 * Tc),
%   with (nu2 * Tc) measured in control steps.
%
%   See also RUNCMAS, SIMULATEROUND, ANALYZECONSENSUS.

arguments
    params  (1,1) struct
    results (1,1) struct
end

nRounds = numel(results.traj);
nStep   = round(params.controller.freqInner / params.controller.freqOuter);
kappaC  = params.comm.kappaC;

perf.totalEnergy   = 0;
perf.totalOverhead = 0;
perf.totalMessages = 0;
times = nan(1, nRounds);
duties = nan(1, nRounds);

for r = 1:nRounds
    traj = results.traj{r};
    Nr   = size(traj.P, 2);

    step = results.consensus(r).step;
    if isnan(step), step = nStep; end          % no consensus: charge the full round

    overhead = Nr * kappaC * step;             % Kbits
    messages = Nr * step;                       % C-MAS message count
    energy   = traj.roundEnergy;                % mAh

    perf.round(r).N        = Nr;
    perf.round(r).step     = step;
    perf.round(r).time     = results.consensus(r).time;
    perf.round(r).duty     = results.consensus(r).duty;
    perf.round(r).overhead = overhead;
    perf.round(r).messages = messages;
    perf.round(r).energy   = energy;

    perf.totalEnergy   = perf.totalEnergy   + energy;
    perf.totalOverhead = perf.totalOverhead + overhead;
    perf.totalMessages = perf.totalMessages + messages;
    times(r)  = results.consensus(r).time;
    duties(r) = results.consensus(r).duty;
end

perf.avgConsensusTime = mean(times, 'omitnan');
perf.avgDutyCycle     = mean(duties, 'omitnan');

% joint connectivity (over all sliding windows of every round)
if isfield(results, 'conn') && ~isempty(results.conn)
    perf.jointConnected = all([results.conn.allConnected]);
    perf.minLambda2     = min([results.conn.minLambda2]);
else
    perf.jointConnected = NaN;  perf.minLambda2 = NaN;
end

fprintf('--- C-MAS performance -----------------------------------\n');
fprintf('Total energy     : %.2f mAh\n', perf.totalEnergy);
fprintf('Total overhead   : %.1f Kbits\n', perf.totalOverhead);
fprintf('Total messages   : %d\n', perf.totalMessages);
fprintf('Avg convergence  : %.1f s\n', perf.avgConsensusTime);
fprintf('Avg duty cycle   : %.1f %%\n', perf.avgDutyCycle);
fprintf('Joint connected  : %d  (min lambda2 = %.3f)\n', perf.jointConnected, perf.minLambda2);
fprintf('---------------------------------------------------------\n');
end
