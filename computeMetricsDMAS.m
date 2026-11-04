function perf = computeMetricsDMAS(params, results)
%COMPUTEMETRICSDMAS  Aggregate D-MAS performance metrics from a run.
%
%   perf = COMPUTEMETRICSDMAS(params, results) summarizes the RESULTS returned
%   by RUNDMAS into per-round and total metrics. It mirrors COMPUTEMETRICS but
%   uses the DISTRIBUTED communication-cost model: every agent exchanges with
%   each of its |N_i| neighbours at every control step over the convergence
%   time T_c (= step control steps), so
%
%     O_D-MAS = sum_{t=1}^{nu2 Tc} sum_i |N_i(t)| kappaD   (Kbits)
%     M_D-MAS = sum_{t=1}^{nu2 Tc} sum_i |N_i(t)|           (messages)
%
%   The neighbour set is taken from the round's communication graph and held
%   fixed over the round (sum_i |N_i| = results.comm(r).sumDeg), so per round
%     overhead = kappaD * step * sumDeg,   messages = step * sumDeg.
%   Energy, convergence time and duty cycle are computed exactly as for C-MAS.
%
%     perf.round(r)   .N .step .time .duty .overhead .messages .energy
%                     .consIters .gap .connected
%     perf.totalEnergy, perf.totalOverhead, perf.totalMessages
%     perf.avgConsensusTime, perf.avgDutyCycle, perf.avgConsIters, perf.avgGap
%
%   See also RUNDMAS, COMPUTEMETRICS, ANALYZECONSENSUS.

arguments
    params  (1,1) struct
    results (1,1) struct
end

nRounds = numel(results.traj);
nStep   = round(params.controller.freqInner / params.controller.freqOuter);
kappaD  = params.comm.kappaD;

perf.totalEnergy   = 0;
perf.totalOverhead = 0;
perf.totalMessages = 0;
times   = nan(1, nRounds);
duties  = nan(1, nRounds);
citers  = nan(1, nRounds);
gaps    = nan(1, nRounds);

for r = 1:nRounds
    traj = results.traj{r};
    Nr   = size(traj.P, 2);

    step = results.consensus(r).step;
    if isnan(step), step = nStep; end           % no consensus: charge the full round

    sumDeg   = results.comm(r).sumDeg;          % sum_i |N_i| (static over the round)
    overhead = kappaD * step * sumDeg;          % Kbits  (Eq. DMAS_Overhead)
    messages = step * sumDeg;                    % messages (Eq. DMAS_Message)
    energy   = traj.roundEnergy;                 % mAh

    perf.round(r).N         = Nr;
    perf.round(r).step      = step;
    perf.round(r).time      = results.consensus(r).time;
    perf.round(r).duty      = results.consensus(r).duty;
    perf.round(r).overhead  = overhead;
    perf.round(r).messages  = messages;
    perf.round(r).energy    = energy;
    perf.round(r).consIters = results.info(r).consIters;
    perf.round(r).gap       = results.info(r).gap;
    perf.round(r).connected = results.comm(r).connected;

    perf.totalEnergy   = perf.totalEnergy   + energy;
    perf.totalOverhead = perf.totalOverhead + overhead;
    perf.totalMessages = perf.totalMessages + messages;
    times(r)  = results.consensus(r).time;
    duties(r) = results.consensus(r).duty;
    citers(r) = results.info(r).consIters;
    gaps(r)   = results.info(r).gap;
end

perf.avgConsensusTime = mean(times,  'omitnan');
perf.avgDutyCycle     = mean(duties, 'omitnan');
perf.avgConsIters     = mean(citers, 'omitnan');
perf.avgGap           = mean(gaps,   'omitnan');

% joint connectivity (over all sliding windows of every round)
if isfield(results, 'conn') && ~isempty(results.conn)
    perf.jointConnected = all([results.conn.allConnected]);
    perf.minLambda2     = min([results.conn.minLambda2]);
else
    perf.jointConnected = NaN;  perf.minLambda2 = NaN;
end

fprintf('--- D-MAS performance -----------------------------------\n');
fprintf('Total energy      : %.2f mAh\n', perf.totalEnergy);
fprintf('Total overhead    : %.1f Kbits\n', perf.totalOverhead);
fprintf('Total messages    : %d\n', perf.totalMessages);
fprintf('Avg convergence   : %.1f s\n', perf.avgConsensusTime);
fprintf('Avg duty cycle    : %.1f %%\n', perf.avgDutyCycle);
fprintf('Avg consensus iter: %.1f\n', perf.avgConsIters);
fprintf('Avg C-MAS gap     : %.3f m\n', perf.avgGap);
fprintf('Joint connected   : %d  (min lambda2 = %.3f)\n', perf.jointConnected, perf.minLambda2);
fprintf('---------------------------------------------------------\n');
end
