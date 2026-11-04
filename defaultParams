function params = defaultParams()
%DEFAULTPARAMS  Default configuration for the energy-efficient rendezvous study.
%
%   params = DEFAULTPARAMS() returns a struct holding every simulation
%   parameter, grouped into self-describing sub-structs:
%
%       params.sim         general simulation settings (N, dim, seed, ...)
%       params.arena       workspace geometry (width, height)
%       params.types       1xK struct array of heterogeneous agent types
%       params.controller  aggregation-controller gains and potential weights
%       params.comm        communication-model parameters
%       params.region      region-of-interest and denied-region polygons
%
%   The two heterogeneous agent types are modelled on the Bitcraze CrazyFlie
%   platform family:
%       Type 1  Scout          CrazyFlie 2.1            (nano quadcopter)
%       Type 2  Carrier/Relay  CrazyFlie 2.1 Brushless  (payload + relay)
%
%   All parameters may be overridden after calling this function, e.g.
%       params = defaultParams();  params.sim.N = 50;
%
%   See also GENERATEINSTANCE.

% ------------------------------------------------------------------ general
params.sim.N        = 30;      % number of agents
params.sim.dim      = 3;       % spatial dimension (2 or 3)
params.sim.seed     = 1;       % random-number-generator seed (reproducibility)
params.sim.homogeneous = false;% true => all agents are Type 1
params.sim.denied   = true;    % true => a denied region is present
params.sim.forceDeniedCentroid = true; % true => guarantee the weighted
                               % centroid lies in the denied region, i.e. an
                               % instance that activates H-CRM-D (used for the
                               % scalability study; leave false for a natural
                               % fixed-denied instance, e.g. the Sec. VI-B demo)

% ------------------------------------------------------------------- arena
params.arena.width  = 30;      % horizontal extent of the flight arena, m
params.arena.height = 10;      % vertical extent of the flight arena, m

% ------------------------------------------------- heterogeneous agent types
% Power model:  w_i(v) = powerA * ||v||^2 + powerB   [W]
% Type 1 (powerA,powerB) is calibrated from CrazyFlie 2.1 data [OWD22];
% Type 2 is scaled for its higher mass/payload and should be refined from logs.
t1 = struct( ...
    'name',      'Scout (CrazyFlie 2.1)', ...
    'powerA',    3.3, ...     % quadratic power coefficient, W/(m/s)^2
    'powerB',    5.5, ...     % hover (static) power, W
    'vMin',      0.0, ...     % minimum speed, m/s
    'vMax',      1.0, ...     % maximum speed, m/s
    'aMin',     -0.5, ...     % minimum acceleration, m/s^2
    'aMax',      0.5, ...     % maximum acceleration, m/s^2
    'commRange', 8.0, ...     % effective coordination radius, m
    'battery',   250);        % battery capacity, mAh

t2 = struct( ...
    'name',      'Carrier/Relay (CrazyFlie 2.1 Brushless)', ...
    'powerA',    4.0, ...
    'powerB',    7.0, ...
    'vMin',      0.0, ...
    'vMax',      2.0, ...
    'aMin',     -1.0, ...
    'aMax',      1.0, ...
    'commRange', 12.0, ...
    'battery',   350);

params.types = [t1, t2];       % 1xK struct array; type index == agentType

% -------------------------------------------------------- controller gains
params.controller.alpha   = 2.5;   % attraction gain
params.controller.beta    = 50;    % repulsion gain
params.controller.gamma   = 0.2;   % repulsion range
params.controller.delta   = 20;    % tracking gain
params.controller.zeta    = 100;   % obstacle-avoidance gain
params.controller.psi     = 0.8;   % velocity-smoothing factor in [0,1]
params.controller.epsilon = 2.0;   % rendezvous tolerance, m

% potential weights (attraction, tracking, obstacle, energy)
params.controller.wAttr   = 0.2;
params.controller.wTrack  = 0.8;
params.controller.wObstc  = 18.0;
params.controller.wEnergy = 2.0;

% control frequencies and horizon
params.controller.freqOuter = 0.005;                 % rendezvous recompute rate, Hz
params.controller.freqInner = 10;                   % aggregation control rate, Hz
params.controller.simTime   = 2 / params.controller.freqOuter;  % total duration, s

% joint-connectivity window: the union of the communication graphs over any
% jointWindow consecutive control steps is assumed connected. The window in
% seconds is T = jointWindow / freqInner (here 5 s). The constant should be
% at least the worst-case (initial, most-dispersed) number of steps for
% separated components to link up; validate it empirically once the dynamics
% are run. It must stay well below the rendezvous period 1/freqOuter.
params.controller.jointWindow = 460;   % window length, control steps
params.controller.logGraph    = true; % log the per-step communication graph in
% simulateRound for the joint-connectivity
% check; set false to save memory for
% large N (disables the lambda2 metric)

% rendezvous / consensus detection
params.controller.consensusEps    = 5;  % squared-distance tolerance for consensus, m^2
params.controller.consensusConsec = 5;    % consecutive steps to declare consensus
params.controller.speedFactor     = 0.85; % speed cap factor applied before the bound

% --------------------------------------------------------- communication
params.comm.kappaC = 0.5;      % C-MAS packet size, Kbits
params.comm.kappaD = 1.0;      % D-MAS packet size (incl. region sharing), Kbits

% --------------------------------------------- region-of-interest & denied
% Polygons are defined on the original 10 m template and scaled to the arena,
% so the geometry is preserved when the arena size changes.
roiScale = params.arena.width / 10;
deniedScale = params.arena.width / 15;

roiTemplate = [ ...
    7.2556  8.5120;
    2.5237  7.8361;
    1.3872  6.8761;
    1.8878  1.5429;
    6.4391  3.6164;
    8.4342  5.2487;
    8.6665  6.8612;
    7.8229  7.9595];

deniedTemplate = [ ...
    4.5  2.2;
    5.0  5.0;
    7.2  4.5;
    7.5  3.0];

params.region.roi    = roiScale * roiTemplate;      % M x 2 vertices, m
params.region.denied = deniedScale * deniedTemplate;   % Md x 2 vertices, m
end
