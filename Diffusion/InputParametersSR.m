function [BasicParameters, Balls, ComputationSet, Cleft, SizeCube, Adhesion] = InputParametersSR()
%INPUTPARAMETERSSR
% -------------------------------------------------------------------------
%
% Input file
%   - Reads key–value pairs from 'statisticSR.txt' using the format:
%       key = value
%     Lines starting with '%' are treated as comments.
%
% Output structs
%   - BasicParameters: general run controls and counters
%   - Balls:           geometry and distribution of spherical obstacles
%   - ComputationSet:  diffusion constants and integration step sizes
%   - Cleft:           synaptic cleft geometry
%   - SizeCube:        simulation box size
%   - Adhesion:        kinetics/probabilities for astroglial binding
%
% Notes
%   - BasicParameters.BeginTime / EndTime are kept as FRACTIONS of
%     BasicParameters.TotalTime (as in your original scripts).
%   - SaveFiles remains a character flag ('y' to enable), for compatibility.
% -------------------------------------------------------------------------

    % ---- Read parameter file ------------------------------------------------
    fname = 'statisticSR.txt';
    fid = fopen(fname, 'r');
    if fid < 0
        error('InputParametersSR:FileNotFound', 'Cannot open "%s".', fname);
    end
    s = textscan(fid, '%s%f', 'CommentStyle', '%', 'Delimiter', '=');
    fclose(fid);

    if isempty(s) || isempty(s{1})
        error('InputParametersSR:EmptyFile', 'No readable key–value pairs found in "%s".', fname);
    end

    % Convenience masks (case-insensitive matching of keys)
    idxnp                         = strncmpi('np', s{1}, 2);
    idxTrials                     = strncmpi('Trials', s{1}, 6);
    indMaxProbAdhesive            = strncmpi('MaxProbAdhesive', s{1}, 15);
    indTimeINsideAdhesiveZone     = strncmpi('TimeINsideAdhesiveZone', s{1}, 23);
    indProbabilityofAstrocytes    = strncmpi('ProbabilityofAstrocytes', s{1}, 24);
    indNumbersphere               = strncmpi('Numbersphere', s{1}, 12);
    indUnboundProb                = strncmpi('UnboundProb', s{1}, 11);

    % ---- Parse numeric parameters from file --------------------------------
    % Required keys (will error if missing)
    try
        np                       = s{2}(idxnp);                      % # of particles
        Trials                   = s{2}(idxTrials);                  % # of trials (runs)
        MaxProbAdhesive          = s{2}(indMaxProbAdhesive);         % cap probability for binding increment
        TimeINsideAdhesiveZone   = s{2}(indTimeINsideAdhesiveZone);  % ms time constant in binding kernel
        ProbabilityofAstrocytes  = s{2}(indProbabilityofAstrocytes); % fraction of spheres flagged as astroglial
        Numbersphere             = s{2}(indNumbersphere);            % # of spheres to place
        UnboundProb              = s{2}(indUnboundProb);             % probability to remain unbound
    catch
        error('InputParametersSR:MissingKeys', ...
              'One or more required keys are missing from "%s".', fname);
    end

    % ---- Light-touch validation (warn rather than assert for robustness) ----
    warnIf(~isscalar(np) || np<=0 || mod(np,1)~=0, ...
    'np should be a positive integer.');

warnIf(~isscalar(Trials) || Trials<=0 || mod(Trials,1)~=0, ...
    'Trials should be a positive integer.');

warnIf(any(MaxProbAdhesive<0 | MaxProbAdhesive>1), ...
    'MaxProbAdhesive should be in [0,1].');

warnIf(any(UnboundProb<0 | UnboundProb>1), ...
    'UnboundProb should be in [0,1].');

warnIf(any(ProbabilityofAstrocytes<0 | ProbabilityofAstrocytes>1), ...
    'ProbabilityofAstrocytes should be in [0,1].');

warnIf(~isscalar(Numbersphere) || Numbersphere<1 || mod(Numbersphere,1)~=0, ...
    'Numbersphere should be a positive integer.');

warnIf(TimeINsideAdhesiveZone<0, ...
    'TimeINsideAdhesiveZone should be >= 0 ms.');

    % -------------------------------------------------------------------------
    % GEOMETRY (nanometres)
    % -------------------------------------------------------------------------
    SizeCube.Size     = 4000;         % Simulation box side: 4 µm = 4000 nm
    SizeCube.AxisSize = SizeCube.Size/2;

    % -------------------------------------------------------------------------
    % BASIC PARAMETERS
    % -------------------------------------------------------------------------
    BasicParameters.SaveFiles        = 'y';   % 'y' → write auxiliary outputs (kept for compatibility)
    BasicParameters.NumberOfParticles= np;    % particles released
    BasicParameters.charge           = 0;     % not used in current model
    BasicParameters.TotalTime        = 8;     % ms, total simulation duration
    BasicParameters.Numbersphere     = Numbersphere;
    BasicParameters.UnboundProb      = UnboundProb;

    % These are interpreted as FRACTIONS of TotalTime elsewhere in your code:
    BasicParameters.BeginTime        = 0.8;   % fraction of TotalTime (e.g., 0.8*8 ms)
    BasicParameters.EndTime          = 0.9;   % fraction of TotalTime
    BasicParameters.Trials           = Trials;
    BasicParameters.Control          = 2;     % if 1 → write final particle distribution

    % -------------------------------------------------------------------------
    % SPHERES ("Balls"): neuronal & astroglial processes (overlapping)
    % Radii are in nanometres, matching EM-observed cross-sections.
    % -------------------------------------------------------------------------
    Balls.name                   = 'Balls';
    Balls.MinBallRadius          = 50.0;     % nm (≈ 50–300 nm range for cortical neuropil)
    Balls.MaxBallRadius          = 350.0;    % nm
    Balls.BettaControl           = 0.8;      % target occupancy β ≈ 0.8 (i.e., α ≈ 0.2 ECS)
    Balls.OverlapingBalls        = 'y';      % allow overlapping spheres
    Balls.ProbabilityofAstrocytes= ProbabilityofAstrocytes; % fraction flagged as astroglial

    % Guard against swapped radii in the input
    if Balls.MaxBallRadius < Balls.MinBallRadius
        warning('MaxBallRadius < MinBallRadius; swapping values.');
        tmp = Balls.MaxBallRadius;
        Balls.MaxBallRadius = Balls.MinBallRadius;
        Balls.MinBallRadius = tmp;
    end

    % -------------------------------------------------------------------------
    % DIFFUSION / INTEGRATION
    %   D_free = 0.5 µm^2/ms = 5e5 nm^2/ms
    %   dt     = 0.0001 ms  (10× smaller than the original 0.001 ms)
    %   StepRad = sqrt(6 * D * dt)  → nm per step
    % -------------------------------------------------------------------------
    ComputationSet.DifCoefINIT     = 5e5;        % nm^2/ms (free-space glutamate diffusivity)
    ComputationSet.GTimeStep       = 0.0001;     % ms  (↓ 10× integration step as requested)
    ComputationSet.StepRad         = sqrt(6 * ComputationSet.DifCoefINIT * ComputationSet.GTimeStep); % nm
    ComputationSet.StepLocalMovement= 0.1;       % dimensionless jitter factor near surfaces (×StepRad)

    % -------------------------------------------------------------------------
    % SYNAPTIC CLEFT GEOMETRY (nm)
    %   220-nm diameter hemispheres around a 20-nm cleft (centre-aligned)
    %   RadiusofCleft refers to hemispherical radius from the centre.
    % -------------------------------------------------------------------------
    Cleft.RadiusofCleft     = 100;    % nm (≈ 200 nm diameter; adjust to 110 for exactly 220 nm)
    Cleft.CleftWidth        = 10;     % nm (half-height; total gap ≈ 20 nm)
    Cleft.DistanceAstrocyte = 10;     % nm (min distance from cleft margin to sphere centres)

    % -------------------------------------------------------------------------
    % ASTROGLIAL BINDING ("adhesion") KINETICS
    %   - TimeOfDeAdhesion: typical immobilisation window (ms) if used deterministically
    %   - MaxProbAdhesive:  cap for per-step binding probability increment
    %   - TimeINsideAdhesiveZone: time constant Ψ (ms) in P(t)=1-exp(-t/Ψ)
    % -------------------------------------------------------------------------
    Adhesion.TimeOfDeAdhesion       = 4;                 % ms (often superseded by stochastic dwell)
    Adhesion.MaxProbAdhesive        = MaxProbAdhesive;   % 0..1
    Adhesion.TimeINsideAdhesiveZone = TimeINsideAdhesiveZone; % ms

    % -------------------------------------------------------------------------
    % Helper: warning wrapper (keeps the main body clean)
    % -------------------------------------------------------------------------
    function warnIf(cond, msg)
        if cond
            warning(['InputParametersSR:Validation: ' msg]);
        end
    end
end
