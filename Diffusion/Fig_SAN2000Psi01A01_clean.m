%==========================================================================
%  Fig_SAN2000Psi01A01  —  Glutamate diffusion + stochastic astroglial uptake
%==========================================================================
%
%  PURPOSE
%  -------
%  Monte-Carlo random-walk simulation of glutamate molecules diffusing
%  through a model neuropil represented as a dense cloud of (possibly
%  overlapping) spheres. A fraction of those spheres are flagged as
%  astroglial processes; when a molecule enters an astroglial sphere it can
%  be transiently bound/taken up with a saturating probability. The script
%  produces:
%     * radial distance histograms of FREE and BOUND molecules vs. time
%       (DistanceFree*.txt / DistanceBound*.txt),
%     * an apparent diffusion coefficient (StatisticOut.txt, Diffusion_Set*.txt),
%     * the extracellular-space occupancy fraction Betta (and Alpha = 1-Betta),
%     * optional spatial snapshots of the molecular cloud (PD_*.txt).
%
%  UNITS
%  -----
%     Length : nanometres (nm)
%     Time   : milliseconds (ms)
%
%  DEPENDENCIES
%  ------------
%     InputParametersSR.m  (reads statisticSR.txt and returns the parameter
%                           structs used below).
%
%  REPRODUCIBILITY NOTE
%  --------------------
%  This is a faithful, professionally reformatted version of the original
%  script. ONLY comments, sectioning, and indentation have changed — every
%  executable statement, every arithmetic expression, and the exact order
%  and count of random draws (rand/randn) are preserved, so results are
%  numerically identical to the original. A few "quirks" are kept on purpose
%  and flagged inline (e.g. the literal 3.14 instead of pi, and a uint32
%  distance buffer); changing them WOULD change the numbers.
%
%  Original header: "new version with saturation of uptake glia, 01/02/2022"
%==========================================================================

clear;
clc;
rng('shuffle')   % Seed RNG from the system clock -> each run differs (by design)

%% ------------------------------------------------------------------------
%  1. LOAD PARAMETERS
%  ------------------------------------------------------------------------
%  All physical / numerical constants live in statisticSR.txt and are parsed
%  into structs by InputParametersSR().
[BasicParameters, Balls, ComputationSet, Cleft, SizeCube, Adhesion] = InputParametersSR();

DistanceBins        = 100;   % # radial bins for the distance histograms
TimeBins            = 100;   % # time bins over the simulation window
NonAdhesiveDistance = 50;    % nm — molecules closer than this to the centre
                             %      cannot be taken up (near the release site)

%% ------------------------------------------------------------------------
%  2. UNPACK FREQUENTLY USED PARAMETERS
%  ------------------------------------------------------------------------
np           = BasicParameters.NumberOfParticles;                 % # glutamate molecules released
charge       = BasicParameters.charge;                            % (kept for file headers; unused in model)
TotalTime    = BasicParameters.TotalTime;                         % ms — total simulated diffusion time
Numbersphere = BasicParameters.Numbersphere;                      % # spheres forming the neuropil
BeginTime    = BasicParameters.BeginTime * BasicParameters.TotalTime; % ms — start of D-fit window
EndTime      = BasicParameters.EndTime   * BasicParameters.TotalTime; % ms — end   of D-fit window
Trials       = BasicParameters.Trials;                            % # independent runs for statistics

fileIDStat = fopen('StatisticOut.txt','w');   % summary statistics output file

%% ------------------------------------------------------------------------
%  3. OCCUPANCY-PROBE AND GEOMETRY SETTINGS
%  ------------------------------------------------------------------------
NumberControlParticel = 100000;                 % # probe points used to estimate Betta (ECS fraction)
XTest = zeros(1, NumberControlParticel);        % (legacy preallocation; not used below)
DistanceX = SizeCube.Size;                      % nm — side length of the cubic simulation box

% Visualisation extents (only affects the live plot, not the physics)
AxisSize    = SizeCube.AxisSize;                % half-box size for the y/z axes
SizeOfaxisX = 1;                                % x-axis limit factor: Xlim = SizeOfaxisX * AxisSize

% Deterministic dwell time (in integration steps) if uptake were non-stochastic.
DigitalTimeofAdhesion = round(Adhesion.TimeOfDeAdhesion / ComputationSet.GTimeStep);

%% ------------------------------------------------------------------------
%  4. PER-MOLECULE STOCHASTIC DWELL TIMES
%  ------------------------------------------------------------------------
%  Each molecule draws a random "release/de-adhesion" countdown. With
%  probability ProbabilityofUptaken it gets 0 (i.e. it can be taken up
%  immediately); otherwise it is assigned a sampled dwell time from a
%  clipped Gaussian (mean 4000, std 2000, range [1000, 7000] steps).
r = rand(np,1);            % uniform draws, one per molecule (selects the branch below)
v = zeros(np,1);           % (legacy; unused)

ProbabilityofUptaken = 1 - BasicParameters.UnboundProb;   % 1 - P(stay unbound)

% --- Sampled dwell-time pool (Gaussian, then rounded and clipped) ---------
sample1 = randn(1, round(np));        % standard normal draws
sample1 = sample1 * 2000 + 4000;      % rescale: mean 4000, std 2000
sample1 = round(sample1);             % integer step counts
sample1(sample1 < 1000) = 1000;       % clip lower tail
sample1(sample1 > 7000) = 7000;       % clip upper tail

% (An optional second molecular population existed in the original and is
%  intentionally left disabled to preserve behaviour.)
% sample2 = randn(1, round(np/2)); ... (commented out)

sample = sample1;

% Assign each molecule its dwell time (0 = immediately uptakable).
DigitalTimeofAdhesionStochastic = zeros(np,1);
for i = 1:np
    if r(i) < ProbabilityofUptaken
        DigitalTimeofAdhesionStochastic(i) = 0;
    else
        DigitalTimeofAdhesionStochastic(i) = sample(i);
    end
end

%% ------------------------------------------------------------------------
%  5. TIME-GRID AND ACCUMULATORS
%  ------------------------------------------------------------------------
NumberTimeInteraction = round(TotalTime / ComputationSet.GTimeStep); % total integration steps
Betta = zeros(1, Trials);   % per-trial extracellular-space occupancy fraction

%% ========================================================================
%  6. MAIN TRIAL LOOP
%  ========================================================================
for IterTrial = 1:Trials   % repeat the whole experiment for statistics

    CreatMatrixAverage = 0.01;   % running marker controlling histogram-column placement
    IndexMatrix        = 1;      % current column index in the average matrices
    Time = 0:ComputationSet.GTimeStep:TotalTime;

    figure(1)

    % --- Initialise sphere arrays --------------------------------------
    [x,y,z] = sphere;                 % unit-sphere mesh (for plotting only)
    r  = zeros(Numbersphere,2);       % col 1 = radius (nm), col 2 = astroglia flag (1/0)
    xd = zeros(Numbersphere,1);       % sphere centre x (nm)
    yd = zeros(Numbersphere,1);       % sphere centre y (nm)
    zd = zeros(Numbersphere,1);       % sphere centre z (nm)

    % --- Diffusion-coefficient accumulators ----------------------------
    DifCoefTotal  = zeros(1);
    DifCoefTotalX = zeros(1);

    TotalNumberinteration = 1;        % (outer repeat kept at 1, as in original)

    % ====================================================================
    %  6a. GEOMETRY REALISATION (one pass: TotalNumberinteration == 1)
    % ====================================================================
    for NumberIteretion = 1:TotalNumberinteration
        TempTime = NumberIteretion * ComputationSet.GTimeStep;

        % --- Sphere #1: a tiny seed sphere placed away from the release
        %     point so it does not block molecules at the origin. -------
        r(1) = 3; r(1,2) = 0;   % radius = 3 nm, flagged non-astroglial
        xd(1) = 1000;           % centre at (1000,1000,1000) nm
        yd(1) = 1000;
        zd(1) = 1000;
        surf((r(1)*x+xd(1)), (r(1)*y+yd(1)), (r(1)*z+zd(1)))   % draw seed sphere
        hold on

        % --- Per-molecule state matrix "Stak" ------------------------------
        %     col 1 : uptake/binding counter (1 = free; >1 = bound, counting up)
        %     col 2 : (auxiliary marker)
        %     col 3 : "alive/released" flag (>0 means the molecule is active)
        Stak = ones(np, 1, 1);
        Stak(:,2) = 10;
        Stak(:,3) = 10;

        % --- Molecule position buffers (grow in time; col = molecule) ------
        x_par = zeros(1,np);
        y_par = zeros(1,np);
        z_par = zeros(1,np);

        TimeOfAdhesive = zeros(1, np);   % per-molecule time spent inside an astroglial sphere

        % --- Scalar accumulators for this geometry pass --------------------
        Distance      = zeros(1);
        TempDistance  = 0;
        TempDistanceX = 0;
        TempDistanceY = 0;
        TempDistanceZ = 0;

        DifCoef        = zeros(1);
        NumberParticles= zeros(1);
        DifCoefX       = zeros(1);
        DifCoefY       = zeros(1);
        DifCoefZ       = zeros(1);

        ParameterControl = 0.5;          % min inter-sphere gap reference (nm)
        ControlParam     = 0;
        SpaceBalls       = (4/3)*3.14* r(1,1)^3;   % running occupied volume (NOTE: literal 3.14 kept on purpose)

        % ----------------------------------------------------------------
        %  6b. POPULATE THE BOX WITH SPHERES (spheres 2..Numbersphere)
        % ----------------------------------------------------------------
        for j = 2:Numbersphere
            test = j;

            if strcmp(Balls.OverlapingBalls,'y')
                % ---- Overlapping spheres allowed: just keep the cleft clear ----
                while 1
                    DistanceBetweenBalls = 2*ParameterControl + zeros(1,j,'uint32'); % (uint32 buffer kept as-is)
                    xd(j) = DistanceX*(rand()-0.5);
                    yd(j) = DistanceX*(rand()-0.5);
                    zd(j) = DistanceX*(rand()-0.5);
                    r(j,1) = (Balls.MinBallRadius + rand()*Balls.MaxBallRadius);
                    if rand() < Balls.ProbabilityofAstrocytes
                        r(j,2) = 1;     % astroglial process
                    else
                        r(j,2) = 0;     % neuronal/other process
                    end
                    % Accept only if the sphere stays outside the synaptic cleft margin.
                    if sqrt(xd(j)^2+yd(j)^2+zd(j)^2) - Cleft.RadiusofCleft - r(j,1) >= Cleft.DistanceAstrocyte
                        break
                    end
                end
            else
                % ---- Non-overlapping spheres: also enforce min pairwise gap ----
                while 1
                    % Inner loop: propose a sphere outside the cleft margin
                    while 1
                        DistanceBetweenBalls = 2*ParameterControl + zeros(1,j,'uint32');
                        xd(j) = DistanceX*(rand()-0.5);
                        yd(j) = DistanceX*(rand()-0.5);
                        zd(j) = DistanceX*(rand()-0.5);
                        r(j,1) = (Balls.MinBallRadius + rand()*Balls.MaxBallRadius);
                        if rand() < Balls.ProbabilityofAstrocytes
                            r(j,2) = 1;
                        else
                            r(j,2) = 0;
                        end
                        if sqrt(xd(j)^2+yd(j)^2+zd(j)^2) - Cleft.RadiusofCleft - r(j,1) >= Cleft.DistanceAstrocyte
                            break
                        end
                    end

                    % Check pairwise separation against all previously placed spheres.
                    % NOTE: r(j) uses linear indexing into the 2-col matrix (kept as-is);
                    % the uint32 buffer above clamps negative gaps to 0.
                    for CheakIter = 1:j-1
                        DistanceBetweenBalls(CheakIter) = ...
                            sqrt((xd(j)-xd(CheakIter))^2 + (yd(j)-yd(CheakIter))^2 + (zd(j)-zd(CheakIter))^2) ...
                            - (r(j)+r(CheakIter));
                    end
                    if min(DistanceBetweenBalls) >= 0.1   % 0.1 nm minimum gap
                        break
                    end
                end
            end

            % Accumulate occupied volume (literal 3.14 kept on purpose).
            SpaceBalls = SpaceBalls + (4/3)*3.14* r(j,1)^3;

            % --- Optional per-sphere rendering (disabled for speed) --------
            % Drawing every sphere is correct but dramatically slows the run.
            % surf((r(j,1)*x+xd(j)), (r(j,1)*y+yd(j)), (r(j,1)*z+zd(j)), ...
            %      'FaceColor','red','EdgeColor','none')
            % xlim([-SizeOfaxisX*AxisSize SizeOfaxisX*AxisSize])
            % ylim([-AxisSize AxisSize]); zlim([-AxisSize AxisSize])
            % colormap([0 1 0])
        end % sphere-placement loop
        hold on

        % --- Optionally save the sphere distribution to disk ---------------
        if strcmp(BasicParameters.SaveFiles,'y')
            Filename = sprintf('Balls distribution %s.txt', ...
                datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-ss'));
            writematrix([round(xd), round(yd), round(zd), round(r)], Filename);
        end

        % --- Total filled fraction of the box (excluding the cleft volume) -
        FilledSpace = SpaceBalls / ((DistanceX)^3 - (4/3)*3.14*Cleft.RadiusofCleft^3);

        % ----------------------------------------------------------------
        %  6c. ESTIMATE EXTRACELLULAR-SPACE FRACTION (Betta) BY PROBING
        % ----------------------------------------------------------------
        %  Scatter NumberControlParticel probe points outside the cleft; a
        %  point is "inside tissue" (-1) if it falls within ANY sphere.
        PositivePatricle = ones(1, NumberControlParticel);
        xTest = linspace(0,0,NumberControlParticel);
        yTest = linspace(0,0,NumberControlParticel);
        zTest = linspace(0,0,NumberControlParticel);

        for ic = 1:NumberControlParticel
            % Draw a probe point outside the cleft margin.
            while 1
                xTest(ic) = DistanceX*(rand()-0.5);
                yTest(ic) = DistanceX*(rand()-0.5);
                zTest(ic) = DistanceX*(rand()-0.5);
                if sqrt(xTest(ic)^2+yTest(ic)^2+zTest(ic)^2) - Cleft.RadiusofCleft > Cleft.DistanceAstrocyte
                    break
                end
            end

            % Flag as inside-tissue if it lies within any placed sphere.
            for j = 2:Numbersphere
                VectorControl = sqrt((xTest(ic)-xd(j))^2 + (yTest(ic)-yd(j))^2 + (zTest(ic)-zd(j))^2) - r(j,1);
                if VectorControl < 0
                    PositivePatricle(ic) = -1;
                end
            end
        end

        % Betta = fraction of probe points inside tissue; Alpha = ECS fraction.
        InPar = sum(PositivePatricle(:) < 0);
        Betta(IterTrial) = InPar / NumberControlParticel;
        AlphaSpace       = 1 - Betta(IterTrial);

        % ----------------------------------------------------------------
        %  6d. PREPARE THE LIVE PLOT AND HISTOGRAM MATRICES
        % ----------------------------------------------------------------
        daspect([1 1 1])
        p_plot = scatter3(x_par(1,:), y_par(1,:), z_par(1,:), 20, 'filled');
        axis([-SizeOfaxisX*AxisSize SizeOfaxisX*AxisSize  -AxisSize AxisSize  -AxisSize AxisSize]);

        TimeIter = 0;
        InputPar = [np; Numbersphere];
        AverageTimeMatrix     = zeros(DistanceBins, TimeBins+1);  % BOUND-molecule distance histograms vs time
        AverageTimeMatrixFree = zeros(DistanceBins, TimeBins+1);  % FREE-molecule  distance histograms vs time

        % --- Snapshot schedule: capture the cloud at a few target times ----
        SnapshotTimes = [0.0125, 0.0375, 0.125, 0.375];  % fractions of TotalTime (~0.1,0.3,1,3 ms at 8 ms)
        t_targets_ms  = SnapshotTimes * TotalTime;        % target times (ms)
        saved         = false(size(t_targets_ms));        % one-shot save flags

        dt_ms  = ComputationSet.GTimeStep;                % integration step (ms)
        t_prev = 0;                                       % previous step's time (ms)

        % ================================================================
        %  6e. TIME INTEGRATION (random walk + uptake)
        % ================================================================
        for i = 1:NumberTimeInteraction

            % Trigger a (second) release event at t = 10 steps' worth of ms.
            if i == round(10/ComputationSet.GTimeStep)
                Stak(:,3) = 8;
            end

            TimeIter = TimeIter + 1;

            % ---- Periodic logging + radial-histogram accumulation --------
            if strcmp(BasicParameters.SaveFiles,'y')
                TimeOutPercent = 1/TimeBins;
                if TimeIter == round(TimeOutPercent*TotalTime/ComputationSet.GTimeStep)
                    TimeofDiffusion = i * ComputationSet.GTimeStep;
                    formatSpec = 'Time of diffusion  %8.4f  ms and total time is  %8.4f ms\n';
                    fprintf(formatSpec, TimeofDiffusion, TotalTime);
                    TimeIter = 0;

                    if (i == round(CreatMatrixAverage*TotalTime/ComputationSet.GTimeStep))
                        % Radial distance of each molecule from the origin.
                        SQBound = zeros(1, BasicParameters.NumberOfParticles);
                        SQFree  = zeros(1, BasicParameters.NumberOfParticles);
                        TempNew = (x_par(i,:).^2 + y_par(i,:).^2 + z_par(i,:).^2).^0.5;
                        TempNew(TempNew < Cleft.RadiusofCleft) = 0;   % ignore inside-cleft positions
                        ttt = Stak(:,2);

                        % Split into BOUND (Stak>1) and FREE (Stak<=1) populations.
                        SQBound(Stak(:,1) > 1)   = TempNew(Stak(:,1) > 1);
                        SQFree(Stak(:,1) < 1.1)  = TempNew(Stak(:,1) < 1.1);

                        % Radial bin edges (nm).
                        DistanceBinsNew = (Cleft.RadiusofCleft + round(DistanceX/DistanceBins)) ...
                                          : round(DistanceX/DistanceBins) ...
                                          : (DistanceX + Cleft.RadiusofCleft + round(DistanceX/DistanceBins));
                        TimeBinsNew = 0:1:TimeBins;

                        % Histogram both populations over distance.
                        [NUmm,     HSQ]     = histcounts(SQBound, DistanceBinsNew);
                        HSQ = HSQ(2:end);
                        [NUmmFree, HSQFree] = histcounts(SQFree,  DistanceBinsNew);
                        HSQFree = HSQFree(2:end);

                        % First valid time -> store bin centres in col 1 and counts in col 2;
                        % afterwards append a new counts column per logged time.
                        if CreatMatrixAverage < TimeOutPercent + 0.001
                            AverageTimeMatrix(1:end,1)     = HSQ/1000;     % bin centres (um)
                            AverageTimeMatrix(1:end,2)     = NUmm;         % bound counts
                            AverageTimeMatrixFree(1:end,1) = HSQFree/1000; % bin centres (um)
                            AverageTimeMatrixFree(1:end,2) = NUmmFree;     % free counts
                        else
                            IndexMatrix = IndexMatrix + 1;
                            AverageTimeMatrix(1:end,IndexMatrix+1)     = NUmm;
                            AverageTimeMatrixFree(1:end,IndexMatrix+1) = NUmmFree;
                        end
                        CreatMatrixAverage = CreatMatrixAverage + TimeOutPercent;
                    end
                end
            end % logging / histogram block

            TempNumber = 0;   % # molecules counted toward the D estimate this step

            % ------------------------------------------------------------
            %  Per-molecule update
            % ------------------------------------------------------------
            for j = 1:np

                if Stak(j,3) > 0   % only active (released) molecules move
                    TempPosition = sqrt(x_par(i,j)^2 + y_par(i,j)^2 + z_par(i,j)^2);

                    if TempPosition > Cleft.RadiusofCleft
                        % ---- FREE-SPACE diffusion (outside the cleft) -------
                        if Stak(j,1) == 1
                            % Free molecule: full random-walk step in 3D.
                            x_par(i+1,j) = x_par(i,j) + ComputationSet.StepRad*2*(rand()-0.5);
                            y_par(i+1,j) = y_par(i,j) + ComputationSet.StepRad*2*(rand()-0.5);
                            z_par(i+1,j) = z_par(i,j) + ComputationSet.StepRad*2*(rand()-0.5);
                        else
                            % Bound molecule: frozen in place this step.
                            x_par(i+1,j) = x_par(i,j);
                            y_par(i+1,j) = y_par(i,j);
                            z_par(i+1,j) = z_par(i,j);
                        end

                        % ---- Sphere interaction (uptake + soft reflection) --
                        for j_Sphere = 1:Numbersphere
                            % (Optional decay of sticky spheres over a timescale
                            %  is available in the original but left disabled.)

                            % If the proposed position lands inside a sphere:
                            if sqrt((x_par(i+1,j)-xd(j_Sphere))^2 + ...
                                    (y_par(i+1,j)-yd(j_Sphere))^2 + ...
                                    (z_par(i+1,j)-zd(j_Sphere))^2) < r(j_Sphere,1)

                                % Possible astroglial uptake (only away from the release site).
                                if TempPosition > NonAdhesiveDistance
                                    if r(j_Sphere,2) == 1
                                        % Saturating binding probability:
                                        %   P = MaxProb * (1 - exp(-t_in / Psi))
                                        TimeOfAdhesive(j) = TimeOfAdhesive(j) + 1;
                                        if Adhesion.MaxProbAdhesive * ...
                                           (1 - exp(-(ComputationSet.GTimeStep*TimeOfAdhesive(j))/Adhesion.TimeINsideAdhesiveZone)) > rand()
                                            Stak(j,1) = Stak(j,1) + 1;   % become bound
                                            TimeOfAdhesive(j) = 0;
                                        end
                                    end
                                end

                                % Soft local jitter instead of hard reflection:
                                % the molecule "waits" near the sphere until it
                                % randomly steps clear (avoids costly multi-sphere
                                % reflection checks in dense packings).
                                x_par(i+1,j) = x_par(i,j) + ComputationSet.StepLocalMovement*ComputationSet.StepRad*2*(rand()-0.5);
                                y_par(i+1,j) = y_par(i,j) + ComputationSet.StepLocalMovement*ComputationSet.StepRad*2*(rand()-0.5);
                                z_par(i+1,j) = z_par(i,j) + ComputationSet.StepLocalMovement*ComputationSet.StepRad*2*(rand()-0.5);
                            end
                        end % sphere-interaction loop

                        % Advance the bound-state counter and release after the
                        % per-molecule stochastic dwell time elapses.
                        if Stak(j,1) > 1
                            Stak(j,1) = Stak(j,1) + 1;
                        end
                        if Stak(j,1) == DigitalTimeofAdhesionStochastic(j)
                            Stak(j,1) = 1;   % released back to free state
                        end

                        % Open-boundary rule: only molecules still inside the box
                        % contribute to the mean-squared-displacement / D estimate.
                        if (sqrt((x_par(i+1,j))^2) < DistanceX/2) || ...
                           (sqrt((y_par(i+1,j))^2) < DistanceX/2) || ...
                           (sqrt((z_par(i+1,j))^2) < DistanceX/2)
                            TempDistance  = TempDistance  + (x_par(i+1,j)^2 + y_par(i+1,j)^2 + z_par(i+1,j)^2)/(6*i*ComputationSet.GTimeStep);
                            TempDistanceX = TempDistanceX + (x_par(i+1,j)^2)/(2*i*ComputationSet.GTimeStep);
                            TempNumber    = TempNumber + 1;
                        end

                    else
                        % ---- INSIDE-CLEFT diffusion (reduced/confined steps) -
                        if (sqrt((x_par(i,j))^2) < Cleft.RadiusofCleft)
                            x_par(i+1,j) = x_par(i,j) + ComputationSet.StepRad*2*(rand()-0.5);
                        else
                            x_par(i+1,j) = x_par(i,j);
                        end
                        if (sqrt((y_par(i,j))^2) < Cleft.RadiusofCleft)
                            y_par(i+1,j) = y_par(i,j) + ComputationSet.StepRad*2*(rand()-0.5);
                        else
                            y_par(i+1,j) = y_par(i,j);
                        end
                        if (sqrt((z_par(i,j))^2) < Cleft.CleftWidth)
                            % Thin cleft: small z-jitter only.
                            z_par(i+1,j) = z_par(i,j) + 0.1*ComputationSet.StepRad*2*(rand()-0.5);
                        else
                            % Outside the cleft half-width: nudge back toward the gap.
                            if z_par(i,j) < 0
                                z_par(i+1,j) = z_par(i,j) + 0.1*Cleft.CleftWidth;
                            end
                            if z_par(i,j) > 0
                                z_par(i+1,j) = z_par(i,j) - 0.1*Cleft.CleftWidth;
                            end
                        end

                        % Same open-boundary contribution rule as above.
                        if (sqrt((x_par(i+1,j))^2) < DistanceX/2) || ...
                           (sqrt((y_par(i+1,j))^2) < DistanceX/2) || ...
                           (sqrt((z_par(i+1,j))^2) < DistanceX/2)
                            TempDistance  = TempDistance  + (x_par(i+1,j)^2 + y_par(i+1,j)^2 + z_par(i+1,j)^2)/(6*i*ComputationSet.GTimeStep);
                            TempDistanceX = TempDistanceX + (x_par(i+1,j)^2)/(2*i*ComputationSet.GTimeStep);
                            TempNumber    = TempNumber + 1;
                        end
                    end

                    % "Do not fly out of the pipe" guard (currently a no-op:
                    %  the corrective lines are intentionally commented out).
                    if (sqrt(x_par(i,j)^2 + y_par(i,j)^2 + z_par(i,j)^2)) > SizeCube.Size/2
                        % x_par(i+1,j)=x_par(i-5,j); ... (disabled)
                    end

                else
                    % Inactive molecule: hold position.
                    x_par(i+1,j) = x_par(i,j);
                    y_par(i+1,j) = y_par(i,j);
                    z_par(i+1,j) = z_par(i,j);
                end

                % Far-escaped molecules are parked even further away so they
                % no longer interact with the populated region.
                if (sqrt(x_par(i,j)^2 + y_par(i,j)^2 + z_par(i,j)^2)) > 4500
                    x_par(i+1,j) = x_par(i,j) + 5000;
                    y_par(i+1,j) = y_par(i,j) + 5000;
                    z_par(i+1,j) = z_par(i,j) + 5000;
                end
            end % per-molecule loop

            % ------------------------------------------------------------
            %  Snapshot capture (time-crossing detection): save the cloud
            %  once, the first step at/after each target time.
            % ------------------------------------------------------------
            t_curr = i * dt_ms;
            to_save = find(~saved & (t_prev < t_targets_ms) & (t_targets_ms <= t_curr));

            for k = 1:numel(to_save)
                ss   = to_save(k);
                t_ms = t_targets_ms(ss);   % label by target time

                % Encode binding state: bound -> 0, free -> 1.
                StakNew = Stak(:,1);
                StakNew(StakNew > 1) = 0;
                A = [x_par(i,:); y_par(i,:); z_par(i,:); StakNew'];

                Filename = sprintf('PD_%.4fms_%s.txt', t_ms, ...
                    datetime('now','TimeZone','local','Format','yyyy-MM-dd_HH-mm-ss'));

                fileID = fopen(Filename, 'w');
                if fileID < 0
                    error('File write failed: %s', Filename);
                end
                fprintf(fileID, '%10.5f %10.5f %10.5f %10.0f \n', A);
                fclose(fileID);

                saved(ss) = true;   % mark as captured
            end

            t_prev = t_curr;
            % --- end snapshot capture ---

            % --- Update the live scatter plot ------------------------------
            set(p_plot,'XData',x_par(i,:));
            set(p_plot,'YData',y_par(i,:));
            set(p_plot,'ZData',z_par(i,:));

            % --- Per-step diffusion-coefficient bookkeeping ----------------
            NumberParticles(i+1) = TempNumber;
            TempNumberParticla   = TempNumber;
            DifCoef(i+1)  = (TempDistance  / TempNumber)/10^6;   % isotropic D estimate (um^2/ms)
            DifCoefX(i+1) = (TempDistanceX / TempNumber)/10^6;   % x-only D estimate
            TempNumber    = 0;
            TempDistance  = 0;
            TempDistanceX = 0;
            pause(0.0001);   % yield to the graphics event queue (keeps the plot live)

            % --- Optional dump of the final particle distribution ----------
            if (TempNumberParticla < np) || (i == NumberTimeInteraction - 1)
                if (BasicParameters.Control == 1)
                    A = [x_par(i,:); y_par(i,:); z_par(i,:)];
                    InputParameters = [charge; Numbersphere];
                    fileID = fopen('FinalDistribution.txt','w');
                    fprintf(fileID,'%10.5f %10.5f\n', InputParameters);
                    fprintf(fileID,'%10.5f %10.5f %10.5f\n', A);
                    Control = 2;
                    fclose(fileID);
                end
            end

        end % time-integration loop

        % Accumulate this geometry pass into the trial totals.
        DifCoefTotal  = DifCoefTotal  + DifCoef;
        DifCoefTotalX = DifCoefTotalX + DifCoefX;

    end % geometry-realisation loop (NumberIteretion)

    % --------------------------------------------------------------------
    %  6f. SAVE PER-TRIAL DISTANCE HISTOGRAMS
    % --------------------------------------------------------------------
    Filename = sprintf('DistanceBound %s.txt', ...
        datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-ss'));
    fileID = fopen(Filename,'w');
    writematrix(AverageTimeMatrix, Filename)
    fclose(fileID);

    Filename = sprintf('DistanceFree %s.txt', ...
        datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-ss'));
    fileID = fopen(Filename,'w');
    writematrix(AverageTimeMatrixFree, Filename)
    fclose(fileID);

end % MAIN TRIAL LOOP

%% ========================================================================
%  7. APPARENT DIFFUSION COEFFICIENT OVER THE FIT WINDOW
%  ========================================================================
Begin = fix(BeginTime/ComputationSet.GTimeStep);
End   = fix(EndTime/ComputationSet.GTimeStep);
MeanDiffusion = mean(DifCoefTotalX(Begin : End));
STDDiffusion  = std (DifCoefTotalX(Begin : End));

figure(2)
% (Time-series plots of distance / D were used during development; left
%  commented to match the original output behaviour.)

% --- Write summary statistics (mean D, std D, last Betta) ----------------
DiffResults = [MeanDiffusion; STDDiffusion];
fprintf(fileIDStat, '%10.5f %10.5f %10.5f\n', DiffResults, Betta(IterTrial));
fclose(fileIDStat);

% --- Write the full D(t) trace plus particle counts ----------------------
A = [Time; DifCoefTotalX/TotalNumberinteration; NumberParticles];
InputParameters = [charge; Numbersphere];
Filename = sprintf('Diffusion_Set_ms_%s.txt', ...
    datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-ss'));
fileID = fopen(Filename,'w');
fprintf(fileID,'%10.5f %10.5f\n', InputParameters);
fprintf(fileID,'%10.5f %10.5f %10.5f\n', A);
fclose(fileID);

%==========================================================================
%  End of script
%==========================================================================
