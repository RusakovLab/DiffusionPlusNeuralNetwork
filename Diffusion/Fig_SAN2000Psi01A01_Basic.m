% new version with saturayion of uptake glia
% 01/02/2022
%*******************************************

clear;
clc;
rng('shuffle') % set seeds the random number generator based on the current time

% Convert the input parameters to the numerical value
[BasicParameters, Balls, ComputationSet, Cleft, SizeCube, Adhesion]=InputParametersSR();
DistanceBins = 100;
TimeBins = 100;
NonAdhesiveDistance = 50; % nm

%********************************************************
% InputParameters

np = BasicParameters.NumberOfParticles;              % Number of p
charge = BasicParameters.charge;             % Charge of particles
TotalTime = BasicParameters.TotalTime;                  % Time of computation in ms 0.21
Numbersphere = BasicParameters.Numbersphere;          % b=0.3 13000 b=0.8 40000
BeginTime = BasicParameters.BeginTime*BasicParameters.TotalTime;       % Offset  to compute the diffusion coefficient in ms 0.015
EndTime = BasicParameters.EndTime*BasicParameters.TotalTime;         % End time to compute the diffusion coefficient in ms 0.02
Trials = BasicParameters.Trials;                  % Number of Runs to compute the statistics
fileIDStat = fopen('StatisticOut.txt','w');       % The name of file to save the statistical results


%**************************************************************************
NumberControlParticel=100000; % The number of particles to control the  space for diffusion, parameter betta
XTest=zeros(1,NumberControlParticel);
DistanceX = SizeCube.Size;  % in nm, the size modulo the spatial region, x-coordinate, in which the balls are placed./ // Betta 0.3 1700

%***********************************************
% Size of visualization area
AxisSize=SizeCube.AxisSize; %maximum size axis y and z for visualization of diffusion in micrometers
SizeOfaxisX = 1; % Xlim= SizeOfaxisX * AxisSize

% *************************************************************************
% Adhesive parameters
DigitalTimeofAdhesion = round(Adhesion.TimeOfDeAdhesion/ComputationSet.GTimeStep);


r = rand(np,1);
v = zeros(np,1);
%**************************************************************************
ProbabilityofUptaken = 1-BasicParameters.UnboundProb; %0.5; 1 - is NOT UNBOUND

%**************************************************************************
% Set the random number generator seed for reproducibility
% rng(1);
% Generate a normally distributed sample with mean 0 and standard deviation 1
sample1 = randn(1, round(np));
% Scale the sample to have mean 4000 and standard deviation 2000
sample1 = sample1 * 2000 + 4000;
% Round the sample to the nearest integer
sample1 = round(sample1);
% Clip the sample to the desired range of 1000 to 7000
sample1(sample1 < 1000) = 1000;
sample1(sample1 > 7000) = 7000;
%**************************************************************************

%**************************************************************************
% % Set the random number generator seed for reproducibility
% % rng(1);
% % Generate a normally distributed sample with mean 0 and standard deviation 1
% sample2 = randn(1, round(np/2));
% % Scale the sample to have mean 4000 and standard deviation 2000
% sample2 = sample2 * 2000 + 4000+10000;
% % Round the sample to the nearest integer
% sample2 = round(sample2);
% % Clip the sample to the desired range of 1000 to 7000
% sample2(sample2 < 1000+10000) = 1000+10000;
% sample2(sample2 > 7000+10000) = 7000+10000;
%**************************************************************************

%sample = [sample1 sample2];
%sample = cat(2, sample1, sample2);
sample = sample1;

DigitalTimeofAdhesionStochastic = zeros(np,1);
for i = 1:np
    if r(i) < ProbabilityofUptaken
        DigitalTimeofAdhesionStochastic(i) = 0;
    else
        DigitalTimeofAdhesionStochastic(i) = sample(i); %DigitalTimeofAdhesion;
    end
end
%**************************************************************************

NumberTimeInteraction=round(TotalTime/ComputationSet.GTimeStep);  % Number of step time for diffusion calculations
Betta=zeros(1,Trials);



for IterTrial = 1:Trials % number of trials to compute the mean value of diffusion coefficient
    CreatMatrixAverage = 0.01;
    IndexMatrix = 1;
    Time=0:ComputationSet.GTimeStep:TotalTime;
    figure(1)
    % determination of initial matrices for balls
    [x,y,z] = sphere;
    r=zeros(Numbersphere,2);
    xd=zeros(Numbersphere,1);yd=zeros(Numbersphere,1);zd=zeros(Numbersphere,1);
    % determination of initial vector of diffusion coefficient
    DifCoefTotal=zeros(1);
    DifCoefTotalX=zeros(1);
    % seed number
    TotalNumberinteration=1;

    % rout of computaion in time
    for NumberIteretion = 1:  TotalNumberinteration
        TempTime=NumberIteretion*ComputationSet.GTimeStep;
        % setting the first ball so as not to occupy the center of the release of molecules
        r(1) =3; r(1,2) =0;   % small initial radius 1 is an astroglia
        xd(1)=1000;   % position x
        yd(1)=1000;
        zd(1)=1000;
        surf((r(1)*x+xd(1)), (r(1)*y+yd(1)),(r(1)*z+zd(1)))  % sphere centered at origin
        hold on

        Stak=ones(np, 1, 1);
        Stak(:,2)=10;
        Stak(:,3)=10;

        % second release
        %         for IterNP=1:round(np/2)
        %             Stak(IterNP,3)=10;
        %         end;
        x_par = zeros(1,np);
        y_par = zeros(1,np);
        z_par = zeros(1,np);
        TimeOfAdhesive = zeros(1, np);
        Distance=zeros(1);
        TempDistance=0;
        TempDistanceX=0;
        TempDistanceY=0;
        TempDistanceZ=0;

        DifCoef=zeros(1);
        NumberParticles=zeros(1);
        DifCoefX=zeros(1);
        DifCoefY=zeros(1);
        DifCoefZ=zeros(1);
        ParameterControl=0.5;
        ControlParam=0;
        SpaceBalls=(4/3)*3.14* r(1,1)^3;
        %    while  abs((100*(FSpace(IterTrial)-Balls.BettaControl)/Balls.BettaControl))>5
        for j=2:Numbersphere
            test=j;
            %             ************************************************************
            %              Overlapping ball distibution
            if strcmp(Balls.OverlapingBalls,'y')
                while 1
                    DistanceBetweenBalls = 2*ParameterControl+zeros(1,j,'uint32');
                    xd(j)=DistanceX*(rand()-0.5);
                    yd(j)=DistanceX*(rand()-0.5);
                    zd(j)=DistanceX*(rand()-0.5);
                    r(j,1) =(Balls.MinBallRadius + rand()*Balls.MaxBallRadius);
                    if rand() < Balls.ProbabilityofAstrocytes
                        r(j,2) =1; else
                        r(j,2) =0;
                    end
                    if sqrt(xd(j)^2+yd(j)^2+zd(j)^2) - Cleft.RadiusofCleft - r(j,1)>= Cleft.DistanceAstrocyte

                        %ParameterControl here the 0.1 is a minimum distance between balls
                        break
                    end
                end
            else
                %**************************************************************************
                %             Non overlapping balls
                %**************************************************************************
                while 1 % Conditions of non-overlapping balls
                    %
                    while 1
                        DistanceBetweenBalls = 2*ParameterControl+zeros(1,j,'uint32');
                        xd(j)=DistanceX*(rand()-0.5);
                        yd(j)=DistanceX*(rand()-0.5);
                        zd(j)=DistanceX*(rand()-0.5);
                        r(j,1) =(Balls.MinBallRadius + rand()*Balls.MaxBallRadius);
                        if rand() < Balls.ProbabilityofAstrocytes
                            r(j,2) =1; else
                            r(j,2) =0;
                        end
                        if sqrt(xd(j)^2+yd(j)^2+zd(j)^2) - Cleft.RadiusofCleft - r(j,1) >= Cleft.DistanceAstrocyte
                            %ParameterControl here the 0.1 is a minimum distance between balls
                            break
                        end
                    end

                    % Conditions of non-overlapping balls
                    for CheakIter=1:j-1
                        DistanceBetweenBalls(CheakIter)=sqrt((xd(j)-xd(CheakIter))^2 + (yd(j)-yd(CheakIter))^2 + (zd(j)-zd(CheakIter))^2)...
                            - (r(j)+r(CheakIter));

                    end
                    if min(DistanceBetweenBalls) >= 0.1
                        %ParameterControl here the 0.1 is a minimum distance between balls
                        break
                    end
                end % create the ball distribution
            end

            %             %*************************************************************
            %
            SpaceBalls=SpaceBalls+(4/3)*3.14* r(j,1)^3;
            %    Visualization of balls. For routine computation should be commented on.
            % However the visualization tooks some time and drasticaly decrease the time of computation.
            %**************************************************************************

            %             surf((r(j,1)*x+xd(j)), (r(j,1)*y+yd(j)),(r(j,1)*z+zd(j)),'FaceColor','red','EdgeColor','none')
            %             xlim([-SizeOfaxisX*AxisSize SizeOfaxisX*AxisSize])
            %             ylim([-AxisSize AxisSize])
            %             zlim([-AxisSize AxisSize])
            %             colormap([0 1 0])
            %             NumberSpherePrint=j;
            %**************************************************************************
        end % End line of filling the space by a balls
        hold on

        %******************************Write Balls file
        if strcmp(BasicParameters.SaveFiles,'y')
            Filename = sprintf('Balls distribution %s.txt', datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-ss'));
            writematrix([round(xd),  round(yd), round(zd), round(r)],Filename);
        end
        % *****************************


        % Modified at 04/10/2020

        FilledSpace=SpaceBalls/((DistanceX)^3 - (4/3)*3.14*Cleft.RadiusofCleft^3);


        % The calculation of betta, parameters of the  filled space by balls
        PositivePatricle=ones(1,NumberControlParticel);
        xTest=linspace(0,0,NumberControlParticel);
        yTest=linspace(0,0,NumberControlParticel);
        zTest=linspace(0,0,NumberControlParticel);

        for ic=1:NumberControlParticel
            while 1
                xTest(ic)=DistanceX*(rand()-0.5);
                yTest(ic)=DistanceX*(rand()-0.5);
                zTest(ic)=DistanceX*(rand()-0.5);
                if sqrt(xTest(ic)^2+yTest(ic)^2+zTest(ic)^2) - Cleft.RadiusofCleft > Cleft.DistanceAstrocyte
                    %ParameterControl here the 0.1 is a minimum distance between balls
                    %TempDistance(i) = sqrt(xTest(ic)^2+yTest(ic)^2+zTest(ic)^2);
                    break
                end
            end

            for j=2:Numbersphere
                VectorControl = sqrt((xTest(ic)-xd(j))^2 + (yTest(ic)-yd(j))^2 + (zTest(ic)-zd(j))^2) - r(j,1);
                if   VectorControl <0
                    PositivePatricle(ic)=-1;
                end

            end
        end
        % Calculation Betta and Alpha
        InPar=sum(PositivePatricle(:) < 0);

        %   Betta
        Betta(IterTrial)=InPar/NumberControlParticel;
        %************************************************************
        %   Alpha
        AlphaSpace = 1 - Betta(IterTrial);
        %************************************************************

        % The line of end of computation of filled space parameter Betta

        % Test1 = (100*abs(Betta(IterTrial)-Balls.BettaControl)/Balls.BettaControl);
        %      end % while
        %******************************************************************
        daspect([1 1 1])
        % Iteration to store positions of particles

        p_plot=scatter3(x_par(1,:),y_par(1,:),z_par(1,:),20,'filled');
        axis([-SizeOfaxisX*AxisSize SizeOfaxisX*AxisSize  -AxisSize AxisSize -AxisSize AxisSize]);
        TimeIter=0;
        InputPar = [np; Numbersphere];
        AverageTimeMatrix = zeros(DistanceBins,TimeBins+1);
        AverageTimeMatrixFree = zeros(DistanceBins,TimeBins+1);
        %         Computation in Time
% === Snapshot settings (fractions of TotalTime) ===
SnapshotTimes = [0.0125, 0.0375, 0.125, 0.375];   % ≈ 0.1, 0.3, 1, 3 ms if TotalTime=8 ms
t_targets_ms  = SnapshotTimes * TotalTime;         % target times in ms
saved         = false(size(t_targets_ms));         % one-time save flags

dt_ms  = ComputationSet.GTimeStep;                 % time step in ms
t_prev = 0;                                        % previous time (ms)

        for i=1:NumberTimeInteraction
          
            % 10 is a time of second release
            if i == round(10/ComputationSet.GTimeStep)
                Stak(:,3)=8;
            end
            TimeIter=TimeIter+1;
            %     Visualization of calculation in time
            if strcmp(BasicParameters.SaveFiles,'y')
                TimeOutPercent = 1/TimeBins;
                if TimeIter == round(TimeOutPercent*TotalTime/ComputationSet.GTimeStep)
                    TimeofDiffusion = i * ComputationSet.GTimeStep;
                    formatSpec = 'Time of diffusion  %8.4f  ms and total time is  %8.4f ms\n';
                    fprintf(formatSpec,TimeofDiffusion, TotalTime);
                    TimeIter=0;
                    if (i == round(CreatMatrixAverage*TotalTime/ComputationSet.GTimeStep))
                        SQBound=zeros(1,BasicParameters.NumberOfParticles);
                        SQFree=zeros(1,BasicParameters.NumberOfParticles);
                        TempNew = (x_par(i,:).^2 + y_par(i,:).^2+ z_par(i,:).^2 ).^0.5;
                        TempNew(TempNew < Cleft.RadiusofCleft)=0;
                        ttt= Stak(:,2);
                        SQBound(Stak(:,1)>1)=TempNew(Stak(:,1)>1);
                        SQFree(Stak(:,1)<1.1)=TempNew(Stak(:,1)<1.1);
                        DistanceBinsNew=(Cleft.RadiusofCleft+round(DistanceX/DistanceBins)):round(DistanceX/DistanceBins):(DistanceX+Cleft.RadiusofCleft+round(DistanceX/DistanceBins));

                        TimeBinsNew = 0:1:TimeBins;
                        [NUmm,HSQ] = histcounts(SQBound,DistanceBinsNew);
                        HSQ = HSQ(2:end);
                        [NUmmFree,HSQFree] = histcounts(SQFree,DistanceBinsNew);
                        HSQFree = HSQFree(2:end);
                        %PrintV=[HSQ; NUmm];
                        if CreatMatrixAverage < TimeOutPercent+0.001
                            %(rows, colms)
                            AverageTimeMatrix(1:end,1) = HSQ/1000;
                            AverageTimeMatrix(1:end,2) = NUmm;
                            AverageTimeMatrixFree(1:end,1) = HSQFree/1000;
                            AverageTimeMatrixFree(1:end,2) = NUmmFree;
                        else
                            IndexMatrix=IndexMatrix+1;
                            AverageTimeMatrix(1:end,IndexMatrix+1) = NUmm;
                            AverageTimeMatrixFree(1:end,IndexMatrix+1) = NUmmFree;
                        end
                        CreatMatrixAverage = CreatMatrixAverage+TimeOutPercent;
                    end
                                                        % previous time (start at 0 ms)


              

                end
            end % visualization

            TempNumber=0;


            for j=1:np % iteraction np

                if Stak(j,3) > 0
                    TempPosition = sqrt(x_par(i,j)^2+y_par(i,j)^2+z_par(i,j)^2);
    				%****************************************************************************************************************************
                    if TempPosition > Cleft.RadiusofCleft %
                        %  if (sqrt((x_par(i,j))^2)  > 250^2) || (sqrt((y_par(i,j))^2)  > 50^2) || (sqrt((z_par(i,j))^2)  > 50^2)
                        if Stak( j, 1) == 1
                            x_par(i+1,j)=x_par(i,j)+ ComputationSet.StepRad*2*(rand()-0.5);
                            y_par(i+1,j)=y_par(i,j)+ ComputationSet.StepRad*2*(rand()-0.5);
                            z_par(i+1,j)=z_par(i,j)+ ComputationSet.StepRad*2*(rand()-0.5);

                        else
                            x_par(i+1,j)=x_par(i,j);
                            y_par(i+1,j)=y_par(i,j);
                            z_par(i+1,j)=z_par(i,j);
                        end
                        for j_Sphere=1:Numbersphere
                            % this procedure reduces the number of sticky balls with time scale timescale =
                            %                          TimeScale = 1000;
                            %                          if r(j_Sphere,2)==1
                            %                              if (1-exp(-(ComputationSet.GTimeStep/TimeScale))) > rand()
                            %                                  r(j_Sphere,2) = 0;
                            %                              end
                            %                          end
                            % Next condition determines the location of a particle that should have entered the ball
                            % but instead remains outside of it. The particle stays in place until it randomly moves away from the ball.
                            % Reflecting a particle from a ball does not work well because the reflection can lead a particle into another ball,
                            % requiring the programme to check the particle's new position each time, resulting in a very large number of checks (infinite calculations)
                            % in the case of a high density of balls.

                            if sqrt((x_par(i+1,j)-xd(j_Sphere))^2+(y_par(i+1,j)-yd(j_Sphere))^2+(z_par(i+1,j)-zd(j_Sphere))^2) < r(j_Sphere,1)
                                if TempPosition > NonAdhesiveDistance
                                    if r(j_Sphere,2)==1
                                        %************************************************
                                        TimeOfAdhesive(j) = TimeOfAdhesive(j) + 1;
                                        if Adhesion.MaxProbAdhesive*(1-exp(-(ComputationSet.GTimeStep*TimeOfAdhesive(j))/Adhesion.TimeINsideAdhesiveZone)) > rand()
                                            Stak(j, 1)=Stak(j, 1)+1;
                                            TimeOfAdhesive(j)=0;
                                        end
                                    end
                                end
                                % ComputationSet.StepLocalMovement is a step of
                                % small walking near the ball
                                x_par(i+1,j)=x_par(i,j)+ ComputationSet.StepLocalMovement*ComputationSet.StepRad*2*(rand()-0.5);
                                y_par(i+1,j)=y_par(i,j)+ ComputationSet.StepLocalMovement*ComputationSet.StepRad*2*(rand()-0.5);
                                z_par(i+1,j)=z_par(i,j)+ ComputationSet.StepLocalMovement*ComputationSet.StepRad*2*(rand()-0.5);
                            end
                        end % end of sphere control

                        if Stak(j, 1)  >1
                            Stak(j, 1)=Stak(j, 1)+1;
                        end
                        if Stak(j, 1)  == DigitalTimeofAdhesionStochastic(j) %DigitalTimeofAdhesion
                            Stak(j, 1)=1;
                        end
                        % conditions an open end, which allows to exclude particles
                        % from consideration of the diffusion coefficient if they fly outside the boundary

                        if (sqrt((x_par(i+1,j))^2)  < DistanceX/2) || (sqrt((y_par(i+1,j))^2)  < DistanceX/2) || (sqrt((z_par(i+1,j))^2)  < DistanceX/2)
                            TempDistance = TempDistance +  (x_par(i+1,j)^2+y_par(i+1,j)^2+z_par(i+1,j)^2)/(6*i*ComputationSet.GTimeStep);
                            TempDistanceX= TempDistanceX + (x_par(i+1,j)^2)/(2*i*ComputationSet.GTimeStep);
                            TempNumber=TempNumber+1;
                        end
                    else
                        if (sqrt((x_par(i,j))^2) < Cleft.RadiusofCleft)
                            x_par(i+1,j)=x_par(i,j)+ ComputationSet.StepRad*2*(rand()-0.5);
                        else
                            x_par(i+1,j)=x_par(i,j);
                        end
                        if  (sqrt((y_par(i,j))^2)  < Cleft.RadiusofCleft)
                            y_par(i+1,j)=y_par(i,j)+ ComputationSet.StepRad*2*(rand()-0.5);
                        else
                            y_par(i+1,j)=y_par(i,j) ;
                        end
                        if (sqrt((z_par(i,j))^2)  < Cleft.CleftWidth)
                            z_par(i+1,j)=z_par(i,j)+ 0.1*ComputationSet.StepRad*2*(rand()-0.5);
                        else
                            if z_par(i,j) < 0
                                z_par(i+1,j)= z_par(i,j) + 0.1*Cleft.CleftWidth;
                            end
                            if z_par(i,j) > 0
                                z_par(i+1,j)= z_par(i,j) - 0.1*Cleft.CleftWidth;
                            end
                        end
                        if (sqrt((x_par(i+1,j))^2)  < DistanceX/2) || (sqrt((y_par(i+1,j))^2)  < DistanceX/2) || (sqrt((z_par(i+1,j))^2)  < DistanceX/2)
                            TempDistance = TempDistance +  (x_par(i+1,j)^2+y_par(i+1,j)^2+z_par(i+1,j)^2)/(6*i*ComputationSet.GTimeStep);
                            TempDistanceX= TempDistanceX + (x_par(i+1,j)^2)/(2*i*ComputationSet.GTimeStep);
                            TempNumber=TempNumber+1;
                        end
                    end
                    % condition that particles do not fly out of the pipe

                    if (sqrt(x_par(i,j)^2+y_par(i,j)^2+z_par(i,j)^2)) > SizeCube.Size/2
                        %                         x_par(i+1,j)=x_par(i-5,j);
                        %                         y_par(i+1,j)=y_par(i-5,j);
                        %                         z_par(i+1,j)=z_par(i-5,j);


                    end
                else
                    x_par(i+1,j)=x_par(i,j);
                    y_par(i+1,j)=y_par(i,j);
                    z_par(i+1,j)=z_par(i,j);

                end
                if (sqrt(x_par(i,j)^2+y_par(i,j)^2+z_par(i,j)^2)) > 4500
                    x_par(i+1,j)=x_par(i,j)+5000;
                    y_par(i+1,j)=y_par(i,j)+5000;
                    z_par(i+1,j)=z_par(i,j)+5000;


                end
            end % iteraction in np
% === Snapshot capture system (time-crossing detection) ===
t_curr = i * dt_ms;

to_save = find(~saved & (t_prev < t_targets_ms) & (t_targets_ms <= t_curr));

for k = 1:numel(to_save)
    ss   = to_save(k);
    t_ms = t_targets_ms(ss);  % label snapshots by target time

    % bound -> 0, free -> 1
    StakNew = Stak(:,1);
    StakNew(StakNew > 1) = 0;
    A = [x_par(i,:); y_par(i,:); z_par(i,:); StakNew'];

    % Safe filename (no risky characters)
    Filename = sprintf('PD_%.4fms_%s.txt', t_ms, ...
        datetime('now','TimeZone','local','Format','yyyy-MM-dd_HH-mm-ss'));

    fileID = fopen(Filename, 'w');
    if fileID < 0
        error('File write failed: %s', Filename);
    end
    fprintf(fileID,'%10.5f %10.5f %10.5f %10.0f \n', A);
    fclose(fileID);

    saved(ss) = true;  % one-time save
end

t_prev = t_curr;
% === End snapshot capture ===


            set(p_plot,'XData',x_par(i,:));
            set(p_plot,'YData',y_par(i,:));
            set(p_plot,'ZData',z_par(i,:));
            NumberParticles(i+1)=TempNumber;
            TempNumberParticla=TempNumber;
            DifCoef(i+1)=(TempDistance/TempNumber)/10^6;
            DifCoefX(i+1)=(TempDistanceX/TempNumber)/10^6;
            TempNumber=0;
            TempDistance=0;
            TempDistanceX=0;
            pause(0.0001);

            if (TempNumberParticla < np) || (i == NumberTimeInteraction - 1)
                if  (BasicParameters.Control == 1)
                    A = [x_par(i,:);  y_par(i,:); z_par(i,:)];
                    InputParameters = [charge; Numbersphere];
                    fileID = fopen('FinalDistribution.txt','w');
                    fprintf(fileID,'%10.5f %10.5f\n', InputParameters);
                    fprintf(fileID,'%10.5f %10.5f %10.5f\n', A);
                    Control = 2;
                    fclose(fileID);
                end
            end



        end % in time

        DifCoefTotal=DifCoefTotal + DifCoef;
        DifCoefTotalX=DifCoefTotalX + DifCoefX;



    end % Trial Number


    %Filename = sprintf('DistanceBound %s.txt', datestr(now,' mmmm dd, yyyy HH-MM-SS AM'));
    Filename = sprintf('DistanceBound %s.txt',datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-ss'));
    fileID = fopen(Filename,'w');
    %writematrix(AverageTimeMatrix,Filename,'Delimiter','tab')
    writematrix(AverageTimeMatrix,Filename)
    fclose(fileID);

    Filename = sprintf('DistanceFree %s.txt',datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-ss'));
    fileID = fopen(Filename,'w');
    writematrix(AverageTimeMatrixFree,Filename)
    fclose(fileID);
end
% computation of diffusion coefficient
Begin= fix(BeginTime/ComputationSet.GTimeStep);
End=fix(EndTime/ComputationSet.GTimeStep);
MeanDiffusion = mean( DifCoefTotalX(Begin : End));
STDDiffusion = std( DifCoefTotalX(Begin : End));
% end of computation
figure(2)
%     plot(Time, Distance)
%     plot(Time, DifCoefTotalX/TotalNumberinteration)
%     cmap = hsv(np); % Creates a np-by-3 set of colors from the HSV colormap

DiffResults = [MeanDiffusion; STDDiffusion];
fprintf(fileIDStat,'%10.5f %10.5f %10.5f\n', DiffResults, Betta(IterTrial));

fclose(fileIDStat);
A = [Time;  DifCoefTotalX/TotalNumberinteration; NumberParticles];
InputParameters = [charge; Numbersphere];
Filename = sprintf('Diffusion_Set_ms_%s.txt', datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-ss'));
fileID = fopen(Filename,'w');
fprintf(fileID,'%10.5f %10.5f\n', InputParameters);
fprintf(fileID,'%10.5f %10.5f %10.5f\n', A);
fclose(fileID);










