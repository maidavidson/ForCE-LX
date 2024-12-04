%% Setup file for ForCE_LS
% Test for longshore transport and structures of variable length - constant
% waves

%% Force model version name (controls which programe version runs)
 model.ver='ForCE_LS';

%% Input/output datafile names (base directory is the project name directory, in which two folders must exit: ontputData & inputData). Filenames can be user defined.
 filename.waveData='./inputData/ConstWaves'; % wave data path/file name
 filename.seaLevel='./inputData/seaLevelData0'; % pAnalyticalath/file containing sea level change data
 filename.shoreline='./inputData/StraightShoreline'; % Shoreline data file
 filename.output='./outputData/TestConstantWaves.mat'; %
 filename.structures='./inputData/SynStructureBWs2'; % Structures data file

%% Model calibration parameters
 model.k1=0.00;  % Dimensionless cross-shore sediment transport coefficent 
 model.k2=0.01;  % Longshore transport coefficient 
 model.a1=0;     % Shoreline feedback cross-shore
 model.a2=0;     % Shoreline feedback long-shore
 model.b1=1;     % Constant term - cross-shore
 model.b2=0;     % Constant term - longshore

%% Choose numerical scheme
 % F=0 (explicit, accuate but low stability), F=1 (implicit, stable but
 % less acurate), F=0.5 (default) (best stability and accuacy)
 model.F=0.5;
 
%% Model temporal and spatial resolution
 model.dt=1; % [days] Define model time-step
 model.ds=20; % [m] Model spatial resolution
 
%% User defined dates (check data limits!)
 model.startDate=datenum(0,1,1,0,0,0); % Leave undefined =[] to start with first survey date
 model.nYears=5; % [years] Number of simulation years
 
%% Model site details
 model.site='Analytic Test';
 model.sedimentDensity=2600; % [kg/m^3]
 model.waterDensity=1025; % [kg/m^3]
 model.poreSpace=0.3; % [Dimensionless]
 model.MSR=0;
 model.gamma=0.78; % Breaker coefficient
 model.D50=35;% mm
 model.A=0.21*model.D50.^0.48;
 
%% Threshold distance for transmission coefficient calcs
 model.thresholdH=0.; % Wave height threshold for sediment motion (default =0m)
 model.shadowSwitch=0; % 1=shadowing on, 0=offsholdTheta=inf; % Shoreline evolution angle before re-digitistion of shoreline

%% Define constants
 constant.days2seconds=24*60*60; % days to seconds conversion factor
 constant.g=9.81; % [m/s^2] gravity
 constant.k1=1/((model.sedimentDensity-model.waterDensity)*model.poreSpace*constant.g); % Sediment transport coefficent 

%% Model variables 
 model.nPointsYear=round(365.25/model.dt);  % Number of points in a year
 model.Npts=model.nYears*model.nPointsYear;

%% Plotting and saving
 model.plotStep=7; % Model plot-step [days]
 model.saveStep=7; % Model save-step [days]

%% Transmission around obstacles
model.kTr=1; % Multiplyer for the effective length (Leffective) of the obstacle
model.TrExp=1; % If definned this uses: TrCoeff=exp(-(Leffective/Xsurf).^0.5) else TrCoeff = 1 - (Leffective/Xsurf)

