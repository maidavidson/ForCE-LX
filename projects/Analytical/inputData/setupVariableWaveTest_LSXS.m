%% Setup file for ForCE_LS
%% This is a variabe wave test including longshore and cross-shore transport

%% Force model version name (controls which programe version runs)
 model.ver='ForCE_LS';

%% Input/output data file names (base directory is the project name directory, in which two folders must exit: ontputData & inputData). Filenames can be user defined.
 filename.waveData='./inputData/VariableWaves';  % wave data path/file name
 filename.seaLevel='./inputData/seaLevelData0';     % pAnalyticalath/file containing sea level change data
 filename.shoreline='./inputData/StraightShoreline';% Shoreline data file
 filename.output='./outputData/TestVariableWaves.mat';       % Output Filename
 filename.structures='./inputData/SynStructureBWs2';% Structures data file

%% Model calibration parameters
 model.k1=0.01;  % Dimensionless cross-shore sediment transport coefficent 
 model.k2=0.01;  % Longshore transport coefficient 
 model.a1=-10;   % Shoreline feedback cross-shore
 model.a2=1e-30; % Shoreline feedback long-shore
 model.b1=1;     % Constant term - cross-shore
 model.b2=0;     % Constant term - longshore

%% Choose numerical scheme
 % F=0 (explicit, accuate but low stability), F=1 (implicit, stable but
 % less acurate), F=0.5 (default) (best stability and accuacy)
 model.F=.1;
 
%% Model temporal and spatial resolution
 model.dt=1; % [days] Define model time-step
 model.ds=20; % [m] Model spatial resolution
 
%% User defined dates (check data limits!)
 model.startDate=datenum(0,1,1,0,0,0); % Leave undefined =[] to start with first survey date
 model.nYears=9; % [years] Number of simulation years
 
%% Model site details
 model.site='Analytic Test';
 model.sedimentDensity=2600; % [kg/m^3]
 model.waterDensity=1025; % [kg/m^3]
 model.poreSpace=0.3; % [Dimensionless]
 model.MSR=3;
 model.gamma=0.78; % Breaker coefficient
 model.D50=3.5;% mm
 model.A=0.21*model.D50.^0.48; 

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