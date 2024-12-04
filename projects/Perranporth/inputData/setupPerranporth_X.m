%% Setup file for ForCE_LX

%% Force model version name (controls which programe version runs)
 model.ver='ForCE_LX';

%% Input/output datafile names (base directory is the project name directory, in which two folders must exit: ontputData & inputData). Filenames can be user defined.
 filename.waveData='./inputData/perranWavesMetOffice';  % wave data path/file name
 filename.seaLevel='./inputData/seaLevelData0';         % path/file containing sea level change data
 filename.shoreline='./inputData/perranShoreline';      % Shoreline data file
 filename.structures='./inputData/perranStructures';    % Shoreline data file
 filename.output='./outputData/perranTest_X';             % Output filename

%% Model calibration parameters
 model.k1=0.04;  % Dimensionless cross-shore sediment transport coefficent 
 model.k2=0.0;  % Longshore transport coefficient 
 model.a1=-60;   % Shoreline feedback cross-shore
 model.a2=0;     % Shoreline feedback long-shore
 model.b1=1;     % Constant term - cross-shore
 model.b2=1;     % Constant term - longshore

 %% Choose numerical scheme
 % F=0 (explicit, accuate but low stability), F=1 (implicit, stable but
 % less acurate), F=0.5 (default) (best stability and accuacy)
 model.F=0.5;

%% Model temporal and spatial resolution
 model.dt=1; % [days] Define model time-step
 model.ds=20; % [m] Model spatial resolution
 
%% User defined dates (check data limits!)
 model.startDate=datenum(2007,1,1,0,0,0); % Must be > wave data start date
 model.nYears=30; % [years] Number of simulation years
 
%% Model site details
 model.site='Perranporth, UK';
 model.sedimentDensity=2600; % [kg/m^3]
 model.waterDensity=1025; % [kg/m^3]
 model.poreSpace=0.6; % [Dimensionless]
 model.MSR=6;
 model.gamma=0.78; % Breaker coefficient
 model.D50=0.35;% mm
 model.A=0.21*model.D50.^0.48;
 model.kTr=1; % Effective length booster for obstructions.


%% Define constants
 constant.days2seconds=24*60*60; % days to seconds conversion factor
 constant.g=9.81; % [m/s^2] gravity
 constant.k1=1/((model.sedimentDensity-model.waterDensity)*model.poreSpace*constant.g); % Sediment transport coefficent 

 %% Plotting and saving
 model.plotStep=30; % Model plot-step [days]
 model.saveStep=7;  % Model save-step [days]
 model.axisEqual=1; % Makes map plot axis equal

