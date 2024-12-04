%% Setup file for ForCE_LS

%% Define filenames for input data
% Note that...
% 1) The wave data file should contain a structure called 'wave' with the
% following fields relevant to the offshore boundary depth of the model... 
% wave.dates - matlab datnum [days]
% wave.H, - significant wave height
% wave.Tp - peak wave period
% wave.Dir - wave directions rel. to shore-normal
% wave.dt - wave samplint interval
% 2) the sea level file should contain a structure called 'seaLevel' with
% the following fields...
% seaLevel.dates - matlab datnum [days]
% seaLevel.z - water level [m] rel. MSL
% seaLevel.dt - sea level sampling interval [days]


%% Force model version name (controls which programe version runs)
 model.ver='ForCE_LS';

%% Input/output datafile names (base directory is the project name directory, in which two folders must exit: ontputData & inputData). Filenames can be user defined.
 filename.waveData='./inputData/synWavesVariableNarra'; % wave data path/file name
 filename.tide='./inputData/tideDataWithSurge'; % path/file containing tide data
 filename.seaLevel='./inputData/seaLevelData0'; % pAnalyticalath/file containing sea level change data
 filename.shoreline='./inputData/synNarraShoreline'; % Shoreline data file
 filename.output='./outputData/TestLXa1_0_b2_0'; %
 %filename.structures='./inputData/opFileStructuresAnalytic'; % Structures data file

%% Model calibration parameters
 model.k1=0.01;  % Dimensionless cross-shore sediment transport coefficent 
 model.k2=0.01;  % Longshore transport coefficient 
 model.a1=0;   % Shoreline feedback cross-shore
 model.a2=0;     % Shoreline feedback long-shore
 model.b1=1;     % Constant term - cross-shore
 model.b2=0;     % Constant term - longshore

%% Choose numerical scheme
 % F=1 (explicit, accuate but low stability), F=0 (implicit, stable but
 % less acurate), F=0.5 (default) (best stability and accuacy)
 model.F=1;
 
%% Model temporal and spatial resolution
 model.dt=1; % [days] Define model time-step
 model.ds=20; % [m] Model spatial resolution
 
%% User defined dates (check data limits!)
 model.startDate=datenum(0,1,1,0,0,0); % Leave undefined =[] to start with first survey date
 model.nYears=30; % [years] Number of simulation years
 
%% Model site details
 model.site='Analytic Test';
 model.sedimentDensity=2600; % [kg/m^3]
 model.waterDensity=1025; % [kg/m^3]
 model.poreSpace=0.6; % [Dimensionless]
 model.MSR=3;
 model.gamma=0.78; % Breaker coefficient
 model.D50=.35;% mm
 model.A=0.21*model.D50.^0.48;

 % Define model boundary conditions at the edges of the domain (optional).
 % NB can also be achieved with structures that shadow this region. Here
 % the characteristic length scales are defined. If unspecified or zero,
 % these are open boundaries. Seting to inf is a zero flux boundary.
 model.L(1)=inf;
 model.L(2)=inf;
 
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
 model.axisEqual=1; % Define to make x-y axis scaling the same