%% ForCE_LX program 
% Mark Davidson 13/6/2024 


%% Load input files [see setup file for details of the required file contents]
    load(filename.waveData);         % Load wave data 
    load(filename.shoreline);        % Load initial shoreline
    load(filename.seaLevel);         % Load SLR data 

%% If there are structures, then load
    if isfield(filename,'structures')
        load(filename.structures);         % Load structures
        model.structures=1;
        for i =1:structures(1).Nstruct
            structures(i).xo=mean(structures(i).x);
            structures(i).yo=mean(structures(i).y);
        end
    end


%% (semi) - Impicit or explicit schemes
    if isfield(model,'F')
        F=model.F;
    else
        F=0.1; % Default Crank-Nicholson, Lax-Wendroff scheme
    end

%% Initial sea level
    SLR0=interp1(seaLevel.dates,seaLevel.z,model.startDate,'linear','extrap');
    
%% Define initial shoreline x and y values
% re-define model.ds to fit required shoreline length
     model.Ns=length(shoreline.x);
    [model]=find_segments(shoreline.Xhead,0,model);

%% Seperate mobile beach from non-mobile coast. 
    [model]=getBeach(shoreline,model);
    [model,hard]=getNonMobile(shoreline,model);

%% Make hard coastline string to compute shadows...
if hard(1).Nseg>1
        model.x_hard=[];
        model.y_hard=[];
        if model.shadowSwitch==1
        for i=1:hard(1).Nseg
            hard(i).Ns=length(hard(i).x);
            [x,y,hard(i).ds,hard(i).Ns]=...
                initialiseShoreline(hard(i).x,hard(i).y,hard(i).Ns,10);
            model.x_hard=[model.x_hard,x];
            model.x_hard=[model.x_hard,NaN];
            model.y_hard=[model.y_hard,y];
            model.y_hard=[model.y_hard,NaN];
            clear x y
        end
        end
end

% Resample shoreline acording to model.ds
    [model.x,model.y,model.dsMod,model.Ns]=initialiseShoreline(model.x,model.y,model.Ns,model.ds);
    model.x0=model.x;model.y0=model.y;
    model.dn=zeros(size(model.x)); % Initialise shoreline normal displacement
    model.nx=model.dn;
    model.ns=model.dn;
    model.dnSLR=model.dn;

% See how many wave gauges    
    wave(1).N=length(wave);

% Find nearest gauge to segment (UTM WGS84)
    if wave(1).N>2
        for ig=1:wave(1).N
            distToGauge(ig)=hypot(mean(model.x0)-wave(ig).x(1), mean(model.y)-wave(ig).y(1));
        end
        inearest=find(distToGauge==min(distToGauge));
    else
        inearest=1;
    end

% Record wave buoy depth for this segment
    model.waveDepth=wave(inearest).depth;
    clear distToGauge

% Find shoreline normal direction
    [model.theta]=findNormal(smooth(model.x0,30),smooth(model.y0,30));
    model.theta=model.theta(:)';
    model.theta0=model.theta;

% Compute weights for wave nodes    
    [model.weights,wave]=computeWeights(smooth(model.x0,30),smooth(model.y0,30),wave);

% Compute the average depth of closure surfzone widths and slopes for all locations
    [model]=computeAverageValues(wave,model);

% Compute average power components
    [model]=computeAveragePowerComponents(wave,model,wave(1).dates);
    DirArray=model.DirMean;
    Paverage=mean(model.PoMean);
    QoAv=constant.k1*Paverage;
    model.Pavn=model.PnMean;
    model.Pavs=model.PsMean;
    model.Pav=model.PoMean;
    Pav1=model.Pavn;
    Pav2=model.Pavs;
    Berm=model.Berm;
    % NB the following has been modified to give local values (15/04/2024)
    Xsurf=model.Xsurf;
    dc=model.dc;
    tanb=model.tanb;
    Dbar=mean(model.Pav./Xsurf);

% Initialise fluxes    
    model.qs=zeros(size(model.x));
    model.qlhb=0;
    model.qrhb=0;
    model.qlhbav=0;
    model.qrhbav=0;
    % Kalman filter shoreline perturbations
    dnx_corr_all=zeros(size(model.x));
    dns_corr_all=zeros(size(model.x));

    
% If there are calibrations, then load file
    Nval=length(model.startDate:model.startDate+model.nYears*365);
    if isfield(filename,'calibration')
        nState=length(KF(1).H);
        load(filename.calibration);         % Load claibration file if specified
        [nSurveys,nProfiles]=size(calibration.MCLArray);
        %R=KF.Noise.^2;
        KF(1).used=1;
        disp('Calibration file loaded')
        KF(1).R=(KF(1).Noise*ones(1,calibration.nProfiles)).^2; % Noise Matrix
        State=zeros(nState,Nval);
        State(:,1)=[0, 0, model.k1, model.k2, model.a1 model.a2 model.b1 model.b2]'    ;   % State Matrix
        P=KF(1).P;
        KF(1).timeMeasurements=calibration.allDates;
        tm=find(calibration.allDates>=model.startDate,1,"first");
        KF(1).trigger=calibration.allDates(tm);
       

 % Define J and P for all locations
        colourFix=[rand(nProfiles,1),rand(nProfiles,1),rand(nProfiles,1)];
        for j=1:nProfiles
            KF(j).Flag=0;
            KF(j).P=P;
            KF(j).J=zeros(nState,nState);
            if isfield(model,'syn')
                x1(j)=calibration.x0(j);
                y1(j)=calibration.y0(j);
            else % Field UTM profile positions
                x1(j)=calibration.UTMx0(j);
                y1(j)=calibration.UTMy0(j);
            end
            distToProf=hypot(model.x0-x1(j), model.y0-y1(j));
            jj=find(distToProf==min(distToProf));
            KF(j).index=jj;
            KF(j).distToProf=distToProf(jj);
            clear jj distToProf 
            critDist=hypot(60,model.ds);
            if KF(j).distToProf>critDist
              KF(j).Flag=0
            else;  
              KF(j).Flag=1; % Flag to indicate if profile is good
            end
        end

% initialise parameters for state param plot - below initialise model date
        nStatePlot=nState-2;
        evolve_upLim = zeros(nStatePlot, 1);
        evolve_lowLim = zeros(nStatePlot, 1);
        evolve_date=zeros(1,1);
        stateName={'XS','LS','k_1','k_2','a_1','a_2','b_1','b_2'};  

        for j=1:nSurveys
            jj=find(~isnan(calibration.MCLArray(j,:)));
            nCal(j)=length(jj);
            profileIndex{j}=jj;
        end

    else
        KF.used=0;
    end  

% Compute shadows lookup table for speed % Not recommended
    if isfield(model,'shadowSwitch')
        shadowSwitch=model.shadowSwitch;
    else
        shadowSwitch=0;
    end

% Initialise model date 
   date(1)=model.startDate;
   date(2)=model.startDate;
   startDate=model.startDate;
   nYears=model.nYears;
   dt=model.dt*constant.days2seconds;
   model.date(1)=startDate;
   nextPlot=startDate; % plot timing
   nextSave=startDate; % File save interval
   iSave=0; model.xArray=[];

%% Main program loop starts here [time-loop] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k=1;

while date(k)<=startDate+(nYears*365.25)

% Increment time counter
    k=k+1;

% Increment the date
    date(k)=date(k-1)+dt/constant.days2seconds;
   
% Break if out of wave data 
    if date(k)>wave(1).dates(end); break ;end
   
% Update shore-normal vales
    [model.theta]=findNormal(model.x,model.y);

% Compute the average depth of closure surfzone widths and slopes for all locations
    [model]=computeAverageValues(wave,model);

% Get sea level rise
    SLR=interp1(seaLevel.dates,seaLevel.z,date(k),'linear','extrap');

% Compute average power components
    [model]=computeAveragePowerComponents(wave,model,date(k));
    DirArray=model.DirMean;
    Paverage=mean(model.PoMean);
    QoAv=constant.k1*Paverage;
    Pn=model.PnMean;
    Ps=model.PsMean;
    Po=model.PoMean;
    XsurfNow=model.Xsurf; % gets the instantaneous surfzone width at the current time-step (used for tansmission calcs)
   % dc=model.dc; % Commented to keep constant in time

% Compute shadows on the basis of the non-shoaled wave direction
    % Initialise arays
       model.shadow=ones(size(model.x));
       shadow=ones(size(model.x));

      % Shadows from hard segments (replaced by structures now) - May be delete? 
    if shadowSwitch==1 && hard(1).Nseg>1
    %[model.shadow]=lookupShadow(model,Dir,i);
    [ xS,yS,shadow] = ...
        find_shadows(model.x,model.y,model.x_hard,model.y_hard,deg2rad(mode(model.DirAbs)),1 );
    model.shadow=interp1(1:length(shadow)-0,double(~shadow),1:model.Ns,'nearest','extrap');
    end

    
    % Now do the structures
    if isfield(model,'structures')
        for in= 1: structures(1).Nstruct
        % Find nearest shoreline point to structure
        distVal=hypot(structures(in).xo-model.x,structures(in).yo-model.y);
        iVal=find(distVal==min(distVal));
        [ xS,yS,shadow2] = ...
        find_shadows(model.x,model.y,structures(in).x(:)',structures(in).y(:)',deg2rad(model.DirAbs(iVal)),1 );
        offset=1;
        shadowStruct=interp1(offset:length(shadow2)+0,double(~shadow2),1:model.Ns,'nearest','extrap');
        jn=find(shadowStruct==0);
        if ~isempty(structures(in).L)
        [structures(in).transCoeff,Leffective1]=calcTransmission(structures(in).L,XsurfNow(iVal(1)),model.dn(iVal(1)),model);
        end
        shadowStruct(jn)=structures(in).transCoeff;
        %shadowStruct(jn)= ones(length(jn),1)-hanning(length(jn)).*structures(in).transCoeff;
        model.shadow=model.shadow.*shadowStruct; clear shadowStruct
        end
    end

% Find waves approaching at more than 90 degrees. NB High angle waves upto 90 degrees are ok.
    j=find(abs(DirArray)>=90);Ps(j)=0;Po(j)=0;  

% longshore flux
    ds=model.dsMod;
    
    if i==1;dns_ds=zeros(size(model.ns));end

    model.Pavs=model.a2.*Dbar.*model.ns + model.b2.*Pav2;

    qs = constant.k1.*model.k2.* (Ps-model.Pavs);
    Qo = 0.5*constant.k1.*model.k2.* model.PoMean(:)';

% Compute dissequilibrium
% The following is a stream-function between ShoreFor and Yates to be
% optimsed for a nd b later
     model.Pavn=model.a1.*Dbar.*model.nx + model.b1.*Pav1;
     Diss=(model.Pavn-Pn);
     
% Compute cross-shore sediment flux
    qn = constant.k1.*model.k1.*Diss;

% Apply shadows
    qs=qs.*model.shadow;

% Nullify high angle fluxes
    qs(j)=0;qn(j)=0;clear j
    model.qs=qs; % Record longshore flux after application of shadows

% Source Fuction definition (Accounts for coastine curvature in implicit scheme)    
    S=zeros(size(qs));
    D=Qo*dt./(dc.*ds.^2); % Diffusson coeff.
        for is=2:model.Ns-1 
        % NB Perhaps do a check on delta(DirArray) here?
        deltaTheta=DirArray(is+1) - DirArray(is-1);
        if deltaTheta>270;deltaTheta=deltaTheta-360;end
         dtheta=tand(deltaTheta);
         S(is)=-0.5*F.*ds*D(is).*dtheta;         
        end 
        S=S(:);

% Boundary Flux %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       [TrCoeff1,Leffective1]=calcTransmission(model.L(1),XsurfNow(1),model.dn(1),model);
        qlhb=TrCoeff1*model.qs(2);
       [TrCoeff2,Leffective2]=calcTransmission(model.L(2),XsurfNow(end),model.dn(end),model); 
        qrhb=TrCoeff2*model.qs(end-1);
        model.TrCoeff=[TrCoeff1, TrCoeff2];
% end boundary flux calc.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute normal displacemet due to LS transport
    [ns, dns_expl]=LaxWendroff(model.ns(:),abs(Qo(:)),qs(:),qlhb,qrhb,ds,dt,dc(:),F,S);
   dns_ds=gradient(model.shadow.*ns,ds);
    Ptemp=Pav2;
    Ptemp=[Ptemp(1)*TrCoeff1,Ptemp,Ptemp(end)*TrCoeff2];
    dPn_ds=gradient(Ptemp,ds);
    clear Ptemp
    dPn_ds(1)=[];
    dPn_ds(end)=[];
   
% Compute normal displacementdue due to cross-shore transport
    dnx=dt.*qn./(Xsurf.*dc); % Wave Power diss

% Sea level rise normal diplacement / equilibrium relaxation model
    timeConstant = (Xsurf.*dc)./(abs(Qo)./Xsurf); % [s]
    model.dnSLR = (dt ./ timeConstant) .* ((SLR0-SLR)./tanb - model.dnSLR); 

% Kalman filter in here
if KF(1).used==1 % do KF stuff

j1=1; % Counter for each data point used (j1) 
Rn=[];
for ikf=1:calibration.nProfiles 
    iProf=KF(ikf).index;
    % 1. State (NB evolutioon ofState is not correct for terms 1-2
    State1=model.nx(iProf)+dnx(iProf);
    State2=ns(iProf);
    State(:,k)=[State1 State2 model.k1 model.k2 model.a1 model.a2 model.b1 model.b2]';

    % 2. Jacobian Matrix
    C=constant.k1/dc(iProf);
    J = ...
            [1 + C*State(3,k)*State(5,k)*Dbar*dt/Xsurf(iProf), 0, -C*dt*(Pn(iProf)-State(5,k)*Dbar*State1-model.b1.*Pav1(iProf))/Xsurf(iProf) , 0, C*State(3,k)*dt*Dbar*State1/Xsurf(iProf), 0, State(3,k)*C*dt*Pav1(iProf)/Xsurf(iProf), 0;...
             0, 1-C*model.shadow(iProf)*State(4,k)*State(6,k)*Dbar*dt/ds, 0, dns_expl(iProf)./State(4,k), 0, -model.shadow(iProf)*State(4,k)*C*dt*Dbar*dns_ds(iProf), 0, -model.shadow(iProf)*C*State(4,k)*dt*dPn_ds(iProf);...
             0, 0, 1, 0, 0, 0, 0, 0;...
             0, 0, 0, 1, 0, 0, 0, 0;...
             0, 0, 0, 0, 1, 0, 0, 0;...
             0, 0, 0, 0 ,0, 1, 0, 0;...
             0, 0, 0, 0 ,0, 0, 1, 0;...
             0, 0, 0, 0 ,0, 0, 0, 1];  
 
        KF(ikf).J=J;

    % 3. Covariance    
        KF(ikf).P= KF(ikf).J * KF(ikf).P * KF(ikf).J' + KF(1).Q;
 
        profileNow=profileIndex{tm}(j1);
 
 if date(k)>=KF(1).trigger && tm+1 < nSurveys && ikf==profileNow && date(k)>KF(1).startDate && date(k)<KF(1).endDate   

% record old state before adjustment
    State_old=State(:,k);

%% Insert Kalman update here...
    % 4. Compute the Kalman Gain
    K = KF(ikf).P * KF(1).H' * (KF(1).H * KF(ikf).P * KF(1).H' + KF(1).R(ikf)).^-1;

if isfield(KF,'opt')
    % 4.5 Compute x,y values
     xp=real(calibration.MCLArray_xy(tm,profileNow));
     yp=imag(calibration.MCLArray_xy(tm,profileNow));
     [Rn(j1), compass_bearing] = distance_to_line_from_point([xp,yp], [model.x(:), model.y(:)]);
     if abs(compass_bearing-model.theta0(iProf))<90
         Rn(j1)=-Rn(j1);
     end
   

    % 5. Update the state vector using the Kalman gain (K) and (x,y) point
    State(:,k) = State(:,k) + K * Rn(j1);

else
    % 5. Update the state vector using the Kalman gain (K) and normal
    % displacement
    Rn(j1)=(calibration.MCLArray(tm,profileNow) - KF(1).H * State(:,k) );
    State(:,k) = State(:,k) + K * Rn(j1);
end

    % 6. Re-calculate the error covariance matrix
    KF(ikf).P = (KF(1).I - K * KF(1).H) * KF(ikf).P;
 

 % 6.5 Impose limits if field 'limits' is defined
 if isfield(KF,'limits')
    State(1,k)=min(max(State(1,k),KF(1).limits(1,1)),KF(1).limits(1,2));
    State(2,k)=min(max(State(2,k),KF(1).limits(2,1)),KF(1).limits(2,2));
    State(3,k)=min(max(State(3,k),KF(1).limits(3,1)),KF(1).limits(3,2));
    State(4,k)=min(max(State(4,k),KF(1).limits(4,1)),KF(1).limits(4,2));
    State(5,k)=min(max(State(5,k),KF(1).limits(5,1)),KF(1).limits(5,2));
    State(6,k)=min(max(State(6,k),KF(1).limits(6,1)),KF(1).limits(6,2));
    State(7,k)=min(max(State(7,k),KF(1).limits(7,1)),KF(1).limits(7,2));
    State(8,k)=min(max(State(8,k),KF(1).limits(8,1)),KF(1).limits(8,2));
 end

     % 7. unpack state matrix
    dnx_corr(j1)=State(1,k)-State_old(1);
    dns_corr(j1)=State(2,k)-State_old(2);
    index_corr(j1)=iProf;
    k1(j1)=State(3,k);
    k2(j1)=State(4,k);
    a1(j1)=State(5,k);
    a2(j1)=State(6,k); 
    b1(j1)=State(7,k);
    b2(j1)=State(8,k); 

    if KF(ikf).Flag==0 || model.shadow(iProf)==0 % Check for shadows or data point too far away from shoreline
    % Reject point
        dnx_corr(j1)=0;
        dns_corr(j1)=0;
        k1(j1)=NaN;
        k2(j1)=NaN;
        a1(j1)=NaN;
        a2(j1)=NaN; 
        b1(j1)=NaN;
        b2(j1)=NaN; 
    end % end if for missing out data

    j1=j1+1;


% Update time of next measurement
if j1>nCal(tm)
    j1=1; % Counter for each data point used (j1) % Emily edit
    tm=tm+1;
    KF(1).trigger=KF(1).timeMeasurements(tm);
    

    % Unpack model corrections for the longshore and cross-shore models and
    % spatially interpolate to other grid points + Find global means
    if exist('dnx_corr')
        dns_corr_all=interp1(index_corr,dns_corr,1:model.Ns,'nearest','extrap');
        dnx_corr_all=interp1(index_corr,dnx_corr,1:model.Ns,'nearest','extrap');
        dns_corr_all=smooth(dns_corr_all,3)';
        dnx_corr_all=smooth(dnx_corr_all,3)';
        model.k1=nanmean(k1);
        model.k2=nanmean(k2);
        model.a1=nanmean(a1);
        model.a2=nanmean(a2);
        model.b1=nanmean(b1);
        model.b2=nanmean(b2);

        % Display updated values
        clc
        disp('Calibration update...')
        disp('.....................')
        disp(['model.k1 = ', num2str(model.k1)])
        disp(['model.k2 = ', num2str(model.k2)])
        disp(['model.a1 = ', num2str(model.a1)])
        disp(['model.a2 = ', num2str(model.a2)])
        disp(['model.b1 = ', num2str(model.b1)])
        disp(['model.b2 = ', num2str(model.b2)])
        disp(['model error = ', num2str(nanmean(Rn)),' [m]'])
        disp(['Next model update = ',datestr(KF(1).trigger)])
        disp('.....................')
        Rn=[];
        clear dnx_corr dns_corr index_corr a1 a2 k1 k2 b1 b2
    else
         dns_corr_all=0*dns_corr_all;
         dnx_corr_all=0*dnx_corr_all;
    end

    % Plot calibration points
    plotfigureOne(model, structures, hard, date(k) )
    figure(1)
    plot(real(calibration.MCLArray_xy(tm,:)),imag(calibration.MCLArray_xy(tm,:)),'r*');

    % Plot State parameters
            fig4 = figure(4);  
            fig4.Units = 'normalized';
            fig4.OuterPosition = [0 0 0.4 1]; % left bottom width height
            conf_level = 0.95; % Set your confidence level
            conf_factor = norminv((1 + conf_level) / 2);

            clf

            upLim = zeros(calibration.nProfiles, length(model.date));
            lowLim = zeros(calibration.nProfiles, length(model.date));

            for l = 1:nStatePlot   % state parameter
                state_values = State(l + 2, 2:length(model.date) + 1);
                for it = 1:calibration.nProfiles
                    % Calculate confidence limits for the current iteration
                    P_matrix = KF(it).P;
                    P_diag = sqrt(abs(diag(P_matrix))); % Extract diagonal elements and take square root
                    upLim(it, :) = state_values + conf_factor * P_diag(l+2);
                    lowLim(it, :) = state_values - conf_factor * P_diag(l+2);
                end

                % Calculate average confidence limits over profiles
                avg_upLim = mean(upLim, 1);
                avg_lowLim = mean(lowLim, 1);

                avg_upLim2(l,:) = avg_upLim(end);
                avg_lowLim2(l,:) = avg_lowLim(end);
            end

            evolve_upLim = [evolve_upLim avg_upLim2];
            evolve_lowLim = [evolve_lowLim avg_lowLim2];
            evolve_date = [evolve_date model.date(end)];

            for l = 1:nStatePlot   % state parameter
                subplot(nStatePlot, 1, l); % final digit = lowercase L not #1
                state_values = State(l + 2, 2:length(model.date) + 1);


                plot(model.date, state_values, 'r','linewidth', 3);
                hold on;

                % Plot evolving average confidence limits
                fill([evolve_date(2:end), fliplr(evolve_date(2:end))], [evolve_upLim(l,2:end), fliplr(evolve_lowLim(l,2:end))], 'red', 'FaceAlpha', 0.3);

                ylabel(stateName(l + 2));
                if l == nStatePlot
                    datetick;
                    xlabel('Date');
                     grid on;
            set(gca, 'FontSize', 16, 'Linewidth', 2);

                else
                    datetick;
                    set(gca, 'XTickLabel', []);
                    grid on;
                    set(gca, 'FontSize', 16, 'Linewidth', 2);

                end
            end

end


 

end % end if

end % end profile loop
 
        State(3,k)=model.k1;
        State(4,k)=model.k2;
        State(5,k)=model.a1;
        State(6,k)=model.a2;
        State(7,k)=model.b1;
        State(8,k)=model.b2;

end % End Kalman Filter (if) 

        % Accumnulate normal displacement components
        model.nx=model.nx+dnx + dnx_corr_all;
        model.ns=ns'+ dns_corr_all;

        model.dn = model.ns + model.nx + model.dnSLR;
        model.x=model.x0+model.dn.*sind(model.theta0);
        model.y=model.y0+model.dn.*cosd(model.theta0);
        dns_corr_all=0*dns_corr_all;
        dnx_corr_all=0*dnx_corr_all;


    clear  P Ps Pn Po ns dnSLR dnx D qn qs


%% File Output Prep %%%%%%%%
    if nextSave<=date(k)
        x=[];
        y=[];
        n=[];
        s=[];
        for i=1:model.Nseg
        x=[x,model.x];
        y=[y,model.y];
        n=[n,model.dn];
        x=[x,NaN];
        y=[y,NaN];
        n=[n,NaN];
        end
    end
    iSave=iSave+1;
    model.date(iSave)=date(k);
    model.xArray(:,iSave)=x;
    model.yArray(:,iSave)=y;
    model.nArray(:,iSave)=n;
    jj=~isnan(x);
    model.s=nan(size(n));
    [model.s(jj)]=computeAlongshoreDist(x(jj),y(jj));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    %% Plot figure 1 Data 
    if nextPlot<=date(k)

    plotfigureOne(model, structures, hard, date(k) )

    % update plotStep
    nextPlot=nextPlot+model.plotStep;


    end % end if for plotting


end % end time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot summary output
    for jj=2:length(model.s)-1; if isnan(model.s(jj));model.s(jj)=(model.s(jj-1)+model.s(jj+1))/2;end;end
    if isnan(model.s(end));model.s(end)=model.s(end-1)+model.ds;end
    figure(2)  
    contourf(model.date,-model.s,model.nArray)
    datetick
    colorbar
    colormap(parula)
    grid on
    box on
    title('Shoreline change with time [m]')
    xlabel('Time')
    ylabel('Alongshore Distance [m]')
    box on
    set(gca,'LineWidth',3,'FontSize',24)

    
% Save final output data
% Record State values in odel array is generated
if KF(1).used==1 
    model.State=State;
end
% Create output directory if not already there
if exist('./outputData','dir')~=7
    mkdir('./outputData')
end

% Test if structures are used and then save
calibration.title='blank structure';
if exist('structures')
    save(filename.output,'model','wave','shoreline','hard','structures',"calibration")
else
    save(filename.output,'model','wave','shoreline','hard','calibration')
end
 
return % End Program >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%%%%%%% Functions....%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model]=find_segments(Xhead,Threshold,model)

% Test for circular boundary
if isfield(model,'Boundary')
    model.Nseg=1;
    model.SegIndex=[1,model.Ns];
    return
end

if Xhead(1)>Threshold;Nseg=0 ;else Nseg=1; istart(1)=1;end
ifin=[];

for i=1:model.Ns-1
    if Xhead(i) > Threshold && Xhead(i+1)<=Threshold
        Nseg=Nseg+1;
        istart(Nseg)=i+1;
    end
    if Xhead(i) <=Threshold && Xhead(i+1)>Threshold
        ifin(Nseg)=i;
    end
end

if length(ifin)<length(istart);ifin(Nseg)=model.Ns;end

model.Nseg=Nseg;
model.SegIndex=[istart(:),ifin(:)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model]=getBeach(shoreline,model)

for iseg=1:model.Nseg
    istart=model.SegIndex(iseg,1);
    ifin=model.SegIndex(iseg,2);
    model(iseg).x=shoreline.x(istart:ifin);
    model(iseg).y=shoreline.y(istart:ifin);
    model(iseg).x0=shoreline.x(istart:ifin);
    model(iseg).y0=shoreline.y(istart:ifin);

    model(iseg).Ns=length(istart:ifin);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,hard]=getNonMobile(shoreline,model)
% Now find non-mobile shoreline
k=0; 

if model.Nseg==1
    hard(1).x=shoreline.x(1:model.SegIndex(1,1));
    hard(1).y=shoreline.y(1:model.SegIndex(1,1));
    hard(2).x=shoreline.x(model.SegIndex(end,end):end);
    hard(2).y=shoreline.y(model.SegIndex(end,end):end);
    model.L(1)=shoreline.Xhead(1);
    model.L(2)=shoreline.Xhead(end);
    if length(hard(2).x)>1 && length(hard(1).x)>1;hard(1).Nseg=2;else hard(1).Nseg=1;end   
return
end



if model.SegIndex(1,1)~=1 
    istart=1;
    ifin=model.SegIndex(1,1);
    hard(1).x=shoreline.x(istart:ifin);
    hard(1).y=shoreline.y(istart:ifin);
    % define headland length scale L for lh boundary
    model.L(1)=max(shoreline.Xhead(istart:ifin));
    k=k+1;
else 
    model.L(1)=0;
end

%% Loop to do internal boundary conditions
for iseg=1:model.Nseg-1

    istart=model.SegIndex(iseg,2);
    ifin=model.SegIndex(iseg+1,1);
    hard(iseg+k).x=shoreline.x(istart:ifin);
    hard(iseg+k).y=shoreline.y(istart:ifin);

    % Need to compute L for internal boundaries
    Ls=hypot(shoreline.x(istart)-shoreline.x(ifin),shoreline.y(istart)-shoreline.y(ifin));
    a=polyarea(shoreline.x(istart:ifin),shoreline.y(istart:ifin));

    if k==0
        model(iseg).L(2)=a/Ls;
        model(iseg+1).L(1)=a/Ls;
    else
        model(iseg).L(2)=a/Ls;
        model(iseg+k).L(1)=a/Ls;
    end
end

Ns=length(shoreline.x);
if model.SegIndex(end,2)~=Ns 
    istart=model.SegIndex(end,2);
    ifin=Ns;
    
    hard(iseg+k+1).x=shoreline.x(istart:ifin);
    hard(iseg+k+1).y=shoreline.y(istart:ifin);
    % define headland length scale L for end boundary
    model(iseg+k).L(2)=max(shoreline.Xhead(istart:ifin));
else
    model(iseg+1).L(2)=0;
end

    hard(1).Nseg=length(hard);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xval,yval,ds,Ns]=initialiseShoreline(x,y,Ns,ds)

% compute alongcoast distance on orginal shoreline
[s1]=computeAlongshoreDist(x(:),y(:));

shorelineLength=max(s1)-min(s1);

% Find new modified ds
if nargin==3
    ds=shorelineLength./(Ns-1);
else
    Ns=floor(shorelineLength/ds);
    ds=shorelineLength./Ns;
end

% recompute s
s2= min(s1):ds:max(s1);
Ns=length(s2);

% reset initial values
if length(unique(x))~=length(x) || length(unique(y))~=length(y)
    %disp('yep')
    x=x+0.001*rand(size(x));
    y=y+0.001*rand(size(y));
    s1=s1+0.001*rand(size(s1));
end


yval=interp1(s1,y,s2,'spline','extrap');
xval=interp1(s1,x,s2,'spline','extrap');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s]=computeAlongshoreDist(x,y) 
x=x(:);y=y(:);
% Compute distance alongshore (S)
dx=x(1:end-1)-x(2:end);
dy=y(1:end-1)-y(2:end);
s=sqrt(dx.^2+dy.^2);s=[0;s];s=cumsum(s);s=s(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta]=findNormal(x,y)
% This version treats boundary values differently by first order diffrence
% of internal values
% Extra mod to deal with just two points by returnin a single value.
 
    N=length(x);
    if N==2
        im=1;
        ip=2;
    else
        im=[1,1:N-2,N-1];
        ip=[2,3:N,N];
    end

    xDiff=x(im)-x(ip);
    yDiff=y(im)-y(ip);
    theta=atan2d(yDiff,-xDiff);

    jj=find(theta<0);
    theta(jj)=theta(jj)+360;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [distance, compass_bearing] = distance_to_line_from_point(point, line_points)

% First find distance to all points
dist=hypot(line_points(:,1)-point(1), line_points(:,2)-point(2));
jj=find(dist==min(dist)); jj=max(jj);

% Compute dx and dy
dx=line_points(jj,1)-point(1);
dy=line_points(jj,2)-point(2);

% Get min distance
distance=dist(jj);

% Compute angle rel. x-axis
theta=atan2d(-dx,-dy);

% convert to compass N
compass_bearing = mod(theta, 360); % Convert mathematical angle to compass bearing (from line to point)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DOC,H12,T12]=Hallermeier(H,T,dt)
% deep water H array [m]
% deep water T array [s]
% wave sampling interval (must be regular)
% DOC - depth of closure [m]
% H12 - 12 hour in one year wave height [m]
% T12 - same for T [s]
% dt is the wave time-step in days

if length(H)<100
    DOC=2*max(H);
    H12=max(H);
    T12=max(T);
    return
end

[H,I]=sort(H,'descend');
T=T(I);
Npts=length(H);
Nyears=(Npts*dt)/365.25;
indexVal=round(Nyears*0.5/365.25);
indexVal=max(indexVal,1);
H12=H(indexVal);
T12=T(indexVal);
DOC= 2.28 * H12 - 68.5 * (H12.^2/ (9.81*T12.^2));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simple Run-up routine
function [R2,setup,swash]=RunUpSimple(Ho)
 setup=0.17*Ho;
 swash=0.35*Ho;
 R2=setup+swash;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R2,setup,swash]=StockdonHolman(Ho,T,beta,swashFact)
g=9.98;
Lo=(g*T.^2)/(2*pi);
setup=swashFact * 0.35 * beta * sqrt(Ho.*Lo);
swash=swashFact * 0.5 * sqrt(Ho .* Lo * (0.563 * beta.^2 + 0.004));
R2=(setup+swash);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phiNew]=upwind_correctionMD(phi)
phimax=45;
n=length(phi);
phiNew=phi;

% Note that we leave boundary values the same and go from 2 to n-1

for i=2:n
    
    if phi(i-1)>0 && phi(i)>0 && ...
            abs(phi(i-1))<=phimax && abs(phi(i))>phimax
        phiNew(i)=phimax;
    end

 end

 for i=n-1:-1:1   
    if  phi(i+1)<0 && phi(i)<0 &&...
            abs(phi(i+1))<=phimax && abs(phi(i))>phimax
        phiNew(i)=-phimax;
    end
 end

 jj=find(abs(phiNew)>90);
 phiNew(jj)=sign(phiNew(jj))*90;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [n,dn]=LaxWendroff(n,Qo,q,qlhb,qrhb,ds,dt,dc,F,S)

%%%% Numerical implicit-explicit solution%%%%%%%
nx=length(n);

% Initialise New Term
 a=zeros(size(q));
 b=a;c=a;d=a;dn=a;

% Diffusson coeff.
 D=Qo*dt./(dc*ds.^2);

 %% define tridiagonal coefficients a-c and known values d
 for i=2:nx-1
     a(i)=-F*D(i-1);
     b(i)= F*2*D(i)+1;
     c(i)=-F*D(i+1);
     dn(i)= -(dt/(2*ds*dc(i))) * (q(i+1)-q(i-1));
     d(i)=n(i) + (1-F) * dn(i) + S(i);
 end 

%  % boundary values
 a(1)=0;   b(1)=1;    c(1)=0;  
 a(nx)=0;  b(nx)=1;   c(nx)=0;
 dn(1) =   - (dt/(2*ds*dc(1))) * (q(2) - qlhb);
 dn(nx)=   - (dt/(2*ds*dc(nx))) * (qrhb - q(nx-1));    
 d(1) =   n(1)     + dn(1);
 d(nx)=   n(nx)    + dn(nx);   
if F==0;return;end
  %% Find an implicit solution over the whole domain in one step by inverting a tridiagonal matrix
  % Note new function at the end
   [n]=implicit(a,b,c,d);
 end

 %% Find solution via matrix inversion 
function [y]=implicit(a,b,c,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implicit solution to tridiagonal equations                            
% Input:                                                                                                               
% a, b, c implicit inversion coefficients.                              
% S: Source function f(x).                                              
% Output:                                                               
% y are approximate solution at grid points                             
% Where: a*y(i-1)+b*y(i)+c*y(i+1)=S     
% Written by Mark Davidson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 S=S(:);a=a(:);b=b(:);c=c(:);
 [n] = length(S);
 A = sparse(n,n);

 A(1,[1 2])       = [b(1) c(1)]; % lhs boundary
for i=2:n-1
 A(i,[i-1 i i+1]) = [a(i) b(i) c(i)];
end   
 A(n,[n-1 n])     = [a(n) b(n)]; % rhs boundary
 
 y = A\S;
 
end

% function [y]=analyticSolution(x,K,tval,theta)
% % Analytical solution
% y=tand(theta).*sqrt((4*K*tval)./pi).*...
%     (exp(-x.^2/(4*K*tval))-((x*sqrt(pi)/(2*sqrt(K*tval)))).*erfc(x./(2*sqrt(K*tval))));
% end
 
%% Running mean algorithm %%%%%%%%%%%%%%%%%%%%%%%
function [nEq]=runningMean(nPointsAveraged,nEq,n)
 nEq = nEq + ((n - nEq)/nPointsAveraged);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ xS,yS,shadowS ] = find_shadows( x,y,x_mc,y_mc,phiw,hard )
%
%  
n=length(x)-1;
if hard==1
    crit=0;
else
    crit=0;
end
if n==0
    xS=[];
    yS=[];
    shadowS=[];
else
    f=180/pi;
    xS=.5*(x(1:n)+x(2:n+1));
    yS=.5*(y(1:n)+y(2:n+1));
    len=5*hypot(max(x_mc)-min(x_mc),max(y_mc)-min(y_mc));
    
    for i=1:n
        xw=[xS(i)+1*sin(phiw),xS(i)+len*sin(phiw)];
        yw=[yS(i)+1*cos(phiw),yS(i)+len*cos(phiw)];
        P1=InterX([x_mc;y_mc],[xw;yw]);
        
        shadowS(i)=size(P1,2)>crit;
        if 0
            figure(10)
            plot(x,y,x_mc,y_mc,xw,yw,'.-',P1(1,:),P1(2,:),'ro','linewidth',2);
            axis equal
            drawnow
            pause
        end
    end
end

end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve 
%   together with any self-intersection points.
%   
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 3.0, 21 Sept. 2010

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

    %...Argument checks and assignment of L2
    %error(nargchk(1,2,nargin));
    if nargin == 1,
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TrCoeff,Leffective]=calcTransmission(L,Xsurf,dn,model)

if isfield(model,'kTr')
    kTr=model.kTr;
else
    kTr=1;
end

Leffective=kTr.*(L-dn);  
if L==0;Leffective=0;end
Leffective=max(Leffective,0);
ratio  = Leffective./Xsurf;


if isfield(model,'TrExp')
    if isfield(model,'p')~=1
        p=2;
    else
        p=model.p;
    end
 TrCoeff=exp(-ratio.^p);
else
 TrCoeff = 1 - ratio;
end

TrCoeff=max(0,TrCoeff);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [weight,wave]=computeWeights(x,y,wave)
wave(1).N=length(wave);

if wave(1).N>1

for i = 1:wave(1).N
    distance(:,i)=hypot(x-wave(i).x(1), y-wave(i).y(1));
end

weight=1./distance.^2;
sumWeight=sum(weight');
weight=weight./sumWeight(:);
else 
    weight=ones(size(x));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [O]=applyWeight(A,weight)
% A is a column matrix, where each column represents a different wave gauge
% A can contain any wave parameter, H, Tp or Dir
% Weight is a row matrix, where each value represents the weighting on a
% particular wave gauge.
% O is the weighted output of A

    % Apply weights via matrix dot product...
O = A * weight';

end


%%%%%%%%
function [E, N] = ll2os( lat, lon )
% [E, N] = LL2OS( LAT, LON ) returns UK National Grid coodinates E, N 
% (easting and northing) for points given by latitude and longitude 
% Convert latitude/longitude => OS National Grid Reference points 
% algorithm and constants adapted from D00659 v2.1 Dec 2010 of
% http://www.ordnancesurvey.co.uk/oswebsite/gps/docs/A_Guide_to_Coordinate_Systems_in_Great_Britain.pdf
% (c) Michael Fourman 2012 

  phi    = lat .* (pi/180); % convert arguments to radians
  lambda = lon .* (pi/180);
  
  a = 6377563.396;
  b = 6356256.909;              % Airy 1830 major & minor semi-axes
  F0 = 0.9996012717;            % NatGrid scale factor on central meridian
  phi0 = 49*pi/180;             % \phi_0 
  lambda0 = -2*pi/180;          % \lambda_0
  N0 = -100000;                 % NatGrid true origin 
  E0 =  400000;                 % northing & easting of true origin, metres
  e2 = 1 - b*b/(a*a);           % eccentricity squared
  n = (a-b)/(a+b);              % C1
  
  e2sin2phi = e2 * sin(phi).^2;
  nu  = a * F0 * (1-e2sin2phi).^(-0.5);
  rho = a * F0 * (1 - e2).*(1-e2sin2phi).^(-1.5);
  
  eta2 = (nu./rho) - 1 ;         %C2
  
  M = b .* F0 .* ( ...
      (1 + n + 1.25 .* (n^2 + n^3)) .* (phi - phi0) - ...
      3 .* (n + n^2 + (0.875 * n^3)) .* sin(phi - phi0) .* cos(phi + phi0) + ...
      1.875 .* (n^2 + n^3) .* sin(2 .* (phi - phi0)).* cos(2 .* (phi + phi0)) - ...
      (35/24) .* n^3 .* sin(3 .* (phi - phi0)) .* cos(3 .* (phi + phi0)) ...
      ) ;                       %C3
  I   = M + N0;
  II  = (nu / 2) .* sin(phi) .* cos(phi);
  III = (nu / 24) .* sin(phi) .* cos(phi).^3 .* (5 - tan(phi).^2 + (9 * eta2));
  IIIA = (nu./720) .* sin(phi) .* cos(phi).^5 .* ...
      (61 - 58 .* tan(phi).^2 + tan(phi).^4);
  IV = nu .* cos(phi);
  V  = (nu ./ 6) .* cos(phi).^3 .* (nu./rho - tan(phi).^2);
  VI = (nu ./ 120) .* cos(phi).^5 .* ...
      (5 - 18 .* tan(phi).^2 + tan(phi).^4 + (14 - 58 .* tan(phi).^2).* eta2);
  N = I + II.* (lambda - lambda0).^2 + IIIA .* (lambda - lambda0).^4 ; %C4
  E = E0 + IV .* (lambda - lambda0) + V .* (lambda - lambda0).^3 + VI .* (lambda - lambda0).^5;  %C5
end

%%%%%%%%%%%%%
% osgb36wgs84.m        ; % (c) dmitry.aleynik@sams.ac.uk 
% Date:    07.01.2011   ; matlab script; !da updated 2018.04.11
% Project: Asimuth EU-FP7: algae bloom monitoring
% based  : on java code Chris Veness 2005-2010 :
% http://www.movable-type.co.uk/scripts/latlong-convert-coords.html

function [lon, lat] = osgb36wgs84(rlon,rlat,key)
%key == 0: o -> w 
%key ~= 0: w -> o 
% ellipsoides partameters       e:
        e.WGS84.a    = 6378137         ; 
        e.WGS84.b    = 6356752.3142    ; 
        e.WGS84.f    = 1/298.257223563 ;
        
        e.Airy1830.a = 6377563.396   ; 
        e.Airy1830.b = 6356256.910   ; 
        e.Airy1830.f = 1/299.3249646 ;

% helmert transform parameters h:
  h.WGS84toOSGB36.tx = -446.448 ; h.WGS84toOSGB36.ty = 125.157; h.WGS84toOSGB36.tz = -542.060 ;% m
  h.WGS84toOSGB36.rx = -0.1502  ; h.WGS84toOSGB36.ry = -0.2470; h.WGS84toOSGB36.rz = -0.8421  ;% sec
  h.WGS84toOSGB36.s  = 20.4894  ; % ppm
  
  h.OSGB36toWGS84.tx = 446.448  ; h.OSGB36toWGS84.ty = -125.157; h.OSGB36toWGS84.tz = 542.060 ;
  h.OSGB36toWGS84.rx = 0.1502   ; h.OSGB36toWGS84.ry =   0.2470; h.OSGB36toWGS84.rz =  0.8421 ;
  h.OSGB36toWGS84.s = -20.4894  ;

         p1.lon    =rlon;
         p1.lat    =rlat;        
         p1.height =5   ;
 if (key==0 ),...         
% convertOSGB36toWGS84(p1);
          p2 = convert(p1, e.Airy1830, h.OSGB36toWGS84, e.WGS84);
 else          
% convertWGS84toOSGB36(p1)
          p2 = convert(p1, e.WGS84, h.WGS84toOSGB36, e.Airy1830);
 end
          lon=p2.lon;
          lat=p2.lat;
end

% end

 function p2 = convert(p1, e1, t, e2) 
  % -- convert polar to cartesian coordinates (using ellipse 1)
  % p1.lon =rlon;  p1.lat =rlat;
  p1.lat = p1.lat*pi/180; 
  p1.lon = p1.lon*pi/180; 

  a = e1.a; b = e1.b ;

  sinPhi    = sin(p1.lat); 
  cosPhi    = cos(p1.lat);
  sinLambda = sin(p1.lon); 
  cosLambda = cos(p1.lon);
          H = p1.height;

   eSq = (a*a - b*b) / (a*a);
   nu = a ./ sqrt(1 - eSq*sinPhi.*sinPhi);

   x1 = (nu+H) .* cosPhi .* cosLambda;
   y1 = (nu+H) .* cosPhi .* sinLambda;
   z1 = ((1-eSq)*nu + H) .* sinPhi;


  % -- apply helmert transform using appropriate params
  
   tx = t.tx; ty = t.ty; tz = t.tz;
   rx = t.rx/3600 * pi/180;  % normalise seconds to radians
   ry = t.ry/3600 * pi/180;
   rz = t.rz/3600 * pi/180;
   s1 = t.s/1e6 + 1;              % normalise ppm to (s+1)

  % apply transform
   x2 = tx + x1*s1 - y1*rz + z1*ry;
   y2 = ty + x1*rz + y1*s1 - z1*rx;
   z2 = tz - x1*ry + y1*rx + z1*s1;


  % -- convert cartesian to polar coordinates (using ellipse 2)

  a = e2.a; b = e2.b;
   precision = 4 / a;  % results accurate to around 4 metres

  eSq   = (a*a - b*b) / (a*a);
   p    = sqrt(x2.*x2 + y2.*y2);
   phi  = atan2(z2, p*(1-eSq));
   phiP = 2*pi;
   
  while (abs(phi-phiP) > precision) ,...
    nu = a ./ sqrt(1 - eSq*sin(phi).*sin(phi));
    phiP = phi;
    phi = atan2(z2 + eSq*nu.*sin(phi), p);
  end
  
   lambda = atan2(y2, x2);
        H = p./cos(phi) - nu;
        
  p2.lat    = phi*180/pi;
  p2.lon    = lambda*180/pi;
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [model]=computeAverageValues(wave,model)
% Function to compute depth of closure (dc), surfzone width (Xsurf) and
% average slope (tanb) for all longshore points given local wave data in
% wave structure.


%% Compute the average depth of closure surfzone widths and slopes for all locations
for i=1:length(wave)
    [Ho]=HoLinear(wave(i).H,wave(i).Tp,wave(i).depth,inf);
    [DOC(i),H12(i),T12(i)]=Hallermeier(Ho,wave(i).Tp,wave(i).dt);
    [R2(i),setup,swash]=RunUpSimple(mean(Ho));
end

% Compute longshore varying values and add to model structure
model.DOC=applyWeight(DOC,model.weights);
model.Berm=applyWeight(R2,model.weights);
model.dc=applyWeight(DOC+R2+model.MSR,model.weights);
model.Xsurf=(model.dc./model.A).^(3/2);
model.tanb=model.dc./model.Xsurf;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model]=computeAveragePowerComponents(wave,model,dates)
if length(dates)>1
    disp(' Computing average power values - may take a while!')
    tic;
end

%% New time-scale with daily sampling
if length(dates)>1
    newdate=dates(1):dates(end);
else
    newdate=dates(1);
end
%% Number of points to be averaged
Ntime=round(min(length(newdate),365.25*5/wave(1).dt));
Ntime=max(1,Ntime);

%% Initialise arrays
HsArray=zeros(length(model.x0),Ntime);DirArray=HsArray;DirArrayAbs=DirArray;Po=HsArray;P=Po;Pn=Po;Ps=Po;DirbArray=DirArray;Hb=Po;hb=Po;

% Interpolate waves data for each gauge (ig)
for k=1:Ntime

for ig=1:wave(1).N
    HsG(ig)=interp1(wave(ig).dates,wave(ig).H,newdate(k),'linear','extrap');
    Tp(k)=interp1(wave(1).dates,wave(1).Tp,newdate(k),'linear','extrap');
    DirG(ig)=interp1(wave(ig).dates,wave(ig).Dir,newdate(k),'linear','extrap'); % NB Need to shoal this array to breaking
end


[HsArray(:,k)]=applyWeight(HsG,model.weights);
[DirArrayAbs(:,k)]=applyWeight(DirG,model.weights);
[Ho]=HoLinear(HsArray(:,k),Tp(k),wave(1).depth,inf);
[Hb(:,k)]=HbKomar(Ho,Tp(k));
hb(:,k)=Hb(:,k)./0.78;

% Define spatial arrays of wave parameters along the coastal segment
% Angle rel. to coast (care to give angle hense transport the right sign here)
DirbArray(:,k)=model.theta(:)-DirArrayAbs(:,k);

% Adjust angle limits
jj=find(DirbArray(:,k) > 180);
DirbArray(jj,k)=DirbArray(jj,k)-360;
jj=find(DirbArray(:,k) <-180);
DirbArray(jj,k)=DirbArray(jj,k)+360;

 % Linear-shoal wave direction
[c]=cLinear(Tp(k),hb(:,k));
[c0]=cLinear(Tp(k),wave(1).depth*ones(size(hb(:,k))));
DirbArray(:,k)=asind(c.*sind(DirbArray(:,k))./c0);

%% Compute P0, Ps and Pn
[c]=cLinear(Tp(k),wave(1).depth);
[wn]=kLinear(Tp(k),wave(1).depth);
[E]=ELinear(HsArray(:,k));
[nval]=nLinear(wn,wave(1).depth);
[P(:,k)]=PLinear(E,c,nval);
Po(:,k)=P(:,k).* cosd(DirbArray(:,k));
Ps(:,k)=Po(:,k).*sind(DirbArray(:,k));
Pn(:,k)=Po(:,k).*cosd(DirbArray(:,k));

end

model.PMean=sum(P,2)/Ntime;model.PsMean=model.PMean';
model.PsMean=sum(Ps,2)/Ntime;model.PsMean=model.PsMean';
model.PnMean=sum(Pn,2)/Ntime;model.PnMean=model.PnMean';
model.PoMean=sum(Po,2)/Ntime;model.PoMean=model.PoMean';
model.DirMean=sum(DirbArray,2)/Ntime;model.DirMean=model.DirMean';
model.HbMean=sum(Hb,2)/Ntime;model.HbMean=model.HbMean';
model.hbMean=sum(hb,2)/Ntime;model.hbMean=model.hbMean';
model.DirAbs=sum(DirArrayAbs,2)/Ntime;model.DirAbs=model.DirAbs';

if k>1; disp(' Done.');toc, end

end

function plotfigureOne(model, structures, hard, date )

    fig1 = figure(1);  
    fig1.Units = 'normalized';
    fig1.OuterPosition = [0.4 0 0.6 0.6]; % left bottom width height
    cla
    hold on
    for i = 1:model.Nseg
        plot(model.x,model.y,  'r','LineWidth',3)
        plot(model.x0,model.y0,'g','LineWidth',1)
    end

    for i = 1:hard(1).Nseg
          plot(hard(i).x,hard(i).y,'k','LineWidth',4)
    end
    % Plot structures if present
    if isfield(model,'structures')
        for i = 1: structures(1).Nstruct
            plot(structures(i).x,structures(i).y,'b',LineWidth=4)
        end
    end
    grid on
    box on
    title(['Site: ',model.site,'. Time = ',datestr(date)])
    legend('Modelled shoreline','Original shoreline','location','best')
    xlabel('Distance [m]')
    ylabel('Distance [m]')
    if isfield(model,'axisEqual')
        axis equal
    end
    box on
    set(gca,'LineWidth',3,'FontSize',24)

end

