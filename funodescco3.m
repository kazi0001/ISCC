function out = funodescco3(ui,Params)
% M-File Single-Column Chromatography Process.
% Bijan Medi, SCBE, NTU, 2011.

% Ex3: Replacing DRF with DrF,new cut time scheduling from sfunvlv17 was
% used. Adding PrF.
% odescc1: Removing mad, rhos from Params
% funodescc1: functional form, global vars removed
% 2: Adding a while loop for Ninj
% 3: Using Structures, rearrenging output, OSCC, ui in non-SI units, Params
% in SI units), tbias used

% System Parameters -------------------------------------------------------
L = Params.L; % m
D = Params.D; % m
eb = Params.eb;
Nz = Params.Nz;
tbias = Params.tbias;
tres = Params.tres;
% Vcol = Params.Vcol;
mad = Params.mad;
% ncol = Params.ncol;

rho = Params.rho; % Kg/m3 Liquid density
Mu = Params.Mu;% Pa.s
dp = Params.dp;%  m
KovA = Params.KovA; % m/s
KovB = Params.KovB; % m/s
phi = Params.phi;

% Isotherm --------------------------------------
HA = Params.HA;
KA = Params.KA;
HB = Params.HB;
KB = Params.KB;
% Discrete arrays are directly passed to the C function
% -----------------------------------------------

% INPUTS ==================================================================

Ncy = Params.Ncy;

cref= Params.cref;
% ntot = Params.ntot; % Total number of columns


% % Adsorption ------------------
yA0 = zeros(1,2*Nz); % Normalized
yB0 = zeros(1,2*Nz); % Normalized
% % -----------------------------



ncF = 0;

% tic
while ncF<3 % Updating Ncy to get necessary number of eluted peaks -------

y0 = [yA0 yB0];

% INPUTS ==================================================================
Vinj = ui(1)*ones(Ncy,1)*1e-9; % ul -> m3
tcy = ui(2)*ones(Ncy,1); % sec
Q = ui(3)*ones(Ncy,1)/1e6/60; % ml/min -> m3/s
cF = ui(4)*ones(Ncy,1); % g/l

dtc1 = ui(5)*ones(Ncy,1); % sec
dtc2 = ui(6)*ones(Ncy,1); % sec
dtc3 = ui(7)*ones(Ncy,1); % sec
% =========================================================================

% -------------------------------------------------------------------------
t0 = 0;
t = [];
y = [];
yA = [];
yB = [];

check = 0; % Initialization, determines just the start of collection.

tsc2 = 0;
t2 = 0;
ncF = 0; % Fraction collector cycle counter
data1 = zeros(Ncy,7); % Stores feed and cut times data
DrF = 0; % Desorbent requirement
PrF = 0; % Productivity
dtcF1 = 0;
dtcF2 = 0;
dtcF3 = 0;
dtcF4 = 0;
% -------------------------------------------------------------------------

% Solving cycle to cycle ==================================================

    
for k=1:Ncy
    
    tspan = [t0:tres:tcy(k)];
    cAin = cF(k)/2;
    cBin = cF(k)/2;
    
    Params1 = [Nz,L,D,eb,Q(k),cref,0,HA,KA,HB,KB,cAin,cBin,Vinj(k),tcy(k),rho,Mu,dp,KovA,KovB];
       
    % ODE SOLVER ----------------------------------------------------------
    [t1,y1] = ode113(@funscc3,tspan,y0,Params.odeoptions,Params1);
    % ---------------------------------------------------------------------
    
    y0 = y1(end,:);
    
    if k==1
        t2 = t1;
        t = t2;
        
        yA = y1(:,1:Nz);
        yB = y1(:,2*Nz+1:3*Nz);
    else
        t2 = t2(end) + t1;
        t = [t;t2(2:end)];
        yA = [yA;y1(2:end,1:Nz)];
        yB = [yB;y1(2:end,2*Nz+1:3*Nz)];
    end
    
    
end
% =========================================================================

for k=1:Ncy
    
    
    % DATA READ =======================
    % Vinj, tcy, Qout, and cF are updated only at the
    % end of each nci th cycle. dtci are stored at nci th location and so they
    % are implemented when nci th cycle arrives.
    
    data1(k,1) = Vinj(k);
    data1(k,2) = tcy(k);
    data1(k,3) = Q(k); % ml/s Sent out directly along with cA and cB, but the stored value is used for DrF
    data1(k,4) = cF(k);
    data1(k,5) = dtc1(k); % Stored for respective cycle
    data1(k,6) = dtc2(k); % Stored for respective cycle
    data1(k,7) = dtc3(k); % Stored for respective cycle
    % =================================
    
end % Ncy

k = 1;
temp=0;
fraction =0;

% tct0=[];

% i1 = 1;
temp1 = 0;
tct1 = [];
sys1 = [];

% i2 = 1;
temp2 = 0;
tct2 = [];
sys2 = [];

% i3 = 1;
temp3 = 0;
tct3 = [];
sys3 = [];

% i4 = 1;
temp4 = 0;
tct4 = [];
sys4 = [];

PB = 0;
PA = 0;
YB = 0;
YA = 0;

% ------------------------
clear y y1 y0
% ------------------------

for i=1:size(t,1)-1
    
    yAout = yA(i,end);
    yBout = yB(i,end);
    
    % Allocating solvent flow rate ----------------------------------------
    % This is for variable solvent flow rate.
    if t(i)<tcy(1)
        k=1;
        ts = 0; % Start time
        te = tcy(k); % End time
        % --------------------------
        Qout = Q(k); % m3/s
        % --------------------------
    elseif temp==0
        %         t(i)-t1>=tcy(k) && t(i)-t1<tcy(k+1) && temp==0
        ts = ts+tcy(k); % Beginning of cycle k+1
        te = te+tcy(k+1); % End of cycle k+1
        k = k+1;
        % --------------------------
        Qout = Q(k); % m3/s
        % --------------------------
        temp = 1;
    elseif t(i)>=te
        temp = 0;
        % --------------------------
        Qout = Q(k+1); % m3/s
        % --------------------------
    else
        % --------------------------
        Qout = Q(k); % m3/s
        % --------------------------
    end
    % -------------------------------------------------------------------------
    
    
    uv = yAout + yBout;
    
    if check==1 && dtcF4 <= 0
        error('Error: Cycle time is smaller than or equal to sum of three cut times.')
    end
    
    
    
    % Arrival of the first cycles and start of Collecting B --------------
    if uv > 0.05 && fraction ==0
        fraction =1;
        tsc2 = t(i);
        
        check = 1; % Marks the start of first cycle
        ncF = 1; % Number of cycles arrived at the fraction collector.
        % ----------------------------------------------------------------
        % First cycle time for the first collection. For the first cycle, these variables can be
        % changed after start up and before the arrival of first cycle.
        
        tcyF = data1(ncF,2); %
        dtcF1 = data1(ncF,5); % sec
        dtcF2 = data1(ncF,6); % sec
        dtcF3 = data1(ncF,7);% sec
        % ----------------------------------------------------------------
        
        dtcF4 = tcyF - (dtcF1 + dtcF2 + dtcF3); % sec
        
    elseif fraction==1 && t(i)-tsc2-dtcF1>=tbias % End of Collecting B
        
        fraction = 2;
        
    elseif fraction==2 && t(i)-tsc2-(dtcF1+dtcF2)>=tbias % End of Collecting A+B
        
        fraction = 3;
        
    elseif fraction==3 && t(i)-tsc2-(dtcF1+dtcF2+dtcF3)>=tbias % End of Collecting A
        
        fraction = 4;
        
    elseif fraction==4 && t(i)-tsc2-(dtcF1+dtcF2+dtcF3+dtcF4)>=tbias % End of Collecting Solvent, Arrival of a new cycle
        
        fraction = 1;
        
        tsc2 = t(i); % Marks the start of each cycle at the arrival of elution profile
        
        % DATA WRITE ==========================================
        % Data is written to vars from the respective cycle
        
        VinjF = data1(ncF,1); % With reference to ncF th cycle
        tcyF = data1(ncF,2); % With reference to ncF th cycle
        QDF = data1(ncF,3); % With reference to ncF th cycle
        cFF = data1(ncF,4); % With reference to ncF th cycle
        dtcF1 = data1(ncF,5); % sec
        dtcF2 = data1(ncF,6); % sec
        dtcF3 = data1(ncF,7);% sec
        dtcF4 = tcyF - (dtcF1 + dtcF2 + dtcF3); % sec
        % This group must be always together
        % =====================================================
        
        % DESORBENT REQUIREMENT ==============================
        DrF(ncF) = (QDF*tcyF + VinjF)/VinjF/cFF; % Desorbent requirement with reference to ncF th cycle mlD/mgF
        % QDF is converted to m3/s
        %=====================================================
        
        % Pressure-related parameters -------------------------------------
        
        DelP(ncF) = phi*QDF/(pi*D^2/4)*L*Mu/(dp^2)/1e5;
        % (bar) Darcy equation from Cox's book (L and D in m)
        %  phi=500, resistance  parameter
        %------------------------------------------------------------------
        
        % Calculating Purity and Recovery =====================================
        % A variable step integrator should be used.
        % m=int(Qout(t)*cA(t),t)
        mA1 = trapez(sys1(:,3).*sys1(:,1),tct1); % kg A collected in Fraction 1
        mB1 = trapez(sys1(:,3).*sys1(:,2),tct1); % kg B collected in Fraction 1
        
        mA2 = trapez(sys2(:,3).*sys2(:,1),tct2); % kg A collected in Fraction 2
        mB2 = trapez(sys2(:,3).*sys2(:,2),tct2); % kg B collected in Fraction 2
        
        mA3 = trapez(sys3(:,3).*sys3(:,1),tct3); % kg A collected in Fraction 3
        mB3 = trapez(sys3(:,3).*sys3(:,2),tct3); % kg B collected in Fraction 3
        
        mA4 = trapez(sys4(:,3).*sys4(:,1),tct4); % kg A collected in Fraction 4
        mB4 = trapez(sys4(:,3).*sys4(:,2),tct4); % kg B collected in Fraction 4
        % =====================================================================
        
        PB(ncF,1) = mB1/(mA1 + mB1)*100;
        PA(ncF,1) = mA3/(mA3 + mB3)*100;
        
        YB(ncF,1) = mB1/(mB1 + mB2 + mB3 + mB4)*100;
        YA(ncF,1) = mA3/(mA1 + mA2 + mA3 + mA4)*100;
        
        % PRODUCTVITY ===============================================
        PrF(ncF) = 1e6*(mA1 + mB1 + mA3 + mB3)/tcyF/mad*60; % (kg->mg) mg F/mg Ad./min
        % ===========================================================
        
        %     i1 = 1;
        sys1 = [];
        tct1 = [];
        temp1 = 0;
        
        %     i2 = 1;
        sys2 = [];
        tct2 = [];
        temp2 = 0;
        
        %     i3 = 1;
        sys3 = [];
        tct3 = [];
        temp3 = 0;
        
        %     i4 = 1;
        sys4 = [];
        tct4 = [];
        temp4 = 0;
        
        % Updating ncF --------------------------------
        ncF = min(ncF + 1,Ncy); % Limits ncF to Ncy
        % ---------------------------------------------
        
    end % uv
    
    
    switch fraction
        
        case 0
            
        case 1
            % Finding Major Step ----------------------------------------------
            % This is intended for picking only major steps, though it is not
            % accurate as there are some minor steps forward.
            % The time series are variable step. So, a variable step integrator
            % is required.
            
            if isempty(tct1)
                tct1(1,1) = t(i);
                sys1 = [sys1;[cref*yAout,cref*yBout,Qout]];
            elseif t(i) > tct1(end,1)
                
                temp1 = temp1+t(i);
                sys1 = [sys1;[cref*yAout,cref*yBout,Qout]];
                
                tct1 = [tct1;t(i)];
                %             i1 = i1+1;
            end
            % -----------------------------------------------------------------
            
            
        case 2
            
            % Finding Major Step ----------------------------------------------
            % This is intended for picking only major steps, though it is not
            % accurate as there are some minor steps forward.
            % The time series are variable step. So, a variable step integrator
            % is required.
            
            if isempty(tct2)
                tct2(1,1) = t(i);
                sys2 = [sys2;[cref*yAout,cref*yBout,Qout]];
            elseif t(i) > tct2(end,1)
                
                temp2 = temp2+t(i);
                sys2 = [sys2;[cref*yAout,cref*yBout,Qout]];
                
                tct2 = [tct2;t(i)];
                %                         i2 = i2+1
            end
            % -----------------------------------------------------------------
            
            
        case 3
            
            % Finding Major Step ----------------------------------------------
            % This is intended for picking only major steps, though it is not
            % accurate as there are some minor steps forward.
            % The time series are variable step. So, a variable step integrator
            % is required.
            
            if isempty(tct3)
                tct3(1,1) = t(i);
                sys3 = [sys3;[cref*yAout,cref*yBout,Qout]];
            elseif t(i) > tct3(end,1)
                
                temp3 = temp3+t(i);
                sys3 = [sys3;[cref*yAout,cref*yBout,Qout]];
                
                tct3 = [tct3;t(i)];
                %             i3 = i3+1;
            end
            % -----------------------------------------------------------------
            
            
        case 4
            
            % Finding Major Step ----------------------------------------------
            % This is intended for picking only major steps, though it is not
            % accurate as there are some minor steps forward.
            % The time series are variable step. So, a variable step integrator
            % is required.
            
            if isempty(tct4)
                tct4(1,1) = t(i);
                sys4 = [sys4;[cref*yAout,cref*yBout,Qout]];
            elseif t(i) > tct4(end,1)
                
                temp4 = temp4+t(i);
                sys4 = [sys4;[cref*yAout,cref*yBout,Qout]];
                
                tct4 = [tct4;t(i)];
                %             i4 = i4+1;
            end
            % -----------------------------------------------------------------
    end % switch
end % for t


% Updating Ncy to get necessary number of eluted peaks -------------------
Ncy = Ncy + 1;
% -------------------------------------------------------------------------
end % While ncF

%plot(t,yB(:,end),'b',t,yA(:,end),'r')

% if ncF<3
%     error('The number of collected cycles ncF must be greater than two')
% end
% trun = toc
% OUTPUTS =================================================================
out = [PA(end) PB(end) YA(end) YB(end) PrF(end) DrF(end) DelP(end)];
% =========================================================================