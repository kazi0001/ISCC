% Smart Fraction collection 2
% Bijan Medi, SCBE, NTU, Oct. 2011.
% Peakdetd2: Random dtc1-3 and Vinj, QD. I/O relations differ from ARX dynamics
% d2: SFC 2: Using slopes for finding the end of cycle
% 2_1: Input reading Vinj(i) corrected to Vinji(ncF)
% dy1: Dynamic mode with random dtc1-3, Vinji(ncF) changed back to Vinji
% and Ouput shifted forward with an initial condition. SpointsF
% artificially added. dcti Plots corrected for Inputs
% 5u1: Five inputs Vinj, QD, dtc1-3
% 2: Modified FPY, LP1, LY1 removed
% 5u5o1: Five inputs Vinj, QD, dtc1-3. Five Outputs: FPD, PA, PB, YA, YB
% 2:
% 3: cFF time varying, imported with other vars
%

clc
clear
clear global
close all

% -----------------------------------------------------------
% Instruction
% -----------------------------------------------------------
% Change tcy
% Change P/Ys
% Change Prs, Drs
% Change LPr2, LDr2
% -----------------------------------------------------------

% load('uvqsdata15.mat') %

% % Fileid = 'lp2rndc1'; % New SS, random Vinj,QD ,dtci [300,1,0.3,2] with cF=34.91 imported
% Fileid = 'lp2test1'; % 
% Fileid = 'lp2crg1'; % RG Change in His from 3.49 to 3.32 and 1.41 to 1.34
% Fileid = 'lp2csp1'; % SP from [0.6 95 95 90 92] to [0.6 98 98 95.0 95.0] at t = 3000
Fileid = 'lp2crg2'; % RG Change in cF from 34.91 to 34.0 and at nci = 38, t=1100

load([Fileid,'.mat']) %

t = UVout.time;
Nt = length(t);
Output = zeros(1,15);

yAout = UVout.signals(1,1).values(:,2);
yBout = UVout.signals(1,1).values(:,3);
Qout = QD/60/1e6; % ml/min->m3/s
Vinj = Vinj/1e9; % uL->m3
% ----------------------------------------------------------

clear UVout

% Weight Factors --------------------------------------------
LPr2 = 1/0.02;
LDr2 = 1/0.2;
% -----------------------------------------------------------


Ptresh = 0.05; % Threshold of peak detection

D = 1*0.01; % cm->m Column diameter
eb = 0.704; % Bed void fraction
tres = 0.01; % Resolution of output time array. This regulates the solver
L = 10*0.01; % cm->m Total length of unit

% -------------------------------
tcy = 30; % s Cycle time
% -------------------------------

rho = 785.8; % Kg/m3 Liquid density. Heptane-Ethanol (65/35 v/v) from Perry at 23
rhos = 2027.03; % Kg/m3 for nonporous solid based on DAICEL data (600 Kg/m3 for porous)
mad = 1e6*pi/4*(D^2)*L*(1-eb)*rhos; % mg of total adsorbent inside the unit (L and D in m)
% -------------------------

% % Isotherm --------------------------------------
% HA = 3.49; % Henry constant component A
% KA = 0.0550; % l/g
%
% HB = 1.41; % Henry constant component B
% KB = 0.0135;% l/g
%
% % -----------------------------------------------

% Steady State ==============================
PAs = 98; % purity of A at SS
PBs = 98; % purity of B at SS
YAs = 95.5; % recovery of A at SS
YBs = 95.1; % recovery of B at SS
Prs = 0.0254; % g/min/g Productivity at SS
Drs = 0.174; % l/g Desorvent requirement at SS
% ==========================================

cref = 1; % g/l

Params.L = L;
Params.D = D;
Params.eb = eb;
Params.tres = tres;

Params.rho = rho;
Params.rhos = rhos;
Params.mad = mad;

Params.cref = cref;
% Params.LPr1 = LPr1;

section = 0; % Refer to notes
nci = 1; % Number of cycles injected
ncF = 1; % Number of cycles collected
Vinji = Vinj(1); % m3 Injection volume based on injection (sampled per cycle time)
cFi = cFpoints(1); % g/L feed conc. based on injection (sampled per cycle time)

sys12 = [];
tct12 = [];
temp12 = 0;
slop = [0 0];
i1 = 1;

SpointsF = zeros(1000,5); % Setpoint values at the time of calculting cycle values
DpointsF = zeros(1000,3); % Raw dtc1-3 values

for i=1:Nt
    
    % The definition of fractions here is different than original fraction
    % collection defined within the optimization problem
    %    Section 0 /   Section 1    \   Section 2   / Section 3 \  Section 4/
    % -------------------------------------------------------------------------
    %  Ini. Solv.|First valley|First Peak|Second valley|Second Peak|First valley
    % -------------------------------------------------------------------------
    
    if (yAout(i) + yBout(i))*cref>=Ptresh && section == 0
        section = 1;
        
    elseif slop(2)==1 && section ==1
        section = 2;
        
    elseif slop(1)==1 && section ==2
        
        section = 3;
        
    elseif slop(2)==1 && section ==3
        
        section = 4;
        
    elseif slop(1)==1 && section ==4
        
        section = 1;
        
        
        
        
        % Effective cycle time --------------------------
        tcyF = tct12(end)-tct12(1);
        % -----------------------------------------------
        
        
        Nparams = [tcyF tcyF tcyF];
        
        SpointsF(ncF,:) = Spoints(i,:); % Setpoint values at the end of each cycle
        DpointsF(ncF,:) = Dpoints(i,:); % Raw data for dtc1-3
        
        Params.VinjF = Vinji(ncF); % m3
        
        % ============================================================
        cFF = cFi(ncF); % Feed concentration at ncF th injected cycle
        % ============================================================
        
        Params.cFF = cFF;
        
        % Elution Profiles -----------------------------
        Params.Data = sys12;
        Params.tct = tct12;
        % ----------------------------------------------
        
        Params.Nparams = Nparams;
        
        
        % Preparing Cut Intervals -----------------------------------------
        x1 = DpointsF(ncF,1:3); % dtc1-3 as inputs
        x(ncF,:) = x1/sum(x1)*(tcyF-0.2)./Nparams; % Normalizing wrt tcyF and wrt Nparams
        % -----------------------------------------------------------------
        %         x(ncF,:).*Nparams;
        
        % Assigning Outputs -------------------------------------
        % Outputs are shifted forward to avoid direct feedthrough
        Output(ncF+1,:) = funsfc2(x(ncF,:),Params);
        % -------------------------------------------------------
        
        % Inputs ======================================
        Input(ncF,1) = Vinj(i)*1e9; % uL, Based on what measured  at the current moement
        Input(ncF,2) = Qout(i)*60*1e6; % mL/min, Based on what measured  at the current moement
        Input(ncF,3:5) = x(ncF,:).*Nparams; % dtc1-3 for ncF cycle
        % =============================================
        
        sys12 = [];
        tct12 = [];
        temp12 = 0;
        i1 = 1;
        
        disp(['tcyF = ', num2str(tcyF)])
        ncF = ncF + 1;
        
    end
    
    switch section
        
        case 0
            
        case {1,2,3,4}
            % Finding Major Step ----------------------------------------------
            % This is intended for picking only major steps, though it is not
            % accurate as there are some minor steps forward.
            % The time series are variable step. So, a variable step integrator
            % is required.
            
            if isempty(tct12)
                tct12(1,1) = t(i);
                sys12 = [sys12;[cref*yAout(i),cref*yBout(i),Qout(i)]];
            elseif t(i) > tct12(end,1)
                
                temp12 = temp12+t(i);
                sys12 = [sys12;[cref*yAout(i),cref*yBout(i),Qout(i)]];
                tct12 = [tct12;t(i)];
                i1 = i1+1; % Time index counter for one cycle
                
                % =======================================================
                %PEAK DETECTION
                % =======================================================
                % Based on testing several points in a series, the decision is
                % made if a peak has arrived or a valley. The last valley is
                % the end of one cycle. This algoritm fails if Vinj is too large.
                
                if i1>6 && rem(i1,1)==0
                    % The test is done for every 1 samples
                    
                    ch1 = 0;
                    ch2 = 0;
                    
                    y1 = sys12(end-5:end,1) + sys12(end-5:end,2);
                    %                 ym = mean(y);
                    %                 y1 = (y - ym)/ym; % Time consuming
                    
                    for k=1:2
                        
                        % Valley detection
                        if y1(end-k)>y1(end-k-1) && y1(end-k-1)>y1(end-k-2) && y1(end-k-2)>y1(end-k-3)
                            ch1 = ch1 + 1;
                            
                        end
                        
                        % Peak detection
                        if y1(end-k)<y1(end-k-1) && y1(end-k-1)<y1(end-k-2) && y1(end-k-2)<y1(end-k-3)
                            ch2 = ch2 + 1;
                        end
                    end
                    
                    if ch1>1
                        slop = [1 0]; % Valley has passed
                    end
                    
                    if ch2>1
                        slop = [0 1]; % Peak has passed
                    end
                    
                    
                end
            end
    end
    
    % ==============================================
    % Saving Vinj, cF values from injector point of view
    if t(i)>=nci*tcy
        Vinji(nci+1,1) = Vinj(i);
        cFi(nci+1,1) = cFpoints(i);
        nci = nci + 1;
    end
    % ==============================================
    
end
% ncF = ncF-1;
% nci = nci-1;

% I/O Data Processing -----------------------------------------------------
% Calculating initial values backward
% Truncating samples; the last y sample becomes useless


% y0|y1|...|yncF-1
% ------------------------
% u1|u2|...|uncF
% 1 |2 |...|ncF
Output(1,:) = Output(2,:); % Dummy initial condition
Output = Output(1:ncF-1,:); % Data has already been shifted one sample forward
% -------------------------------------------------------------------------

PA = Output(:,1);
PB = Output(:,2);
YA = Output(:,3);
YB = Output(:,4);
Pr = Output(:,5);
Dr = Output(:,6);

FPD = LPr2*(Pr-Prs)-LDr2*(Dr-Drs);

% Mains -------------------------------------------------------------------
subplot(2,1,1),stairs([0:ncF-2],Output(:,1),'r-','LineWidth',1.5) % PA
hold on
subplot(2,1,1),stairs([0:ncF-2],Output(:,2),'b-','LineWidth',1) % PB
subplot(2,1,1), stairs([0:ncF-2],SpointsF(1:ncF-1,2),'r--','LineWidth',1) % PAsp
subplot(2,1,1), stairs([0:ncF-2],SpointsF(1:ncF-1,3),'b-.','LineWidth',1) % PBsp
ylabel('Pi (%)'), grid on


subplot(2,1,2), stairs([0:ncF-2],Output(:,3),'r-','LineWidth',1.5) % YA
hold on
subplot(2,1,2), stairs([0:ncF-2],Output(:,4),'b-','LineWidth',1) % YB
subplot(2,1,2), stairs([0:ncF-2],SpointsF(1:ncF-1,4),'r--','LineWidth',1) % YAsp
subplot(2,1,2), stairs([0:ncF-2],SpointsF(1:ncF-1,5),'b-.','LineWidth',1) % YBsp
ylabel('Yi (%)'),grid on
xlabel('Cycles')
h = get(0,'CurrentFigure');  saveas(h,[Fileid,'_Limits'],'emf');
% -------------------------------------------------------------------------

figure
% Cut Intervals and Cycle time --------------------------------------------
subplot(5,1,1),stairs([0:ncF-2],Output(:,7),'r-','LineWidth',1.5) % tcyF
ylabel('tcy (s)'), grid on

subplot(5,1,2),stairs([1:ncF-1],Input(:,3),'b-','LineWidth',1.5) % dtc1
ylabel('dtc1 (s)'), grid on

subplot(5,1,3), stairs([1:ncF-1],Input(:,4),'b-','LineWidth',1.5) % dtc2
ylabel('dtc2 (s)'),grid on

subplot(5,1,4), stairs([1:ncF-1],Input(:,5),'b-','LineWidth',1.5) % dtc3
ylabel('dtc3 (s)'), grid on

subplot(5,1,5), stairs([0:ncF-2],Output(:,11),'b-','LineWidth',1.5) % dtc4
xlabel('Cycles'),ylabel('dtc4 (s)'), grid on
h = get(0,'CurrentFigure');  saveas(h,[Fileid,'_Cuts'],'emf');
% -------------------------------------------------------------------------

figure
% Outputs -----------------------------------------------------------------
subplot(3,1,1),stairs([0:ncF-2],FPD,'r-','LineWidth',1.5) % FPD
ylabel('FPD'), grid on

subplot(3,1,2), stairs([0:ncF-2],Output(:,5),'r-','LineWidth',1.5) % Pr
hold on
ylabel('Pr (g/(min g))'), grid on

subplot(3,1,3), stairs([0:ncF-2],Output(:,6),'b-','LineWidth',1.5) % Dr
ylabel('Dr (L/g)'), grid on
xlabel('Cycles')
h = get(0,'CurrentFigure');  saveas(h,[Fileid,'_Mains'],'emf');
% -------------------------------------------------------------------------

figure
% Inputs ------------------------------------------------------------------
subplot(2,1,1),stairs([1:ncF-1],Input(:,1),'r-','LineWidth',1.5)
ylabel('Vinj (uL)'), grid on

subplot(2,1,2),stairs([1:ncF-1],Input(:,2),'b-','LineWidth',1.5)
xlabel('Cycles'),ylabel('QD (mL/min)'), grid on
h = get(0,'CurrentFigure');  saveas(h,[Fileid,'_Inputs'],'emf');
% -------------------------------------------------------------------------

% Data Processing for Plotting --------------------------------------------
Inplot = [[1:ncF-1]',Input];
Outplot = [[0:ncF-2]',Output(:,1:11),FPD,SpointsF(1:ncF-1,:)];
% -------------------------------------------------------------------------


% clear Output Input Qout PA PB YA YB yAout yBout
save([Fileid,'_IOdata'],'Inplot','Outplot')
% clear Inplot OutplotD
