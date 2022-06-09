function [sys,x0,str,ts]= sfunvlvd4(t,x,ui,flag)

% Smart Fraction collection 2
%  Bijan Medi, SCBE, NTU, Oct. 2011.
% Random dtci for identification
% 2: SFC 2: Using slopes for finding the end of cycle
% 3: Min-SP values read again from simulink space
% 3_1: Modified FPY, LP1, LY1 removed, i index changed to k index
% 4: Inserting SS values instead of Min-SP values. Omitting FPY

% WARNING ====================================================
% Global variables from here are not updated in output function
% ============================================================

switch flag
    
    case 0
        [sys,x0,str,ts] = mdlInitializeSizes;
        
    case 1
        sys = [];
    case {2,9}
        sys = []; % do nothing
        
    case 3
        sys = mdlOutputs(t,x,ui);
        
    otherwise
        error(['unhandled flag = ',num2str(flag)]);
end

% end limintm


function [sys,x0,str,ts] = mdlInitializeSizes


% Solution Parameters
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 10;
sizes.NumInputs      = 20;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
str = [];


% =========================================================================
global check2 tsc2 ncF data1
global section temp12 tct12 sys12 Output
global slop i1
global x Nparams


check2 = 0; % Triggers data collection
tsc2 = 0;


data1 = zeros(1,4); % Stores feed and cut times data

section = 0; % Refer to notes
ncF = 0;
Output = zeros(1,15);
x = zeros(1,3);
Nparams = zeros(1,3);

sys12 = [];
tct12 = [];
temp12 = 0;
slop = [0 0];
i1 = 1;

% =========================================================================

% Initial values
x0 = [];

ts  = [0 1];   % sample time: [period, offset]

% end mdlInitializeSizes

function sys = mdlOutputs(t,y,ui)

global check2 ncF data1
global section temp12 tct12 sys12 Output Params
global x Nparams
global slop i1 ch1 ch2

% Steady State ==============================
PAs = ui(12); % purity of A at SS
PBs = ui(13); % purity of B at SS
YAs = ui(14); % recovery of A at SS
YBs = ui(15); % recovery of B at SS
Prs = ui(16); % g/min/g Productivity at SS
Drs = ui(17); % l/g Desorvent requirement at SS
% ==========================================

% Weight Factors --------------------------------------------

LPr2 = 1/0.02;
LDr2 = 1/0.2;
% -----------------------------------------------------------

Ptresh = 0.05; % Threshold of peak detection or end detection

yAout = ui(1); % normalized
yBout = ui(2); % normalized
nci = ui(3); % Injector Cycle Counter
% Error_cy = ui(4); % tinj>=tcy
Vinj = ui(5)/1e9; % ul -> m3
tcy = ui(6); % s
Qout = ui(7)/60/1e6; % ml/min->m3/s
cF = ui(8);
Ninj = ui(9);
cref = ui(10); % Reference conce. mg/ml
mad = ui(11); % mg

% DATA READ =======================
% Vinj, tcy, Qout, and cF are read from column s-function, which updates the data only at the
% end of each nci th cycle. dtci are stored at nci th location and so they
% are implemented when nci th cycle arrives.
data1(nci,1) = Vinj; % ul
data1(nci,2) = tcy; % s
data1(nci,3) = Qout; % m3/s Sent out directly along with cA and cB, but the stored value is used for DrF
data1(nci,4) = cF; % g/l
% =================================

Params.mad = mad;

Params.cref = cref;

% SMART FRACTION COLLECTION ==============================================

% The definition of fractions here is different than original fraction
% collection defined within the optimization problem
%    Section 0 /   Section 1    \   Section 2   / Section 3 \  Section 4/
% -------------------------------------------------------------------------
%  Ini. Solv.|First valley|First Peak|Second valley|Second Peak|First valley  
% -------------------------------------------------------------------------

% Note: Optimization is carried out at the end of each cycle not every time step.
if (yAout+yBout)*cref>=Ptresh && section == 0
    section = 1;
    
elseif slop(2)==1 && section ==1
    section = 2;
    check2 = 0; % Reset for next data collection
    
elseif slop(1)==1 && section ==2
    
    section = 3;
    
elseif slop(2)==1 && section ==3
    
    section = 4;
    
elseif slop(1)==1 && section ==4
    
    section = 1;
    check2 = 1; % Resetting check2 for new cycle
    
    % Updating ncF --------------------------------
    ncF = min(ncF + 1,Ninj); % Limits ncF to Ninj
    % ---------------------------------------------
    
    % DATA WRITE ==========================================
    % Data is written to vars from the respective cycle
    
    VinjF = data1(ncF,1); % m3 With reference to ncF th cycle
    %     tcyF = data1(ncF,2); % With reference to ncF th cycle
    %     QDF = data1(ncF,3); % m3/s With reference to ncF th cycle
    cFF = data1(ncF,4); % With reference to ncF th cycle
    % This group must be always together
    % =====================================================
    
    % Effective cycle time --------------------------
    tcyF = tct12(end)-tct12(1);
    % -----------------------------------------------
    
    Nparams = [tcyF tcyF tcyF];
    
    Params.VinjF = VinjF;
    Params.cFF = cFF;
    Params.Nparams = Nparams;
    
    % Elution Profiles -----------------------------
    Params.Data = sys12;
    Params.tct = tct12;
    % ----------------------------------------------
    
    % Reading and normalizing cut intervals----------------------------
    x1 = ui(18:20); % dtc1-3 as inputs
    x = x1'/sum(x1)*(tcyF-0.2)./Nparams; % Normalizing wrt tcyF and wrt Nparams
    % -----------------------------------------------------------------
    
    Output = funsfc2(x,Params);
    
    disp(['tcyF = ', num2str(tcyF)])

    
    sys12 = [];
    tct12 = [];
    temp12 = 0;
    i1 = 1;
    ch1 = 0;
    ch2 = 0;
    
end
% =========================================================================

switch section
    
    case 0
        % Do nothing, initial elution time
    case {1,2,3,4}
        % Finding Major Step ----------------------------------------------
        % This is intended for picking only major steps, though it is not
        % accurate as there are some minor steps forward.
        
        if isempty(tct12)
            tct12(1,1) = t;
            sys12 = [sys12;[cref*yAout,cref*yBout,Qout]];
        elseif t > tct12(end,1)
            
            temp12 = temp12+t;
            sys12 = [sys12;[cref*yAout,cref*yBout,Qout]];
            
            tct12 = [tct12;t];
            i1 = i1 + 1;
            
            % =======================================================
            %PEAK DETECTION
            % =======================================================
            % Based on testing several points in a series, the decision is
                % made if a peak has arrived or a valley. The last valley is
                % the end of one cycle. This algoritm fails if Vinj is too large.
            
            if i1>10 && rem(i1,1)==0
                % The test is done for every 1 samples
                
                ch1 = 0;
                ch2 = 0;
                
                y1 = sys12(end-8:end,1) + sys12(end-8:end,2);
                %                 ym = mean(y);
                %                 y1 = (y - ym)/ym; % Time consuming
                
                for k=1:5
                    
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
                    slop = [1 0];
                end
                
                if ch2>1
                    slop = [0 1];
                end
                
                
            end
        end
end
% -----------------------------------------------------------------

PA = Output(1);
PB = Output(2);
YA = Output(3);
YB = Output(4);
Pr = Output(5);
Dr = Output(6);

FPD = LPr2*(Pr-Prs)-LDr2*(Dr-Drs);
% FPY = LP2*((PA-PAmin)/100 + (PB-PBmin)/100) + LY2*((YA-YAmin)/100 + (YB-YBmin)/100)+...
%     0.5*(YA-YAmin)*(YB-YBmin)/1e4;

% FPY = LP2*((PA-PAmin)/100 + (PB-PBmin)/100) + LY2*((YA-YAmin)/100 + (YB-YBmin)/100)+...
%     0.5*(PA-PAmin)*(PB-PBmin)*(YA-YAmin)*(YB-YBmin)/50;


% FPY = LP2*((PA-PAmin)/100 + (PB-PBmin)/100) + LY2*((YA-YAmin)/100 + (YB-YBmin)/100)+...
%     10*max(0,YB-YBmin)/100;

sys(1)= (PA-PAs)*max(0,ncF)/(eps+ncF);  % PAb
sys(2)= (PB-PBs)*max(0,ncF)/(eps+ncF);  % PBb
sys(3)= (YA-YAs)*max(0,ncF)/(eps+ncF);  % YAb
sys(4)= (YB-YBs)*max(0,ncF)/(eps+ncF);  % YBb

sys(5)= Output(5); % PrF
sys(6)= Output(6); % DrF

sys(7) = FPD*max(0,ncF)/(eps+ncF); %
sys(8) = 0; %

sys(9) = check2; % Allows to collect data before reset
sys(10) = ncF; %
