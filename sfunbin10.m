function [sys,x0,str,ts]= sfunbin10(t,x,ui,flag,up)

% Bijan Medi, NTU, SCBE, 2010.
% sfunbin8: Used in conjunction with cfundy1 to rectify problems with global variables in C-sfunctions.
% QF removed. Nz is global. Interrupting program by Error_cy=1 directly from here.
% 9: Adding mad and Pr
% 10: Using cfundy3

% WARNING ====================================================
% Global variable from here are not updated in output function
% ============================================================

switch flag
    
    case 0
        [sys,x0,str,ts] = mdlInitializeSizes(ui);
        % clear all
        
    case 1
        sys = mdlDerivatives(t,x,ui,up);
        
    case {2,9}
        sys = []; % do nothing
        
    case 3
        sys = mdlOutputs(t,x,ui,up);
        
    otherwise
        error(['unhandled flag = ',num2str(flag)]);
end

% end limintm


function [sys,x0,str,ts] = mdlInitializeSizes(ui)

global Nz
% Solution Parameters
Nz = 400; % Number of grid points

sizes = simsizes;
sizes.NumContStates  = 4*Nz;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 13;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
% -----------------------------------

sys = simsizes(sizes);
str = [];
%***********************************

% Initial values
x0 = zeros(1,4*Nz);

ts  = [0 0];   % sample time: [period, offset]

% end mdlInitializeSizes

function dy = mdlDerivatives(t,y,ui,up)
% DERIVATIVES

global Nz L D eb Q cref HA KA HB KB rho Mu dp KovA KovB ycAb1 ycBb1

% C -----------------------------------------------------------------------
dy = cfundy3(t,y,[Nz L D eb Q cref 0 HA KA HB KB rho Mu dp KovA KovB ycAb1 ycBb1]);
1;
% -------------------------------------------------------------------------

% =========================================================================

function sys = mdlOutputs(t,y,ui,up)

global nc Nz L D eb Q cref HA KA HB KB rho Mu dp KovA KovB ycAb1 ycBb1
global tsc1 tcy Vinj Ninj cAin cBin cF

% OUTPUT ================================
% Global variables are globally pased to mdlderivatives. This necessiates
% that first, mdloutput is called and updated. Apparently it is true.

cref = ui(6); %

% Column Parameters
%Nz = up(0); % Number of gridpoints
L = up(2)*0.01; % m Enter in cm
D = up(3)*0.01; % m Enter in cm
eb = up(4); % Bed void fraction

% Isotherms ------------------------
HA = up(5); % 3.49 Henry constant component A
KA = up(6); % l/g 0.0550 Langmuir equilibrium constant component A

HB = up(7); % 1.41 Henry constant component B
KB = up(8); % l/g 0.0135 Langmuir equilibrium constant component B
%-----------------------------------------------------------

% Pressure-related parameters
rho = up(9); % Kg/m3 Liquid density
rhos = up(10); % Kg/m3 Stationary Phase density (non-porous) based on DAICEL data (600 Kg/m3 for porous)
Mu = up(11);% Pa.s
dp = up(12)*1e-6;% um -> m

KovA = up(13); % m/s
KovB = up(14); % m/s

% Adsorbent Properties -------------------------------------
mad = pi/4*(D^2)*L*(1-eb)*rhos*1e6; % mg of adsorbent inside the column (L and D in m)
% ----------------------------------------------------------

% Warning: Parameters below should not be indepently copied to Derivatives function
% INITIALZATION ========================================================
if (t==0)
    Vinj = ui(1)*1e-9; % Volume of injection loop (micL) --> m3
    tcy = ui(2); % Cycle time (sec). Initial value
    Q = ui(3)*1e-6/60; % Enter in mL/min! -> m3/s
    cF = ui(4); % Total feed conc. (mg/ml)
    Ninj = ui(5); % Number of injections
    tsc1 = 0; % Local cycle time s (Initialized with each cycle)
    nc = 1; % Cycle counter
    cAin = cF/2;%
    cBin = cF/2;%
end
% ======================================================================

if (t - tsc1 >= tcy)
    % UPDATE =============================================================
    Vinj = ui(1)*1e-9; % Volume of injection loop (micL) --> m3
    tcy = ui(2); % Cycle time (sec). Initial value
    Q = ui(3)*1e-6/60; % Enter in mL/min! -> m3/s
    cF = ui(4); % Total feed conc. (mg/ml)
    Ninj = ui(5); % Number of injections
    tsc1 = t;
    nc = nc + 1; % Cycle counter
    cAin = cF/2;%
    cBin = cF/2;%
    % ====================================================================
    
    ycAb1 = 0*cAin/cref; % To match the else condition
    ycBb1 = 0*cBin/cref; % To match the else condition
    %     QF = 0; % To match the else condition
    % tinj must be smaller than tcy
end
% ---------------------------------------------
% Should be updated before being used
tinj = Vinj/Q; % Injection interval sec
% ---------------------------------------------

if (tinj>=tcy)
    
    Error_cy = 1; % tinj>=tcy has happened, stop the program
    error('Error: Cycle time is smaller than injection time.')
    
    % Pulse Generator
elseif ((t-tsc1<=tinj) && (nc<=Ninj))
    ycAb1 = cAin/cref;
    ycBb1 = cBin/cref;
    %     QF = Q;
    Error_cy = 0;
else
    ycAb1 = 0*cAin/cref;
    ycBb1 = 0*cBin/cref;
    Error_cy = 0;
    %     QF = 0;
end
% ----------------------------------

% Pressure-related parameters
u0 = Q/(3.1415926*D*D/4); % ^ Superficial
dp7 = 500*u0*L*Mu/(dp*dp)/1e5; % ^ (bar) Darcy equation from Cox's book.
%  phi=500, resistance  parameter
%---------------------------

% Axial Dispersion and Mass Transfer -----------------------------
% Should be after Q

out(1) = y(Nz) + y(3*Nz); % uv
out(2) = y(Nz); %  ycA
out(3) = y(3*Nz); %  ycB
out(4) = dp7; % Pressure drop bar
out(5) = min(nc, Ninj); %  Limits the output to the requested number of injections
out(6) = Error_cy; % tinj>=tcy

out1(1) = Vinj*1e9; % ul
out1(2) = tcy;
out1(3) = Q*1e6*60; % ml/min
out1(4) = cF;
out1(5) = Ninj;
out1(6) = cref;
out1(7) = mad;

sys = [out,out1];

%========================================

% end mdlOutputs