function dy = funscc3(t,y,Params)

% Single Component Equilibrium Dispersive Model + Film Mass Transfer
% Bijan Medi, NTU, SCBE, 2010.
% fun 1: Based on WENO
% funscc1: Removing mad, rhos
% 3: All SI units. Changing t<=tinj to t<tinj, using cfundy3.c (corrected B.C.)

Nz = Params(1);
L = Params(2); %  m
D = Params(3); %  m
eb = Params(4);
Q = Params(5); % m3/s
cref= Params(6);
% input_type = Params(7);
HA = Params(8);
KA = Params(9);
HB = Params(10);
KB = Params(11);
cAin = Params(12);
cBin = Params(13);
Vinj = Params(14);% m3

tcy = Params(15);

rho = Params(16); % Kg/m3 Liquid density
Mu = Params(17);% Pa.s
dp = Params(18);% 

KovA = Params(19); % m/s
KovB = Params(20); % m/s

% ---------------------------------------------
% Should be updated before being used
tinj = Vinj/Q; % Injection interval sec
% ---------------------------------------------

if (tinj>=tcy)
    
%     Error_cy = 1; % tinj>=tcy has happened, stop the program
    error('Error: Cycle time is smaller than injection time.')
    
    % Pulse Generator
elseif t<tinj
    ycAb1 = cAin/cref;
    ycBb1 = cBin/cref;
else
    ycAb1 = 0*cAin/cref;
    ycBb1 = 0*cBin/cref;
end
% ----------------------------------

% DERIVATIVES =============================================================
% All units are in SI
dy = cfundy3(t,y,[Nz L D eb Q cref 0 HA KA HB KB rho Mu dp KovA KovB ycAb1 ycBb1]);
% =========================================================================
