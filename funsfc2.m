function Out = funsfc2(x,Params)

% funsfc1: Smart fraction collection scheme for change in flow rate
% 2: Using structure. All units in SI

% System Parameters -------------------------------------------------------

mad = Params.mad/1e6; % mg->Kg

Data = Params.Data;
tct = Params.tct;

tcyF = tct(end)-tct(1);

Nparams = Params.Nparams;


% % INPUTS ==================================================================
% % All for current collected cycle
VinjF = Params.VinjF; % m3

cFF = Params.cFF; % g/l

dtc1 = x(1)*Nparams(1); % sec
dtc2 = x(2)*Nparams(2); % sec
dtc3 = x(3)*Nparams(3); % sec
dtc4 = tcyF - (dtc1+dtc2+dtc3);

tsc = tct(1); % Start of cycle

[unsed,i1] = min(abs(tct-tsc-dtc1)); % Index of first cut time
[unsed,i2] = min(abs(tct-tsc-(dtc1+dtc2))); % Index of second cut time
[unsed,i3] = min(abs(tct-tsc-(dtc1+dtc2+dtc3))); % Index of third cut time

if (dtc1+dtc2+dtc3)>=tcyF
    error('Error: Cycle time is smaller than or equal to sum of three cut times.')
end

%         % m=int(Qout(t)*cA(t),t)
mA1 = trapz(tct(1:i1),Data(1:i1,3).*Data(1:i1,1)); % kg A collected in Fraction 1
mB1 = trapz(tct(1:i1),Data(1:i1,3).*Data(1:i1,2)); % kg B collected in Fraction 1

mA2 = trapz(tct(i1:i2),Data(i1:i2,3).*Data(i1:i2,1)); % kg A collected in Fraction 2
mB2 = trapz(tct(i1:i2),Data(i1:i2,3).*Data(i1:i2,2)); % kg B collected in Fraction 2

mA3 = trapz(tct(i2:i3),Data(i2:i3,3).*Data(i2:i3,1)); % kg A collected in Fraction 3
mB3 = trapz(tct(i2:i3),Data(i2:i3,3).*Data(i2:i3,2)); % kg B collected in Fraction 3

mA4 = trapz(tct(i3:end),Data(i3:end,3).*Data(i3:end,1)); % kg A collected in Fraction 4
mB4 = trapz(tct(i3:end),Data(i3:end,3).*Data(i3:end,2)); % kg B collected in Fraction 4

PB = mB1/(mA1 + mB1)*100;
PA = mA3/(mA3 + mB3)*100;

YB = mB1/(mB1 + mB2 + mB3 + mB4)*100;
YA = mA3/(mA1 + mA2 + mA3 + mA4)*100;

% PRODUCTVITY ===============================================
PrF = (mA1 + mB1 + mA3 + mB3)/tcyF/mad*60; % (kg/kg/min) mg F/mg Ad./min
% ===========================================================

% DESORBENT REQUIREMENT ==============================
DrF = (trapz(tct,Data(:,3)) + VinjF)/VinjF/cFF; % Desorbent requirement with reference to ncF th cycle mlD/mgF
% int(tct*Q is converted to m3/s
%=====================================================

Out = [PA PB YA YB PrF DrF tcyF dtc1 dtc2 dtc3 dtc4 mA1 mB1 mA3 mB3];