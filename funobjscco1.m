function Obj = funobjscco1(x,Params)
% Single-Column Optimization studies
% Bijan Medi & Kazi Monzure Khoda, SCBE, NTU, 2011.
% funobjsmb1: Formulating objective function, YA/B and DelP not included
% funobjscco1:

global  PAend PBend YAend YBend Drend Prend DelPend


Lpt = Params.Lpt;

% Contraints --------------------------------------------------------------
PAmin = Params.PAmin;
PBmin = Params.PBmin;
YAmin = Params.YAmin;
YBmin = Params.YBmin;
% DelPmax = Params.DelPmax;
% -------------------------------------------------------------------------


% Vectorization ---------------------
Nv = size(x,1); % No. of rows in x
% -----------------------------------


% Decision Variables ======================================================
j=1; % Reserved for vectorized
Vinj_v = x(j,1)*Params.Nparams(1); % Injection volume uL
tcy_v = x(j,2)*Params.Nparams(2); % Cycle time sec
QD_v = x(j,3)*Params.Nparams(3); % Solvent flow rate ml/min
cF_v = x(j,4)*Params.Nparams(4); % Total feed concentration mg/ml

dtc1_v = x(j,5)*Params.Nparams(5); % Cut time 1 sec
dtc2_v = x(j,6)*Params.Nparams(6); % Cut time 2 sec
dtc3_v = x(j,7)*Params.Nparams(7); % Cut time 3 sec


% tinj<tcy ---------------------------------------------
% Adding max value for Vinj = min(Vinj,0.99*QD*tcy)
if Vinj_v>=(1000/60)*0.99*QD_v*tcy_v
    Vinj_v = (1000/60)*0.99*QD_v*tcy_v; % ml/min -> ul/s Change in funganlc accordingly
warning('Injection time must be smaller than cycle time. Vinj is capped at  "%-4.3f" uL.', a)
end
% ------------------------------------------------------

% Running M-file-based SCC model ------------------------------------------
out = funodescco3([Vinj_v;tcy_v;QD_v;cF_v;dtc1_v;dtc2_v;dtc3_v],Params);
% -------------------------------------------------------------------------

PAend = out(1);
PBend = out(2);
YAend = out(3);
YBend = out(4);
Prend = out(5);
Drend = out(6);
DelPend = out(7);

% DelPmax not included

Obj(1)  = 1/(eps+Prend) + Lpt*(max(0,(PAmin-PAend)/PAmin) + max(0,(PBmin-PBend)/PBmin) + max(0,(YAmin-YAend)/YAmin) + max(0,(YBmin-YBend)/YBmin));
Obj(2)  = Drend + Lpt*(max(0,(PAmin-PAend)/PAmin) + max(0,(PBmin-PBend)/PBmin) + max(0,(YAmin-YAend)/YAmin) + max(0,(YBmin-YBend)/YBmin));