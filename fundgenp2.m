function [PAend, PBend, YAend, YBend, Prend, Drend, DelPend] = fundgenp2(i)

% global PAend PBend YAend YBend Prend Drend

% Loading Problem
% FileID = Params2.Fid;
% i = Params2.i;
load('gentemp','FileID');

load([FileID,'output']) %
load([FileID,'problem']) %

x = Popend(i,:);

out = funodescco3(x.*Params.Nparams,Params); % Receives non-normalized values

    PAend = out(1);
    PBend = out(2);
    YAend = out(3);
    YBend = out(4);
    Prend = out(5); 
    Drend = out(6);
    DelPend = out(7);