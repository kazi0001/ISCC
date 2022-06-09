% Single-Column Optimization Studies
% Bijan Medi & Kazi Monzure Khoda, NTU, SCBE, 2011.


% From gamulsccib1
% gamulsmb1: The main code is from odesmbdskv5
% gamulsccd3: Using strucutures, increasing Qmax from 30 to
% 70, decreasing tcymin from 15 to 10, decreasing Vinjmin from 500 to 300 using parallel comp and parfor

% a1:
% a2:

clc
clear
clear global

% clear all

close all


% OUTPUT FILE IDENTIFIER  ----------------------------------
FileID = 'OSCA4'; % Adds an indentifier to the output files
% ----------------------------------------------------------

% --------------------------------------------------------------------------
Reproduce = 0; % New Problem = 0, Cont. from an incomplete run = 1, Cont. from last finished run = 2
% -------------------------------------------------------------------------

switch Reproduce
    
    case 0 % New Run
        
        % Weight Factors --------------------------------------------
        Lpt = 100; % Penalty factor
        LambdaP = []; % Productivity, reserved for GA
        LambdaD = []; % Desorbent, reserved for GA
        TolCon1 = 0.5; % Constraint Tolerance for PA/B and YA/B (%)
        TolCon2 = 1; % Constraint Tolerance for Pressure Drop (bar)
        % -----------------------------------------------------------
        
        % Constraints ==============================
        PAmin = 90; % Minimum purity of A
        PBmin = 90; % Minimum purity of B
        YAmin = 85; % Minimum recovery of A
        YBmin = 85; % Minimum recovery of B
        DelPmax = 40; % Maximum pressure drop bar
        % ==========================================
        
        Params.Lpt = Lpt;
        Params.TolCon1 = TolCon1;
        Params.TolCon2 = TolCon2;
        
        Params.PAmin = PAmin;
        Params.PBmin = PBmin;
        Params.YAmin = YAmin;
        Params.YBmin = YBmin;
        
        Params.DelPmax = DelPmax;
        
        
        
        % Normalizing Values -------------------------------
        Nparams = [700 68 9 26 20 1.1 45];
        % --------------------------------------------------
        Params.Nparams = Nparams;
        
        
        D = 1*0.01; % cm->m Column diameter
        eb = 0.704; % Bed void fraction
        Nz = 400; % Number of grid points for each column
        tres = 0.01; % Resolution of output time array. This regulates the solver
        tend = []; % s simulation time
        tbias = -0.001; % bias time set for switching
        ncol = [1]; % Array of columns in four sections
        ntot = sum(ncol); % Total number of columns
        L = 10*0.01; % cm->m Total length of unit
        Ncy = 3; % Minimum Number of cycles
        
        rho = 785.8; % Kg/m3 Liquid density. Heptane-Ethanol (65/35 v/v) from Perry at 23
        rhos = 2027.03; % Kg/m3 for nonporous solid based on DAICEL data (600 Kg/m3 for porous)
        mad = 1e6*pi/4*(D^2)*L*(1-eb)*rhos; % mg of total adsorbent inside the unit (L and D in m)
        Mu = 6.075e-4; % Pa.s % Heptane-Ethanol (65/35 v/v) from DIPPR at 23
        dp = 20*1e-6; % um->m Particle diameter
        phi=500; % resistance  parameter for Darcy equation
        KovA = 0.0061*0.01; % m/s from Zabka 2008;
        KovB = 0.0061*0.01; % m/s from Zabka 2008;
        % -------------------------
        Vcol = (pi/4)*D^2*(L/ntot); % m3 Volume of one column
        
        
        % Isotherm --------------------------------------
        HA = 3.49; % Henry constant component A
        KA = 0.0550; % l/g
        
        HB = 1.41; % Henry constant component B
        KB = 0.0135;% l/g
        
        Nis = []; % Max 36 points
        cAarr = [];
        cBarr = [];
        qAarr = [];
        qBarr = [];
        qAmat = [];
        qBmat = [];
        % -----------------------------------------------
        
        cref = 1; % g/l
        
        
        % =========================================================================
        odeoptions = odeset('Initialstep',0.001,'RelTol',1e-4,'AbsTol',1e-6);
        Params.odeoptions = odeoptions;
        
        % Lower bounds ------------------------------------------------------------
        % Normalized values
        %     Vinj     tcy      QD    cF   dtc1     dtc2   dtc3
        lb = [300      10        5    10    1        0.2    1 ]./Nparams;
        % -------------------------------------------------------------------------
        
        % Upper bounds ------------------------------------------------------------
        % Normalized values
        %     Vinj    tcy     QD   cF  dtc1    dtc2    dtc3
        ub = [5000    90      70   35   90      90      90]./Nparams;
        % -------------------------------------------------------------------------
        
        % Inequalities ------------------------------------------------------------
        %         Vinj   tcy    QD   cF    dtc1    dtc2      dtc3
        A1 =   -1*[ 0      1     0    0      -1      -1       -1].*Nparams;
        A2 = [      0      0     1    0       0       0        0].*Nparams;
        A = [A1; A2];
        b1 = -0.2;
        b2 = [DelPmax/(phi/(pi*D^2/4)*L*Mu/dp^2/1e6/60/1e5)]; % (L and D in m)
        b = [b1 b2];
        % -------------------------------------------------------------------------
        
        
        load('OSCA1output.mat');
        
        x01 = Popend(3:end,:);
        
        x02 = [2600,19,34.0,34.94,4.3,0.5,13.8]./Nparams;
        x03 = [2574,20.26,28.39,34.91,5.,0.45,14.27]./Nparams;
        
        x0 = [x01;x02;x03];
        
        s0 = [];
        
        % -----------------------------------------------
        
        Params.L = L;
        Params.D = D;
        Params.eb = eb;
        Params.Nz = Nz;
        Params.tres = tres;
        Params.tend = tend;
        Params.tbias = tbias;
        Params.ncol = ncol;
        Params.ntot = ntot;
        
        Params.rho = rho;
        Params.rhos = rhos;
        Params.mad = mad;
        Params.Mu = Mu;
        Params.dp = dp;
        Params.phi = phi;
        Params.KovA = KovA;
        Params.KovB = KovB;
        
        Params.Vcol = Vcol;
        Params.Ncy = Ncy;
        
        % Isotherm --------------------------------------
        Params.HA = HA;
        Params.KA = KA;
        Params.HB = HB;
        Params.KB = KB;
        % -----------------------------------------------
        
        Params.cref = cref;
        % -------------------------------------------
        
        
        % Setting Random Stream ------------------------
        % For a new run
        savedStream=RandStream.getDefaultStream;
        savedState=savedStream.State;
        savedType=savedStream.Type;
        rngstate.state = savedState;
        rngstate.type = savedType;
        % -----------------------------------------------
        
        % OPTIONS =================================================================
        options = gaoptimset(@gamultiobj);
        options = gaoptimset(options,'Generations',50);
        options = gaoptimset(options,'PopulationSize',80);
        options = gaoptimset(options,'InitialPopulation' , x0);
        options = gaoptimset(options,'CrossoverFcn',{@distancecrowding,'genotype'});
        
        % Cross over function ---------------------------------------------
        options = gaoptimset(options,'InitialScores' , s0);
        % -----------------------------------------------------------------
        
        options = gaoptimset(options,'CrossoverFcn',{@crossoverintermediate,1.1});
        options = gaoptimset(options,'CrossoverFraction', 0.3);
        options = gaoptimset(options,'CreationFcn',@gacreationlinearfeasible); % Both for ga and gamul
        options = gaoptimset(options,'PopInitRange' ,[lb;ub]);
        options = gaoptimset(options,'MutationFcn' ,@mutationadaptfeasible); % Both for ga and gamul
        options = gaoptimset(options,'ParetoFraction', 0.5);
        options = gaoptimset(options,'Display' ,'iter');
        options = gaoptimset(options,'PlotFcns' ,{@gaplotpareto,@gaplotrankhist,@gaplotscorediversity,@gaplotspread});
        options = gaoptimset(options,'StallTimeLimit', inf);
        options = gaoptimset(options,'OutputFcn',@fungamulout1_1);
        options = gaoptimset(options,'StallGenLimit', 5);
        options = gaoptimset(options,'TolFun', 1e-3);
        options = gaoptimset(options,'TolCon', 1e-3);
        options = gaoptimset(options,'UseParallel', 'always', 'Vectorized', 'off');
        % =========================================================================
        
        % Functions ----------------------------------
        f1 = @(x)funobjscco1(x,Params); % Change the handle in cases 1,2 accordingly
        %f2 = @(x)funganlc2(x,[Params Params_n]);
        % --------------------------------------------
        
        % Setting up problem ------------------
        problem.fitnessfcn = f1;
        problem.nvars = 7;
        problem.Aineq = A;
        problem.Bineq = b;
        problem.Aeq = [];
        problem.Beq = [];
        problem.lb = lb;
        problem.ub = ub;
        problem.rngstate = rngstate;
        problem.solver = 'gamultiobj';
        problem.options = options;
        % --------------------------------------
        
        RepID = 'N.A.';
        
    case 1 % Cont. from an incomplete run
        
        RepID = 'WDI4N'; % Last saved data
        
        % Initial Condition = Last population and scores -----
        load([RepID])
        load([RepID,'problem']);
        options = problem.options;
        x0 = Popend; %
        s0 = Scorend; %
        options = gaoptimset(options,'InitialPopulation' , x0);
        options = gaoptimset(options,'InitialScores' , s0);
        % ----------------------------------------------------
        
        % Setting Random Stream -------------------------------------------
        % Setting Random Stream from previous run
        set(RandStream.getDefaultStream,'State',problem.rngstate.state); % Comment for new randoms
        % -----------------------------------------------------------------
        
        % Change options here:
        % options = gaoptimset(options,'Generations',3);
        
        % Updating problem ----------------------------------
        problem.options = options;
        % ---------------------------------------------------
        
    case 2 % Continue from last finished run
        
        RepID = 'OSCA3'; % Last saved data
        
        %
        % OSCA1: Starting from 1 good point in QD and popend from OSCA1 with removing 1 points, Lpt = 100
        % Converged!
        % OSCA2: Starting from 1 good point in QD and popend from OSCA1 with removing 1 points, Lpt = 50
        % OSCA3: Starting from 2 good point in QD and popend from OSCA1 with removing 2 points, Lpt = 50
        % OSCA4: Starting from 2 good point in QD and popend from OSCA1 with removing 2 points, Lpt = 100
        % Converged!

        
        
        % Initial Condition = Last population and scores -----
        load([RepID,'output']);
        load([RepID,'problem']);
        options = problem.options;
        x0 = Popend;
        s0 = Scorend;
        
        
        
        
        
        options = gaoptimset(options,'InitialPopulation' , x0);
        options = gaoptimset(options,'InitialScores' , []);
        % ----------------------------------------------------
        
        % Setting Random Stream -------------------------------------------
        % Setting Random Stream from previous run
        set(RandStream.getDefaultStream,'State',problem.rngstate.state); % Comment for new randoms
        % ---------------------------------------------- -------------------
        
        
        % =================================================
        %         options = gaoptimset(options,'Generations',20);
        
        % =================================================
        %         %
        %         Weight Factors --------------------------------------------
        Lpt = 100; % Penalty factor
        %         % -----------------------------------------------------------
        %         %
        Params.Lpt = Lpt;
        f1 = @(x)funobjscco1(x,Params);
        problem.fitnessfcn = f1;
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        % Updating problem ----------------------------------
        problem.options = options;
        problem.rngstate = output.rngstate; % Necessary?
        % ---------------------------------------------------
end

% -------------------------------------------------------------------------

% Saving problem -------------------------------------
save([FileID,'problem'],'problem','Params', 'Nparams','options','lb','ub','PAmin','PBmin','YAmin','YBmin','DelPmax');
% ----------------------------------------------------

global outfileiters outfilevars opt1 opt2 bounds

outfileiters = 'Iters';%  Saving iteration details
outfilevars = 'Vars'; % Saving variables vs. iterations

outfileiters = [FileID,outfileiters,'.txt'];
outfilevars = [FileID,outfilevars,'.txt'];

opt1 = ['RepID = ', RepID, '  Max Generations = ', ...
    num2str(gaoptimget(options,'Generations')), '    PopulationSize  = ',...
    num2str(gaoptimget(options,'PopulationSize')), '    PopInitRange  = [',...
    num2str(min(gaoptimget(options,'PopInitRange'))),';',num2str(max(gaoptimget(options,'PopInitRange'))),']'];

bounds.L = lb;
bounds.U = ub;

% Saving lb, ub
opt2 = ['lb = [',num2str(lb,'%-3.6f  '),']','   ub = [',num2str(ub,'%-3.6f  '),']'];

tic
% SOLVER ==================================================================
[x,fval,exitflag,output,population,scores] = gamultiobj(problem);
% =========================================================================

trun = toc

x
fval
exitflag

% Saving exit details -----------------------------------------------------

% Saving the last data to file --------------------------------------
% Saved in mat format
Popend = population;
Scorend = scores;

save([FileID,'output'],'output','Popend','Scorend')

[PA,PB,YA,YB,Pr,Dr,DelP] = fundgenscco1(Popend,Scorend,problem.fitnessfcn,FileID);
% [PA,PB,YA,YB,PR,DR,DelP] = fundatagen1(Popend,Scorend,FileID,NewID); % For storing in a new file
Paretostruct = funrkfindscco2(FileID);


