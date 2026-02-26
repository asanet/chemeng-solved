function orbitfinder_gui05
%% Description:
%
% Version: alpha-4
% Email: ataide@peq.coppe.ufrj.br
% Universidade Federal do Rio de Janeiro
% COPPE - Chemical Engineering Program
%

V = 'alpha-5';
cf = findobj('Name',['orbitfinder_' V]);
delete(cf);
start_orbitfinder(V);

end

function quitFunction(ho,~)
h = guidata(ho);
cf = findobj('Name',['orbitfinder_' h.Version]);
delete(cf);
end

%% Interface Construction
function start_orbitfinder(V)

%Sets the units of your root object (screen) to pixels
set(0,'units','pixels')

%Obtains this pixel information
Pix_SS = get(0,'screensize');
fs_num = ceil(Pix_SS(3)/192);

%Obtain OS
os = computer;

%Set font
if lower(os(1)) == 'g'
    myFont = 'Monospaced';
else
    myFont = 'fixedWidth';
end

%Aliases
fgc = 'ForegroundColor'; hal = 'HorizontalAlignment'; nzd = 'normalized'; pst = 'Position';
bgc = 'BackgroundColor'; str = 'String'; fa = 'FontAngle'; fw = 'FontWeight'; fs = 'FontSize';
pfgc2 = [.5 .5 .5]; pbc = [1 1 1];

%Create base figure
fh = figure('MenuBar','none','ToolBar','none','Units',nzd,'Resize','off', ...
    'NumberTitle','off',pst,[0.1 0.1 .8 .8],'Name',['orbitfinder_' V], ...
    'Tag','main','Color',[.95 .95 .95],'Visible','off','nextplot','new');

%Create figure menus
mh1 = uimenu('Parent',fh,'Label','File');
mh2 = uimenu('Parent',fh,'Label','Tools');
mh3 = uimenu('Parent',fh,'Label','Help');

uimenu(mh1,'Label','Save','Accelerator','S','Enable','off','Tag','msave');
uimenu(mh1,'Label','Load','Accelerator','L','Tag','mload');
uimenu(mh1,'Label','Reset','Accelerator','R','Tag','mreset');
uimenu(mh1,'Label','Open example file','Tag','mexample');
uimenu(mh1,'Label','Quit      ','Separator','on','Accelerator','Q','Tag','mquit');
uimenu(mh2,'Label','EigenMap','Separator','off','Accelerator','E','Tag','eigenmap');
uimenu(mh2,'Label','BrachDiagram','Separator','off','Accelerator','B','Tag','branchdiagram');
uimenu(mh3,'Label','About','Separator','off','Tag','mabout');
uimenu(mh3,'Label','Edit me','Separator','off','Tag','meditme');


%Create main panel
ph1 = uipanel('Parent',fh,'Units',nzd,pst,[0 0 1 1],'Title','Hopf Finder App',fa,'italic',fw,'bold',fs,fs_num);

%Tab Group
tg1 = uitabgroup('Parent',ph1,'Units',nzd,pst,[0.01 .675 .5 .3]);
tab1 = uitab(tg1,'Title','Setup','TooltipString','Setup problem');
tab2 = uitab(tg1,'Title','Configuration','TooltipString','Configure solver');
tg1.SelectedTab = tab1;

%Tab1
%----------------------------------------------------
%Model
uicontrol('Parent',tab1,'Units',nzd,pst,[0.025 .89 .10 .075],'Style','text',str,'Model',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab1,'Units',nzd,pst,[0.025 .70 .63 .175],'Style','edit',bgc,pbc,str,[cd '/*_model.m'],'Tag','selModel',hal,'left',fa,'italic',fgc,pfgc2,'Enable','on')
uicontrol('Parent',tab1,'units',nzd,pst,[0.830-1/6 .70 1/6 .175],'Style','pushbutton',str,'Edit',fs,fs_num,fw,'bold',fa,'normal','tag','edit_model','Enable','on')
uicontrol('Parent',tab1,'units',nzd,pst,[0.830 .70 1/6 .175],'Style','pushbutton',str,'Browse',fs,fs_num,fw,'bold',fa,'normal','tag','selModel_btn','Enable','on')

%X0
uicontrol('Parent',tab1,'Units',nzd,pst,[0.025 .555 .22 .075],'Style','text',str,'X0',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab1,'units',nzd,pst,[0.025 .37 .22 .175],'Style','edit',str,'[0,0]',fs,fs_num,fw,'normal',fa,'italic','tag','x0',hal,'right',fgc,pfgc2,'Enable','on')

%Par0
uicontrol('Parent',tab1,'Units',nzd,pst,[0.27 .555 .22 .075],'Style','text',str,'Par0',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab1,'units',nzd,pst,[0.27 .37 .22 .175],'Style','edit',str,'[0,0]',fs,fs_num,fw,'normal',fa,'italic','tag','par0',hal,'right',fgc,pfgc2,'Enable','on')

%LB
uicontrol('Parent',tab1,'Units',nzd,pst,[0.515 .555 .22 .075],'Style','text',str,'Lower bound',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab1,'units',nzd,pst,[0.515 .37 .22 .175],'Style','edit',str,'[0,0]',fs,fs_num,fw,'normal',fa,'italic','tag','lb',hal,'right',fgc,pfgc2,'Enable','on')

%UB
uicontrol('Parent',tab1,'Units',nzd,pst,[0.76 .555 .22 .08],'Style','text',str,'Upper bound',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab1,'units',nzd,pst,[0.76 .37 .22 .175],'Style','edit',str,'[0,0]',fs,fs_num,fw,'normal',fa,'italic','tag','ub',hal,'right',fgc,pfgc2,'Enable','on')

%Solver
intList = {'ode23 ','ode45','ode113','ode23t','ode23s','ode15s','ode23tb','ode15i','dasslc'};
ssList = {'levenberg-marquardt','trust-region-dogleg','trust-region-reflective'};

%Time Interval
uicontrol('Parent',tab1,'Units',nzd,pst,[0.025 .22 .22 .075],'Style','text',str,'Time interval',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab1,'units',nzd,pst,[0.025 .03 .22 .175],'Style','edit',str,'[0,0]',fs,fs_num,fw,'normal',fa,'italic','Tag','it',hal,'right',fgc,pfgc2,'Enable','on',fs,fs_num)

%Integrator
uicontrol('Parent',tab1,'Units',nzd,pst,[0.27 .22 .22 .075],'Style','text',str,'DAE solver',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab1,'Units',nzd,pst,[0.27 .12 .22 .085],'Style','popupmenu',str,intList,fs,fs_num,'Value',1,'Tag','odeSolver','Value',6)

%Steady-state solver
uicontrol('Parent',tab1,'Units',nzd,pst,[0.515 .22 .26 .075],'Style','text',str,'Steady-state solver',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab1,'Units',nzd,pst,[0.515 .12 .26 .085],'Style','popupmenu',str,ssList,fs,fs_num,'Value',1,'Tag','ssSolver')

%GetRes
uicontrol('Parent',tab1,'units',nzd,pst,[0.830 .03 1/6 .175],'Style','pushbutton',str,'Get Results',fs,fs_num,fw,'bold',fa,'normal','tag','getRes','Enable','off')

%Tab2
%----------------------------------------------------
%TolFun
uicontrol('Parent',tab2,'Units',nzd,pst,[0.025 .89 .22 .075],'Style','text',str,'TolFun',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.025 .7 .22 .175],'Style','edit',str,eps,fs,fs_num,fw,'normal',fa,'italic','tag','TolFun',hal,'right',fgc,pfgc2,'Enable','on')

%TolX
uicontrol('Parent',tab2,'Units',nzd,pst,[0.27 .89 .22 .075],'Style','text',str,'TolX',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.27 .7 .22 .175],'Style','edit',str,eps,fs,fs_num,fw,'normal',fa,'italic','tag','TolX',hal,'right',fgc,pfgc2,'Enable','on')

%ObjLimit
uicontrol('Parent',tab2,'Units',nzd,pst,[0.515 .89 .22 .075],'Style','text',str,'Objective Limit',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.515 .7 .22 .175],'Style','edit',str,1e-2,fs,fs_num,fw,'normal',fa,'italic','tag','ObjectiveLimit',hal,'right',fgc,pfgc2,'Enable','on')

%MaxTrials
uicontrol('Parent',tab2,'Units',nzd,pst,[0.76 .89 .22 .075],'Style','text',str,'Max Trials',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.76 .7 .22 .175],'Style','edit',str,1,fs,fs_num,fw,'normal',fa,'italic','tag','MaxTrials',hal,'right',fgc,pfgc2,'Enable','on')

%nger
uicontrol('Parent',tab2,'Units',nzd,pst,[0.025 .555 .22 .075],'Style','text',str,'MaxIter',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.025 .37 .22 .175],'Style','edit',str,100,fs,fs_num,fw,'normal',fa,'italic','tag','nger',hal,'right',fgc,pfgc2,'Enable','on')

%npas
uicontrol('Parent',tab2,'Units',nzd,pst,[0.27 .555 .22 .075],'Style','text',str,'# of birds',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.27 .37 .22 .175],'Style','edit',str,100,fs,fs_num,fw,'normal',fa,'italic','tag','npas',hal,'right',fgc,pfgc2,'Enable','on')

%MaxFunEvals
uicontrol('Parent',tab2,'Units',nzd,pst,[0.515 .555 .22 .075],'Style','text',str,'MaxFunEvals',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.515 .37 .22 .175],'Style','edit',str,1e6,fs,fs_num,fw,'normal',fa,'italic','tag','MaxFunEvals',hal,'right',fgc,pfgc2,'Enable','on')

%InitialT
uicontrol('Parent',tab2,'Units',nzd,pst,[0.76 .555 .22 .075],'Style','text',str,'Guess Period',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.76 .37 .22 .175],'Style','edit',str,'[3.1416, 1.4142]',fs,fs_num,fw,'normal',fa,'italic','tag','fac',hal,'right',fgc,pfgc2,'Enable','on')

%Atol
uicontrol('Parent',tab2,'Units',nzd,pst,[0.025 .22 .22 .075],'Style','text',str,'AbsTol',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.025 .03 .22 .175],'Style','edit',str,1e-8,fs,fs_num,fw,'normal',fa,'italic','tag','AbsTol',hal,'right',fgc,pfgc2,'Enable','on')

%Rtol
uicontrol('Parent',tab2,'Units',nzd,pst,[0.27 .22 .22 .075],'Style','text',str,'RelTol',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.27 .03 .22 .175],'Style','edit',str,1e-8,fs,fs_num,fw,'normal',fa,'italic','tag','RelTol',hal,'right',fgc,pfgc2,'Enable','on')

%MeshRef
uicontrol('Parent',tab2,'Units',nzd,pst,[0.515 .22 .22 .075],'Style','text',str,'# of Mesh points',hal,'left',fw,'bold',fs,fs_num)
uicontrol('Parent',tab2,'units',nzd,pst,[0.515 .03 .22 .175],'Style','edit',str,'[100, 100]',fs,fs_num,fw,'normal',fa,'italic','tag','MeshRef',hal,'right',fgc,pfgc2,'Enable','on')

%UsePast
uicontrol('Parent',tab2,'units',nzd,pst,[0.76 .15 .22 .125],'Style','checkbox',str,'Use Past SS',fs,fs_num,fw,'normal',fa,'italic','tag','usePast','Value',1)
uicontrol('Parent',tab2,'units',nzd,pst,[0.76 .0 .22 .125],'Style','checkbox',str,'XXX',fs,fs_num,fw,'normal',fa,'italic','tag','XXX','Value',0,'Enable','off')

%----------------------------------------------
%Create Results panel
ph2 = uipanel('Parent',fh,'Units',nzd,pst,[0.525 0.6625 .455 .3],'Title','Calculation & Results',fa,'italic',fw,'bold',fs,fs_num);

%Hybrid
uicontrol('Parent',ph2,'units',nzd,pst,[0.025 .75 .3 .175],'Style','pushbutton',str,'Hybrid',fs,fs_num,fw,'bold',fa,'normal','tag','pushHybrid','Enable','on','Visible','on')
uicontrol('Parent',ph2,'units',nzd,pst,[0.025 .75 .3 .175],'Style','togglebutton',str,'Stop',fs,fs_num,fw,'bold',fa,'normal','tag','stopHybrid','Enable','off','Visible','off','Value',0)

%Swarm
uicontrol('Parent',ph2,'units',nzd,pst,[0.025 .525 .3 .175],'Style','pushbutton',str,'Swarm',fs,fs_num,fw,'bold',fa,'normal','tag','pushSwarm','Enable','on','Visible','on')
uicontrol('Parent',ph2,'units',nzd,pst,[0.025 .525 .3 .175],'Style','togglebutton',str,'Stop',fs,fs_num,fw,'bold',fa,'normal','tag','stopSwarm','Enable','off','Visible','off','Value',0)

%Fminsearch
uicontrol('Parent',ph2,'units',nzd,pst,[0.025 .275 .3 .175],'Style','pushbutton',str,'Simplex',fs,fs_num,fw,'bold',fa,'normal','tag','pushFmin','Enable','on','Visible','on')
uicontrol('Parent',ph2,'units',nzd,pst,[0.025 .275 .3 .175],'Style','togglebutton',str,'Stop',fs,fs_num,fw,'bold',fa,'normal','tag','stopSimplex','Enable','off','Visible','off','Value',0)

%Levenberg-Marquardt
uicontrol('Parent',ph2,'units',nzd,pst,[0.025 .025 .3 .175],'Style','pushbutton',str,'Levenberg-Marquardt',fs,fs_num,fw,'bold',fa,'normal','tag','pushLM','Enable','on','Visible','on')
uicontrol('Parent',ph2,'units',nzd,pst,[0.025 .025 .3 .175],'Style','togglebutton',str,'Stop',fs,fs_num,fw,'bold',fa,'normal','tag','stopLM','Enable','off','Visible','off','Value',0)

%Results view
uicontrol('Parent',ph2,'units',nzd,pst,[0.35 .025 .625 .9],'Style','edit',str,'',fs,fs_num,'tag','resView','Max',8,'Enable','Inactive',hal,'center','FontName',myFont)

%Plot panel
ph3 = uipanel('Parent',fh,'Units',nzd,pst,[0.01 0.01 .97 .63],'Title','Simulation',fa,'italic',fw,'bold',fs,fs_num);

%axes
axes('Parent',ph3,'Units',nzd,pst,[.075 .15 .675 .8],'Tag','ax','Box','on','YGrid','on','FontName','Palatino',fs,fs_num+2);

%Simulation
uicontrol('Parent',ph3,'units',nzd,pst,[0.76 .875 .2275 .08],'Style','pushbutton',str,'Simulation',fs,fs_num,fw,'bold',fa,'normal','tag','simulate','Enable','on')
uicontrol('Parent',ph3,'units',nzd,pst,[0.76 .875 .2275 .08],'Style','togglebutton',str,'Stop',fs,fs_num,fw,'bold',fa,'normal','tag','stopSimulation','Enable','off','Visible','off','Value',0)

%Progress Bar SIM
uicontrol('Parent',ph3,'Units',nzd,pst,[0.7615 .82 .2255 .05],'Style','frame',fgc,[0 0 0],bgc,[1 1 1],'Tag','pbar_base')
uicontrol('Parent',ph3,'Units',nzd,pst,[0.7615 .82 .2255*eps .05],'Style','frame',fgc,[0 0 0],bgc,[1 0 0],'Tag','pbar')

%Plot shuffle
uicontrol('Parent',ph3,'units',nzd,pst,[0.76 .765 .2275 .05],'Style','pushbutton',str,'Next Plot',fs,fs_num,fw,'bold',fa,'normal','tag','nextPlot','Enable','off')
uicontrol('Parent',ph3,'units',nzd,pst,[0.76 .71 .11 .05],'Style','pushbutton',str,'Export Figure',fs,fs_num,fw,'bold',fa,'normal','tag','exportFigure','Enable','off')
uicontrol('Parent',ph3,'units',nzd,pst,[0.8775 .71 .11 .05],'Style','pushbutton',str,'Plot All',fs,fs_num,fw,'bold',fa,'normal','tag','plotAll','Enable','off')

%Build Objective function
uicontrol('Parent',ph3,'units',nzd,pst,[0.76 .57 .2275 .08],'Style','pushbutton',str,'Build Objective Function',fs,fs_num,fw,'bold',fa,'normal','tag','buildObjective','Enable','on','Visible','on')
uicontrol('Parent',ph3,'units',nzd,pst,[0.76 .57 .2275 .08],'Style','togglebutton',str,'Stop',fs,fs_num,fw,'bold',fa,'normal','tag','stopBuilding','Enable','off','Visible','off','Value',0)

%Progress Bar OF
uicontrol('Parent',ph3,'Units',nzd,pst,[0.7615 .515 .2255 .05],'Style','frame',fgc,[0 0 0],bgc,[1 1 1],'Tag','pbar_base2')
uicontrol('Parent',ph3,'Units',nzd,pst,[0.7615 .515 .2255*eps .05],'Style','frame',fgc,[0 0 0],bgc,[1 0 0],'Tag','pbar2')

%Load Save Mesh
uicontrol('Parent',ph3,'units',nzd,pst,[0.76 .46 .11 .05],'Style','pushbutton',str,'Load Mesh',fs,fs_num,fw,'bold',fa,'normal','tag','loadMesh','Enable','on')
uicontrol('Parent',ph3,'units',nzd,pst,[0.8775 .46 .11 .05],'Style','pushbutton',str,'Save Mesh',fs,fs_num,fw,'bold',fa,'normal','tag','saveMesh','Enable','off')

%Save Load Reset
uicontrol('Parent',ph3,'units',nzd,pst,[0.76 .32 .2275 .08],'Style','pushbutton',str,'Load Results',fs,fs_num,fw,'bold',fa,'normal','tag','load','Enable','on')
uicontrol('Parent',ph3,'units',nzd,pst,[0.76 .235 .2275 .08],'Style','pushbutton',str,'Save Results',fs,fs_num,fw,'bold',fa,'normal','tag','save','Enable','off')
uicontrol('Parent',ph3,'units',nzd,pst,[0.76 .15 .2275 .08],'Style','pushbutton',str,'Reset',fs,fs_num,fw,'bold',fa,'normal','tag','reset','Enable','on')

%Create handles superstructure
h = guihandles(fh); 

%Initialize variables
h.num.model = []; h.num.x0 = []; h.num.par0 = []; h.num.lb = []; h.num.ub = []; h.num.it = [];
h.num.TolFun = eps; h.num.TolX = eps; h.num.ObjectiveLimit = 1e-2; h.num.MaxTrials = 1;
h.num.nger = 100; h.num.npas = 100; h.num.MaxFunEvals = 1e6; h.num.initialT = [pi sqrt(2)];
h.num.AbsTol = 1e-8;  h.num.RelTol = 1e-8; h.getError = -3; h.fs_num = fs_num;
h.plotState = []; h.isMeshLoaded = 0; h.idScreen = Pix_SS(3); h.num.isSimulated = 0;
h.num.MeshRef = [100,100];
h.intList = intList;
h.ssList = ssList;
h.Version = V;
h.isEigenMap = 0;
h.isBranchDiagram = 0;

%Address Callbacks
[h.TolFun.Callback, h.TolX.Callback, h.ObjectiveLimit.Callback, h.MaxTrials.Callback, h.nger.Callback, ...
    h.npas.Callback, h.MaxFunEvals.Callback, h.AbsTol.Callback, h.RelTol.Callback] = deal(@getPositiveNumber);

[h.x0.Callback, h.par0.Callback, h.lb.Callback, h.ub.Callback, h.it.Callback, h.MeshRef.Callback, h.fac.Callback] = deal(@getVector);

h.selModel_btn.Callback = @selBtnPress;

[h.pushSwarm.Callback, h.pushHybrid.Callback, h.pushFmin.Callback, h.pushLM.Callback] = deal(@run);

% [h.nameVary.ButtonDownFcn, h.nameVarx.ButtonDownFcn] = deal(@nameVariables);

h.simulate.Callback = @simulate;

h.plotAll.Callback = @plotAllFunction;

h.nextPlot.Callback = @plotNextFunction;

h.exportFigure.Callback = @exportFigureFunction;

h.buildObjective.Callback = @buildObjectiveFunction;

[h.reset.Callback, h.mreset.Callback]  = deal(@resetBtnPress);

h.mexample.Callback = @openExample;

h.saveMesh.Callback = @saveMeshFunction;

h.loadMesh.Callback = @loadMeshFunction;

h.edit_model.Callback = @editFunction;

[h.save.Callback, h.msave.Callback] = deal(@saveResultsFunction);

[h.load.Callback, h.mload.Callback] = deal(@loadResultsFunction);

h.mquit.Callback = @quitFunction;   

h.selModel.Callback = @selModelFunction;

h.getRes.Callback = @getResFunction;

h.eigenmap.Callback = @eigenMapFunction;

h.branchdiagram.Callback = @oneParCont;

h.meditme.Callback = @editme;


%Save superstructure
guidata(fh,h)

%Set figure visible
fh.Visible = 'on';


end

%% Run Callback Function
function run(ho,~)
%Load handles superstructure
h = guidata(ho);

%Check for setup errors
switch h.getError
    case -1
        errordlg('Some information might be unset or wrong. Check for error','Error','modal')
        return
    case -2 
        errordlg('Invalid Model','Error','modal')
        return
    case {-3,-4}
        errordlg('Problem unsetted','Error','modal')
        return

    otherwise
        
        %Check model setup consistency
        numVar = length(h.num.x0);
        numPar = length(h.num.par0);
        try
            if nargout(h.num.model) == 2
                [numF,MassM] = feval(h.num.model,[],h.num.x0,h.num.par0);
            else
                numF = feval(h.num.model,[],h.num.x0,h.num.par0);
                MassM = eye(numVar);
            end
        catch
            errordlg('Check model! Some inconsistency was identified.')
            return
        end
        
        [Bn,Bm] = size(MassM);
        if Bn ~= Bm
            err = sprintf('Mass matrix must be square.');
            errordlg(err,'Error','modal');
            return
        end

        if length(numF) ~= numVar
            err = sprintf('Model returns a vector of length %d, but the length of initial conditions vector is %d',[length(numF),numVar]);
            errordlg(err,'Error','modal')
            return
        end
        
        if length(numF) ~= Bn
            err = sprintf('Dimensions of mass matrix and system must agree.');
            errordlg(err,'Error','modal');
            return
        end
        
        %Check bound consistency
        if length(h.num.lb) ~= numPar || length(h.num.ub) ~= numPar
            errordlg('Lower and upper bounds dimensions do not agree','Error','modal')
            return
        elseif any(h.num.lb >= h.num.ub)
            errordlg('Inconsistent bounds','Error','modal')
            return
        end
        
        if isempty(h.num.it)
            errordlg('Missing integration interval.','Error','modal')
            return
        elseif length(h.num.it) ~= 2
            errordlg('Please, provide only initial and final integration time.','Error','modal')
            return
        end
            
        if h.num.it(end) <= h.num.it(1)
            errordlg('The last entry in tspan must be greater than the first one','Error','modal')
            return
        end
        
        %Check SS solvability
        testF = h.num.model([],h.num.x0,h.num.par0);
        if any(isnan(testF))
            errordlg('The model return NaN with given initial guess. Steady-State cannot be found.')
            return
        end
end

%Now everything should be OK -> Disable interface
set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
    h.save, h.load, h.reset, h.loadMesh, h.saveMesh, h.buildObjective,  ...
    h.nextPlot, h.plotAll, h.exportFigure, h.getRes],'Enable','off')

set(h.resView,'Enable','Inactive')
set(h.pbar,'Position',[0.762 .515 .225*eps .05],'BackgroundColor',[1 0 0])
set(h.pbar2,'Position',[0.7615 .515 .2255*eps .05],'BackgroundColor',[1 0 0])
cla(h.ax)
drawnow;

% =========================================================================
% Setup problem
% =========================================================================
notSS = 0;
opt_ss = optimoptions(@fsolve,'algorithm',h.ssList{get(h.ssSolver,'Value')},'TolFun',h.num.TolFun,...
                        'TolX',h.num.TolX,'Jacobian','off','Display','off','MaxFunEvals',h.num.MaxFunEvals,...
                        'MaxIter',h.num.nger*h.num.npas);
fun = @(state,par) h.num.model([],state,par);
x0 = h.num.x0(:);
lb = h.num.lb(:);
ub = h.num.ub(:);
par0 = h.num.par0(:);
nS = 0;

% Find DAE index by the charpoly order at initial point
J_ss0 = zeros(numVar,numVar);
I = eye(numVar);
f0 = feval(fun,x0,par0);
for m = 1:numVar
    J_ss0(:,m) = ( feval(fun,x0+1e-6*I(:,m),par0) - f0 )/1e-6;
end
r = rank(MassM);
P0 = labudde(J_ss0,MassM);
nPol = length(P0) - 1;

sol.ExplicitAE = (numVar - r)*(r ~= numVar);
sol.HiddenAE = (r - nPol)*(r ~= numVar);
sol.index =  (r - nPol + 1)*(r ~= numVar);

if nPol < 2
    errordlg('This system appears to be of dimension 1. Periodic orbits can''''t arise in such system! ','Error','modal')
    return
end

switch get(ho,'String')
    case 'Swarm'  
        set(h.stopSwarm,'Enable','on','Visible','on')
        s_initial{1} = sprintf('Running optimization.');
        set(h.resView,'String',s_initial,'HorizontalAlignment','center'); drawnow;
        
        tic
        [sol.par,sol.f,sol.nS] = mySwarm(@objective,par0,lb,ub,h.num.nger,h.num.npas,h.num.ObjectiveLimit,...
                                            h.num.MaxFunEvals,h.resView,h.stopSwarm,h.stopHybrid,s_initial,h.idScreen);
        sol.time = toc;
        set(h.stopSwarm,'Enable','off','Visible','off','Value',0) 
        
    case 'Simplex'
        set(h.stopSimplex,'Enable','on','Visible','on');
        optfmin = optimset('Display','none','TolFun',h.num.TolFun,'TolX',h.num.TolX,...
                                'MaxFunEvals',h.num.MaxFunEvals,'MaxIter',h.num.nger,'OutputFcn',@fsolveOutput);
        s_initial{1} = sprintf('Running optimization.');
        s_initial{2} = sprintf('\n F-Count               Best value'); 
        ig = 2;
        viewIndex = floor(h.idScreen/220);
        set(h.resView,'String',s_initial,'HorizontalAlignment','center'); drawnow;
        tic
        [sol.par,sol.f,~,output] = fminsearch(@objective,par0,optfmin);
        sol.time = toc;
        sol.f = sol.f;
        sol.nS = output.funcCount;
        set(h.stopSimplex,'Enable','off','Visible','off','Value',0);
        
    case 'Levenberg-Marquardt'
        
        set(h.stopLM,'Enable','on','Visible','on');
        optLM = optimoptions(@fsolve,'Display','off','TolFun',h.num.TolFun,'TolX',h.num.TolX,...
                                'MaxFunEvals',h.num.MaxFunEvals,'MaxIter',h.num.nger,...
                                'OutputFcn',@fsolveOutput,'Algorithm','levenberg-marquardt');
        s_initial{1} = sprintf('Running optimization.');
        s_initial{2} = sprintf('\n F-Count               Best value'); 
        ig = 2;
        viewIndex = floor(h.idScreen/220);
        set(h.resView,'String',s_initial,'HorizontalAlignment','center'); drawnow;
        tic
        [sol.par,sol.f,~,output] = fsolve(@objective,par0,optLM);
        sol.time = toc;
        sol.f = sol.f;
        sol.nS = output.funcCount;
        set(h.stopLM,'Enable','off','Visible','off','Value',0);
        
    case 'Hybrid'
        
        set(h.stopHybrid,'Enable','on','Visible','on')
        s_initial{1} = sprintf('Starting Swarm Optimization.');
        s_initial{2} = sprintf('\n F-Count               Best value'); 
        set(h.resView,'String',s_initial,'HorizontalAlignment','center'); drawnow;
        optfmin = optimset('Display','none','TolFun',h.num.TolFun,'TolX',h.num.TolX,...
                                'MaxFunEvals',h.num.MaxFunEvals,'MaxIter',h.num.nger,'OutputFcn',@fsolveOutput);
        
        tic
        [parSwarm,~,nSSwarm] = mySwarm(@objective,par0,lb,ub,h.num.nger,h.num.npas,h.num.ObjectiveLimit,...
                                            h.num.MaxFunEvals,h.resView,h.stopSwarm,h.stopHybrid,s_initial,h.idScreen);
                                        
        s_initial{1} = sprintf('Running optimization.');
        s_initial{2} = sprintf('\n F-Count               Best value'); 
        ig = 2;
        viewIndex = floor(h.idScreen/220);
        set(h.resView,'String',s_initial,'HorizontalAlignment','center'); drawnow;

        [sol.par,sol.f,~,output] = fminsearch(@objective,parSwarm,optfmin); 
        sol.time = toc;
        sol.f = sol.f;
        sol.nS = output.funcCount + nSSwarm;
        
%         Problem.f = @objective;
%         bounds = [lb ub];
%         opt.maxevals = h.num.MaxFunEvals;
%         opt.maxits = h.num.nger;
%         
%         tic
%         [sol.f,sol.par,sol.history] = direct(Problem,bounds,opt); 
%         sol.time = toc;
%         sol.nS = sol.history(end,2);
        set(h.stopHybrid,'Enable','off','Visible','off','Value',0);
        drawnow
        
end

sol.notSS = notSS;
sol.OptType = get(ho,'String');
[sol.x,~,~,~,sol.Jac] = fsolve(fun,x0,opt_ss,sol.par);
[sol.derivative,~] = check_der(sol.par);
sol.eig = roots(labudde(real(sol.Jac),MassM));
if isempty(sol.derivative)
    sol.derivative = NaN;
end

% =========================================================================
% Build periodic orbit
% =========================================================================
if sol.f > 1e-8
    s_hopfFound = sprintf(['Objective Function = %.4e. The steady point may not be a Hopf.\n'...
                                'Try to compute periodic orbit anyway?'],sol.f);
else
    s_hopfFound = sprintf('Objective Function = %.4e. \nProceed to compute periodic orbit?',sol.f);
end

choice = questdlg(s_hopfFound,'Periodic orbit search','Yes','No','Yes');

if strcmp('Yes',choice)
    % Check stability
    s_building{1} = sprintf('\nHopf found! Checking orbit stability ...');
    set(h.resView,'String',s_building,'HorizontalAlignment','center','FontSize',h.fs_num,'Enable','on')
    drawnow

    odeSolver = str2func(h.intList{get(h.odeSolver,'Value')});
    if det(MassM)
        optOde = odeset('RelTol',1e-8,'AbsTol',1e-8);
    else
        optOde = odeset('RelTol',1e-8,'AbsTol',1e-8,'Mass',MassM);
    end

    model = h.num.model;
    [~,y_stable] = feval(odeSolver,model,[0 min(100*h.num.it(end),1000)],sol.x+0.1*norm(sol.x),optOde,sol.par);

    windowFat = 10;
    windowLength = fix(length(y_stable(:,1))/windowFat);

    mm = zeros(windowFat,1);
    for i = 1:windowFat
        mm(i) = max(y_stable((i-1)*windowLength+1:i*windowLength,1));
    end

    mmSense = sort(abs(diff(mm)./mm(1:end-1)*100));
    if mmSense(1) >= 0.01
        sol.stability = 'Unstable';
    else
        sol.stability = 'Stable';
    end
    opt3 = optimset('TolFun',eps,'TolX',eps,'MaxIter',100,'Display','iter');
    opt2 = optimoptions(@fsolve,'Display','iter-detailed','TolFun',h.num.TolFun,'TolX',h.num.TolX,...
                        'MaxFunEvals',h.num.MaxFunEvals,'MaxIter',h.num.nger,'Algorithm','levenberg-marquardt');
    f1 = 1; i = 1;
%     x = [h.num.initialT(2)*sol.x;2*pi/max(imag(sol.eig(imag(sol.eig)~=0)))];
    x = [eye(numVar,1);2*pi/max(imag(sol.eig(imag(sol.eig)~=0)))];


    s_building{1} = sprintf('\nHopf found! Computing periodic solution.');

    maxTries = 10;
    while norm(f1) > 1000*eps
        s_building{2} = sprintf('\n Progress: %d%%.',(i-1)/maxTries*100);
        set(h.resView,'String',s_building,'HorizontalAlignment','center','FontSize',h.fs_num,'Enable','on')
        drawnow
    
        [x,f1] = fminsearch(@myStableOrbit,x,opt3,1);
        [x,f2] = fsolve(@myStableOrbit,x,opt2,0);

        if i >= maxTries || norm(f1) <= 10000*eps || norm(f2-f1) < 100*eps
            break
        end
        i = i+1;
    end
    sol.x_T = x(1:end-1); sol.T = x(end);
end

% =========================================================================
% Display results
% =========================================================================
if h.idScreen == 1920
    s_final{1}  = sprintf('\n*** Opt Finished & Results ***            ');
    s_final{2}  = sprintf('\nFile ............................................. %s',[func2str(h.num.model) '.m']);
    s_final{3}  = sprintf(  'Index ............................................ %d ', sol.index);
    s_final{4}  = sprintf(  'Fval ............................................. %.6e',sol.f);

    s = repmat('%.4e, ',1,length(sol.par)); s(end-1:end) = [];
    s_final{5}  = sprintf([ 'Par .............................................. [' s ']'], sol.par);
    
    s = repmat('%.4e, ',1,length(sol.x)); s(end-1:end) = [];
    s_final{6}  = sprintf([ 'Ss ............................................... [' s ']'], sol.x);
    
    s_final{7}  = sprintf(  'Time ............................................. %.4f s',sol.time);
    s_final{8}  = sprintf(  'FCount ........................................... %d',sol.nS);
%     s_final{9}  = sprintf(  'Derivative ....................................... %.4f',sol.derivative(1));
    s_final{9}  = sprintf(  'NotSS ............................................ %d ',sol.notSS);
    s_final{10} = sprintf(  'OptType .......................................... %s',sol.OptType);
    
    s = repmat('%.4e %+.4ei, ',1,length(sol.eig)); s(end-1:end) = [];
    s_final{11} = sprintf([ 'Eigenvalues ...................................... [' s ']'], [real(sol.eig) imag(sol.eig)]');

    if isfield(sol,'T')
        s_final{12} = sprintf( 'Orbit Stability................................... %s',sol.stability);
        s_final{13} = sprintf( 'Period ........................................... %.4f',sol.T);
        s = repmat('%.4e, ',1,length(sol.x_T)); s(end-1:end) = [];
        s_final{14} = sprintf(['Pp ............................................... [' s ']\n'], sol.x_T);
    end
    
else
    s_final{1} = sprintf('\n*** Opt Finished & Results ***      ');
    s_final{2} = sprintf( '\nFile ................................. %s ',[func2str(h.num.model) '.m']);
    s_final{3} = sprintf( 'Index ................................ %d ', sol.index);
    s_final{4} = sprintf( 'Fval ................................. %.6e ',sol.f);

    s = repmat('%.4e, ',1,length(sol.par)); s(end-1:end) = [];
    s_final{5} = sprintf(['Par .................................. [' s '] '], sol.par);
    
    s = repmat('%.4e, ',1,length(sol.x)); s(end-1:end) = [];
    s_final{6} = sprintf(['Ss ................................... [' s '] '], sol.x);
    s_final{7} = sprintf( 'Time ................................. %.4f s ',sol.time);
    s_final{8} = sprintf( 'FCount ............................... %d ',sol.nS);
%     s_final{9} = sprintf( 'Derivative ........................... %.4f ',sol.derivative(1));
    s_final{9} = sprintf('NotSS ................................ %d ',sol.notSS);
    s_final{10} = sprintf('OptType .............................. %s ',sol.OptType);
    
    s = repmat('%.4e %+.4ei, ',1,length(sol.eig)); s(end-1:end) = [];
    s_final{11} = sprintf([ 'Eigenvalues .......................... [' s '] '], [real(sol.eig) imag(sol.eig)]');
    
    if isfield(sol,'T')
        s_final{12} = sprintf( 'Orbit Stability....................... %s ',sol.stability);
        s_final{13} = sprintf( 'Period ............................... %.4f ',sol.T);
        s = repmat('%.4e, ',1,length(sol.x_T)); s(end-1:end) = [];
        s_final{14} = sprintf(['Pp ................................... [' s '] \n'], sol.x_T);
    end
end

set(h.resView,'String',s_final,'HorizontalAlignment','right','FontSize',h.fs_num,'Enable','on')

%Enable Interface
set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
    h.save, h.load, h.reset, h.buildObjective, h.getRes],'Enable','on')
drawnow;

h.num.sol = sol;
guidata(ho,h)

% =========================================================================
% Auxiliary nested functions
% =========================================================================

    function [f,w0] = check_der(par)
        
        [x_ss,~,~,~,J_ss] = fsolve(fun,x0,opt_ss,par);
        [~,~,~,~,J_ss2] = fsolve(fun,x_ss,opt_ss,par+1e-6);
        
        L = eig(J_ss,MassM);
        L2 = eig(J_ss2,MassM);
        
        R = imag(L)~=0;
                
        f = unique( min( abs( real(L2(R==1) - L(R==1))/1e-6 ) ).*sign(real(L2(R==1) - L(R==1))/1e-6) );
        w0 = unique ( abs( imag( L2(R==1) ) ) );
        
    end

    function f = objective(par)
        
        [x_ss,f_ss,status,~,J_ss] = fsolve(fun,x0,opt_ss,par);
        if status <= 0 
            f = 1000*pi;
            return
%             [x_ss,f_ss,status,~,J_ss] = hfsolve(fun,x0,opt_ss,5,par);
%             if status <= 0
%                 [x_ss,f_ss,~,~,J_ss] = hfsolve(fun,x0,opt_ss,30,par);
%             end
        end
        
        if norm(f_ss) > 1e-6 || any(imag(x_ss)~=0) || any(any(imag(J_ss)~=0))
            warning('Steady-state not found!');
            notSS = notSS + 1;
            f = 10000*pi;
            return
        elseif get(h.usePast,'Value')
            x0 = x_ss;
        end
        
        P = labudde(J_ss,MassM);
        nPol2 = length(P) - 1;
        if nPol ~= nPol2
            warning('Order changed!');
            fprintf('\nStarting n = %d and Current n = %d\n',nPol,nPol2);
            f = 10000*pi/2;
            return
        end
        
        for j = nPol2-1:-1:1
            for k = 1:nPol2-1
                pos = (2*j-k+1);
                if pos > nPol2+1 || pos < 1 
                    Delta(j,k) = 0;
                else
                    Delta(j,k) = P(pos);
                end
            end
        end
        
        con = -min(0,P(end));
        bar = -min(0,min(par - lb)) - min(0,min(ub - par));
        
        f = norm(f_ss) + norm(det(Delta)) + pi*1e6*bar + pi*1e6*con;
    end

    function stop = fsolveOutput(~, optimValues, ~)
        
        stop1 = logical(get(h.stopSimplex,'Value'));
        stop2 = logical(get(h.stopLM,'Value'));
        stop = (stop1 || stop2);
        if ~mod(optimValues.funccount,10)
            ig = ig+1;
            s{ig} = sprintf(' %7.0d %24.6e',[optimValues.funccount+nS,optimValues.fval]);
            if (ig > viewIndex)
                set(h.resView,'String',[s_initial s{ig-viewIndex:ig}]); drawnow;
            else
                set(h.resView,'String',[s_initial s{1:ig}]); drawnow;
            end

        end
    end

    function f = myStableOrbit(x,flag)
        
        y0 = x(1:end-1); T = x(end);
        yp0 = feval(model,[],y0,sol.par);

        [~,y] = feval(odeSolver,model,[0 T],y0,optOde,sol.par);
        
        pen = abs(max(y(:,1)) - min(y(:,1)));
        f = [ (y(1,:) - y(end,:))';
%               norm( feval(model,[],y(1,:),sol.par) - feval(model,[],y(end,:),sol.par) ) ] + 100*double(T<0);
              yp0(1)  ] + 100*double(T<0) + 100*double(pen<1e-2);

        if flag
            f = norm(f)^2;
        end
    end

end

%% Simulate Callback Functions
function simulate(ho,~)

h = guidata(ho);

%Check for setup errors
switch h.getError
    case -1
        errordlg('Some information might be unset or wrong. Check for error','Error','modal')
        return
        
    case -2 
        errordlg('Invalid Model','Error','modal')
        return
        
    case {-3,-4}
        errordlg('Problem unsetted','Error','modal')
        return
end

%Check model consistency
numVar = length(h.num.x0);
try
    if nargout(h.num.model) == 2
        [numF,MassM] = feval(h.num.model,[],h.num.x0,h.num.par0);
    else
        numF = feval(h.num.model,[],h.num.x0,h.num.par0);
        MassM = eye(numVar);
    end
catch
    errordlg('Check model! Some inconsistency was identified.')
    return
end

if length(numF) ~= length(h.num.x0)
    err = sprintf('Model returns a vector of length %d, but the length of initial conditions vector is %d',...
                            [length(numF),length(h.num.x0)]);
    errordlg(err,'Error','modal')
    return
end

%Check bound consistency
if length(h.num.lb) ~= length(h.num.par0(:)) || length(h.num.ub) ~= length(h.num.par0(:))
    errordlg('Lower and upper bounds dimensions do not agree','Error','modal')
    return
elseif any(h.num.lb >= h.num.ub)
    errordlg('Inconsistent bounds','Error','modal')
    return
end

if h.num.it(end) == h.num.it(1)
    errordlg('The last entry in tspan must be different from the first entry','Error','modal')
    return
end

%Disable interface
set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
    h.save, h.load, h.reset, h.loadMesh, h.saveMesh, h.buildObjective, ...
    h.nextPlot, h.plotAll, h.exportFigure],'Enable','off')

set(h.stopSimulation,'Enable','on','Visible','on')
set(h.resView,'Enable','Inactive')
set(h.pbar,'Position',[0.762 .515 .225*eps .05],'BackgroundColor',[1 0 0])
set(h.pbar2,'Position',[0.7615 .515 .2255*eps .05],'BackgroundColor',[1 0 0])
drawnow;

%Solve ODE System
odeSolver = str2func(h.intList{get(h.odeSolver,'Value')});

if det(MassM)
    options = odeset('AbsTol',h.num.AbsTol,'RelTol',h.num.RelTol,'OutputFcn',@theprogress);
else
    options = odeset('AbsTol',h.num.AbsTol,'RelTol',h.num.RelTol,'OutputFcn',@theprogress,'Mass',MassM);
end

y0 = h.num.x0; par0 = h.num.par0;

try
    [t,y] = feval(odeSolver,h.num.model,h.num.it,y0,options,par0);
catch
    err = sprintf(['DAE system cannot be integrated with specified tolerance.\n'...
        'Consider loosen it up or provide better consistent inital conditions!']);
    errordlg(err,'Error','modal')
    %Enable interface
    set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
        h.load, h.reset, h.loadMesh, h.buildObjective, h.resView ...
        h.nextPlot, h.plotAll, h.exportFigure, h.save, h.msave],'Enable','on')
    set(h.stopSimulation,'Enable','off','Visible','off','Value',0)
    return
end

h.num.t = t; h.num.y = y;
set(h.pbar,'BackgroundColor',[0 1 0],'Position',[0.7615 .82 .2255 .05])

%Plot first result
set(h.main,'nextplot','add','CurrentAxes',h.ax)
h1 = plot(t,y(:,1));
set(h.main,'nextplot','new')
set(h1,'LineWidth',1.5)
h.ax.YAxis.TickLabelFormat = '%,.4g';
% h.ax.XAxis.TickLabelFormat = '%';
set(h.ax,'YGrid','on','Box','on','FontSize',h.fs_num+2)
ylabel('State #1','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
xlabel('Time','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)

h.numVar = length(h.num.x0);
h.plotCounterI = 2;
h.plotCounterJ = 0;

%Enable interface
set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
    h.load, h.reset, h.loadMesh, h.buildObjective, h.resView ...
    h.nextPlot, h.plotAll, h.exportFigure, h.save, h.msave],'Enable','on')
set(h.stopSimulation,'Enable','off','Visible','off','Value',0)
h.num.isSimulated = 1; h.isEigenMap = 0;

guidata(ho,h)

    function status = theprogress(t,~,flag,~)

        if isempty(flag)
            if ~mod(round(t/h.num.it(end),2),0.2)  
                set(h.pbar,'Position',[0.7615 .82 .2255*round(t/h.num.it(end),2) .05])
                drawnow;
            end
        end
        if get(h.stopSimulation,'Value')
            status = 1;
        else
            status = 0;
        end
    end
end

function plotNextFunction(ho,~)

h = guidata(ho);

set(h.main,'nextplot','add','CurrentAxes',h.ax)
if ~h.isBranchDiagram
    if h.plotCounterJ == 0
        plot(h.num.t,h.num.y(:,h.plotCounterI),'LineWidth',1.5)
        set(h.main,'nextplot','new')
        ylabel(['State #' num2str(h.plotCounterI)],'FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        xlabel('Time','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        h.ax.YAxis.TickLabelFormat = '%,.4g';
    %     h.ax.XAxis.TickLabelFormat = [];
        set(h.ax,'YGrid','on','Box','on','FontSize',h.fs_num+2)
        h.plotCounterI = h.plotCounterI+1;
        if h.plotCounterI > h.numVar
            h.plotCounterJ = 2;
            h.plotCounterI = 1;
        end
    else
        pp = plot(h.num.y(:,h.plotCounterI),h.num.y(:,h.plotCounterJ),h.num.y(1,h.plotCounterI),... 
            h.num.y(1,h.plotCounterJ),'s',h.num.y(end,h.plotCounterI),h.num.y(end,h.plotCounterJ),'+');
        set(h.main,'nextplot','new')
        set(pp,'LineWidth',1.5  );
        set(pp(2),'MarkerFaceColor',[0.8500    0.3250    0.0980],'MarkerEdgeColor',[0.8500    0.3250    0.0980])
        set(pp(3),'MarkerFaceColor',[0.9290    0.6940    0.1250],'MarkerEdgeColor',[0.9290    0.6940    0.1250])
    %     set(pp(2),'MarkerSize',16)
        ylabel(['State #' num2str(h.plotCounterJ)],'FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        xlabel(['State #' num2str(h.plotCounterI)],'FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        h.ax.YAxis.TickLabelFormat = '%,.4g';
        h.ax.XAxis.TickLabelFormat = '%,.4g';
        set(h.ax,'YGrid','on','Box','on','FontSize',h.fs_num+2)
        h.plotCounterJ = h.plotCounterJ + 1;
        if  h.plotCounterJ > h.numVar
            h.plotCounterI = h.plotCounterI + 1;
            h.plotCounterJ = h.plotCounterI + 1;
        end
        if h.plotCounterI >= h.numVar
            h.plotCounterJ = 0;
            h.plotCounterI = 1;
        end
    end
else

    if h.stateSel >= h.numVar
        h.stateSel = 1;
    else
        h.stateSel = h.stateSel + 1;
    end
    pp = plot(h.par,h.stableC(:,h.stateSel),h.par,h.unstableC(:,h.stateSel));
    set(pp(1),'LineStyle','-','Color','k','LineWidth',1.5)
    set(pp(2),'LineStyle',':','Color','r','LineWidth',1)
    xlabel('Par #1','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
    ylabel(['State #' num2str(h.stateSel)],'FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
    set(h.ax,'YGrid','on','Box','on','FontSize',h.fs_num+2)


    
end

guidata(ho,h)
end

function plotAllFunction(ho,~)

h = guidata(ho);
iF = length(findobj(0,'type','figure'));

for i = 1:h.numVar
    newFig = figured(i+iF);
    plot(h.num.t,h.num.y(:,i),'LineWidth',1.5)
    ylabel(['State #' num2str(i)],'FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
    xlabel('Time','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
    set(newFig,'Units','normalized','Position',[0.3536 0.3454 0.4609 0.5583])
end

counter = 1;
for i = 1:h.numVar
    for j = i+1:h.numVar
        newFig = figured(h.numVar+counter+iF);
        pp = plot(h.num.y(:,i),h.num.y(:,j),h.num.y(1,i),h.num.y(1,j),'s',h.num.y(end,i),h.num.y(end,j),'+');
        set(pp(2),'MarkerFaceColor',[0.8500    0.3250    0.0980],'MarkerEdgeColor',[0.8500    0.3250    0.0980])
        set(pp(3),'MarkerFaceColor',[0.9290    0.6940    0.1250],'MarkerEdgeColor',[0.9290    0.6940    0.1250])
        set(pp,'LineWidth',1.5)
        ylabel(['State #' num2str(j)],'FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        xlabel(['State #' num2str(i)],'FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        set(newFig,'Units','normalized','Position',[0.3536 0.3454 0.4609 0.5583])
        counter = counter + 1;
    end
end

guidata(ho,h)
end

function exportFigureFunction(ho,~)

h = guidata(ho);
axesObject = h.ax;

axes_units = 'normalized';
newFig = figured; delete(newFig.Children(2));  

axesObject2 = copyobj(axesObject,newFig);

set(axesObject2,'Units',axes_units);
set(axesObject2,'Position',[0.1 0.1 .85 .85]);
set(newFig,'Units','normalized','Position',[0.3 0.25 0.4 0.5])

if h.isEigenMap
    cursorMode = datacursormode(newFig);
    hPlot = get(axesObject2,'Children');
    nLines = length(hPlot);
    hDatatip{nLines-2,1} = 0;
    XData{nLines-2,1} = 0;
    YData{nLines-2,1} = 0;
    set(cursorMode,'UpdateFcn',{@eigenCursor,h.eigenMap},'DisplayStyle','datatip')
    k = 1;
    rec = 0;
    for i = 1:nLines
        if length(hPlot(i).XData) > 2
            XData{k} = hPlot(i).XData;
            YData{k} = hPlot(i).YData;
            hDatatip{k} = cursorMode.createDatatip(hPlot(i));
            set(hDatatip{k}, 'Marker','o', 'MarkerSize',10, 'MarkerFaceColor',hPlot(i).Color,...
                     'MarkerEdgeColor','w', 'OrientationMode','auto',...
                     'Position', [hPlot(i).XData(1) hPlot(i).YData(1) 0])
            updateDataCursors(cursorMode)
            drawnow;
            k = k+1;
        end
    end
    set(axesObject2,'FontSize',16,'Position',[0.13 0.11 0.775 0.815])
    rec = 1;
    ind = 1;
%     set(cursorMode, 'enable', 'off')
end

guidata(ho,h)

    function txt = eigenCursor(~,event_obj,data)

        % cursorMode = datacursormode(gcf);
        pos = get(event_obj,'Position');
        I = get(event_obj, 'DataIndex');
        txt = {['Re : ',sprintf('% 2.4g',pos(1))],...
               ['Im : ',sprintf('% 2.4g',pos(2))],...
               ['Par: ',sprintf('% 2.8g',data(I,1))]};

        if rec == 1
            rec = 0;
            ind = I;
            for j = 1:nLines - 2
                hDatatip{j}.Cursor.Position = [XData{j}(ind) YData{j}(ind) 0];
                drawnow;
            end
            rec = 1;
        end
        
    end
end

%% Build ObjF Functions
function buildObjectiveFunction(ho,~)
h = guidata(ho);
cla(h.ax)
%Check for setup errors
switch h.getError
    case -1
        errordlg('Some information might be unset or wrong. Check for error','Error','modal')
        return
    case -2 
        errordlg('Invalid Model','Error','modal')
        return
    case {-3,-4}
        errordlg('Problem unsetted','Error','modal')
        return
    otherwise
        %Check model setup consistency
        try
            numVar = length(h.num.x0);
            if nargout(h.num.model) == 2
                [numF,MassM] = feval(h.num.model,[],h.num.x0,h.num.par0);
            else
                numF = feval(h.num.model,[],h.num.x0,h.num.par0);
                MassM = eye(numVar);
            end
        catch
            errordlg('Check model! Some inconsistency was identified.')
            return
        end
end

if length(numF) ~= length(h.num.x0)
    err = sprintf('Model returns a vector of length %d, but the length of initial conditions vector is %d',...
                        [length(numF),numVar]);
    errordlg(err,'Error','modal')
    return
end

%Check bound consistency
if length(h.num.lb) ~= length(h.num.par0(:)) || length(h.num.ub) ~= length(h.num.par0(:))
    errordlg('Lower and upper bounds dimensions do not agree','Error','modal')
    return
% elseif any(h.num.lb >= h.num.ub)
%     errordlg('Inconsistent bounds','Error','modal')
%     return
end

%Check SS solvability
testF = h.num.model([],h.num.x0,h.num.par0);
if any(isnan(testF))
    errordlg('The model return NaN with given initial guess. Steady-State cannot be found.')
    return
end
     
%Disable interface
set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
    h.save, h.load, h.reset, h.loadMesh, h.saveMesh, h.buildObjective, ...
    h.nextPlot, h.plotAll, h.exportFigure, h.getRes],'Enable','off')

set(h.resView,'Enable','Inactive')
set(h.pbar2,'Position',[0.7615 .515 .2255*eps .05],'BackgroundColor',[1 0 0])
set(h.pbar,'Position',[0.762 .515 .225*eps .05],'BackgroundColor',[1 0 0])
set(h.stopBuilding,'Visible','on','Enable','on')
drawnow;

if ~h.isMeshLoaded
    h.mesh = [];
    numPar = length(h.num.par0);
    if length(h.num.lb) ~= numPar || length(h.num.ub) ~= numPar
        errordlg('Lower and upper bounds dimensions do not agree','Error','modal')
        return
%     elseif any(h.num.lb > h.num.ub)
%         errordlg('Inconsistent bounds','Error','modal')
%         return
    end

    %Setup
    x_ss = h.num.x0;
    fun = @(state,par) h.num.model([],state,par);
    opt_ss = optimoptions(@fsolve,'algorithm',h.ssList{get(h.ssSolver,'Value')},...
                            'TolFun',h.num.TolFun,'TolX',h.num.TolX,'Jacobian','off','Display','off');

    %Build
    if numPar == 1
        meshRef = h.num.MeshRef(1);
        par = linspace(h.num.lb,h.num.ub,meshRef)';
        f = inf(meshRef,1);
        x = zeros(meshRef,1);
        for i = 1:meshRef

            %OBJ
            [x_ss,f_ss,status,~,J_ss] = fsolve(fun,x_ss,opt_ss,par(i));
            if status <= 0
                [x_ss,f_ss,status,~,J_ss] = hfsolve(fun,x_ss,opt_ss,5,par(i));
                if status <= 0
                    [x_ss,f_ss,status,~,J_ss] = hfsolve(fun,x_ss,opt_ss,25,par(i));
                end
            end
            P = labudde(J_ss,MassM);
            nPol2 = length(P) - 1;
            for j = nPol2-1:-1:1
                for k = 1:nPol2-1
                    pos = (2*j-k+1);
                    if pos > nPol2+1 || pos < 1 
                        Delta(j,k) = 0;
                    else
                        Delta(j,k) = P(pos);
                    end
                end
            end
            con = -min(0,P(end));

            if status <=0
                f(i) = nan;
            else
                f(i) = norm(f_ss) + det(Delta)^2 + 1e6*con;
            end

            %OBJ
            
            x = par;
            if ~mod(i,meshRef/20)
                set(h.pbar2,'Position',[0.7615 .515 .2255*i/meshRef .05],'BackgroundColor',[1 0 0]); drawnow;
                if get(h.stopBuilding,'Value')
                    break
                end
            end
        end
        set(h.pbar2,'Position',[0.7615 .515 .2255 .05],'BackgroundColor',[0 1 0]); drawnow;
        set(h.main,'nextplot','add','CurrentAxes',h.ax)
        plot(x,f,'LineWidth',1.5)
        set(h.main,'nextplot','new')
        set(h.ax,'YGrid','on','Box','on','FontSize',h.fs_num+2)
        ylabel('Objective function','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        xlabel('Parameter','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        h.mesh = [x,f];
    elseif numPar == 2
        meshRefx = h.num.MeshRef(1);
        meshRefy = h.num.MeshRef(2);
        x = linspace(h.num.lb(1),h.num.ub(1),meshRefx);
        y = linspace(h.num.lb(2),h.num.ub(2),meshRefy);
        f = inf(meshRefy,meshRefx);

        for i = 1:meshRefx
            for j = 1:meshRefy
                [x_ss,f_ss,status,~,J_ss] = fsolve(fun,x_ss,opt_ss,[x(i),y(j)]);
                if status <= 0
                    [x_ss,f_ss,status,~,J_ss] = hfsolve(fun,x_ss,opt_ss,5,[x(i),y(j)]);
                    if status <= 0
                        [x_ss,f_ss,status,~,J_ss] = hfsolve(fun,2*x_ss.*rand(numVar,1),opt_ss,25,[x(i),y(j)]);
                    end
                end
                
                P = labudde(J_ss,MassM);
                nPol2 = length(P) - 1;
                for m = nPol2-1:-1:1
                    for k = 1:nPol2-1
                        pos = (2*m-k+1);
                        if pos > nPol2+1 || pos < 1 
                            Delta(m,k) = 0;
                        else
                            Delta(m,k) = P(pos);
                        end
                    end
                end
                con = -min(0,P(end));

                if status <=0
                    f(j,i) = nan;
                else
                    f(j,i) = det(Delta)^2;
                end
            
            end
            if ~mod(i,meshRefx/20)
                set(h.pbar2,'Position',[0.7615 .515 .2255*i/meshRefx .05],'BackgroundColor',[1 0 0]); drawnow;
                if get(h.stopBuilding,'Value')
                    break
                end
            end
        end
        set(h.pbar2,'Position',[0.7615 .515 .2255 .05],'BackgroundColor',[0 1 0]); drawnow;
        [X,Y] = meshgrid(x,y);
        set(h.main,'CurrentAxes',h.ax)
        surfc(X,Y,f,'EdgeColor','none')
        zlabel('Objective function','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        ylabel('Par #2','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        xlabel('Par #1','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
        disp('Finish building OF')
        h.mesh(:,:,1) = X; h.mesh(:,:,2) = Y; h.mesh(:,:,3) = f; 
    else
        errordlg('Cannot handle more tha 2 parameters.','Error','modal')
        return
    end
else
        
end
    
%Enable interface
set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
    h.load, h.reset, h.loadMesh, h.buildObjective, ...
    h.exportFigure, h.saveMesh, h.getRes],'Enable','on')

set(h.stopBuilding,'Visible','off','Enable','off','Value',0)
h.isEigenMap = 0;
h.nSSOF = sum(sum(isnan(f)));
fprintf('\n Not SS evaluations = %2d\n',h.nSSOF);
guidata(ho,h)

end

function saveMeshFunction(ho,~)
h = guidata(ho);

defaultName = [func2str(h.num.model) '_' datestr(now,'yyyy-mm-dd,HH:MM:SS')];
[FileName,PathName,FilterIndex] = uiputfile('*.mat','Save mesh',defaultName);

if FilterIndex == 1
    fullName = [PathName FileName];
    mesh_orbitFinder = h.mesh; %#ok<NASGU>
    save(fullName,'mesh_orbitFinder')
else
    errordlg('Data not saved!','Error','modal')
    return
end

end

function loadMeshFunction(ho,~)
h = guidata(ho);

[FileName,PathName,FilterIndex] = uigetfile('*.mat','Load mesh');

if FilterIndex == 1
    vinfo = who('-file',[PathName FileName]);
    if ~ismember('mesh_orbitFinder',vinfo)
        errordlg('Invalid file!','Error','modal');
        return
    else
        tmp = load([PathName FileName],'mesh_orbitFinder');
        h.mesh = tmp.mesh_orbitFinder;
        set(h.pbar2,'Position',[0.7615 .515 .2255 .05],'BackgroundColor',[0 1 0]); drawnow;
        h.getError = -4;
        h.isMeshLoaded = 1;
        
        if length(size(h.mesh)) == 3
        
            set(h.main,'CurrentAxes',h.ax)
            surfc(h.mesh(:,:,1),h.mesh(:,:,2),h.mesh(:,:,3),'EdgeColor','none')
            zlabel('Objective function','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
            ylabel('Par #2','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
            xlabel('Par #1','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
            set(h.pbar2,'Position',[0.7615 .515 .2255 .05],'BackgroundColor',[0 1 0]);
            set(h.exportFigure,'Enable','on');
            drawnow;

        elseif length(size(h.mesh)) == 2

            set(h.main,'nextplot','add','CurrentAxes',h.ax)
            plot(h.mesh(:,1),h.mesh(:,2),'LineWidth',1.5)
            set(h.main,'nextplot','new')
            set(h.ax,'YGrid','on','Box','on','FontSize',h.fs_num+2)       
            ylabel('Objective function','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
            xlabel('Parameter','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
            set(h.pbar2,'Position',[0.7615 .515 .2255 .05],'BackgroundColor',[0 1 0]);
            set(h.exportFigure,'Enable','on');
            drawnow;

        end
    
        
    end
else
    errordlg('Invalid file!','Error','modal')
    return
end

guidata(ho,h)

end

%% Eigen Map Functions
function eigenMapFunction(ho,~)
h = guidata(ho);

cla(h.ax)
%Check for setup errors
switch h.getError
    case -1
        errordlg('Some information might be unset or wrong. Check for error','Error','modal')
        return
    case -2 
        errordlg('Invalid Model','Error','modal')
        return
    case {-3,-4}
        errordlg('Problem unsetted','Error','modal')
        return
    otherwise
        %Check model setup consistency
        try
            numVar = length(h.num.x0);
            if nargout(h.num.model) == 2
                [numF,MassM] = feval(h.num.model,[],h.num.x0,h.num.par0);
            else
                numF = feval(h.num.model,[],h.num.x0,h.num.par0);
                MassM = eye(numVar);
            end
        catch
            errordlg('Check model! Some inconsistency was identified.')
            return
        end
end

if length(numF) ~= length(h.num.x0)
    err = sprintf('Model returns a vector of length %d, but the length of initial conditions vector is %d'...
                    ,[length(numF),numVar]);
    errordlg(err,'Error','modal')
    return
end

%Check bound consistency
if length(h.num.lb) ~= length(h.num.par0(:)) || length(h.num.ub) ~= length(h.num.par0(:))
    errordlg('Lower and upper bounds dimensions do not agree','Error','modal')
    return
elseif any(h.num.lb >= h.num.ub)
    errordlg('Inconsistent bounds','Error','modal')
    return
end

%Check SS solvability
testF = h.num.model([],h.num.x0,h.num.par0);
if any(isnan(testF))
    errordlg('The model return NaN with given initial guess. Steady-State cannot be found.')
    return
end
     
%Disable interface
set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
    h.save, h.load, h.reset, h.loadMesh, h.saveMesh, h.buildObjective, ...
    h.nextPlot, h.plotAll, h.exportFigure, h.getRes],'Enable','off')

set(h.resView,'Enable','Inactive')
set(h.pbar2,'Position',[0.7615 .515 .2255*eps .05],'BackgroundColor',[1 0 0])
set(h.pbar,'Position',[0.762 .515 .225*eps .05],'BackgroundColor',[1 0 0])
set(h.stopBuilding,'Visible','on','Enable','on')
drawnow;

%Setup
x0 = h.num.x0(:);
fun = @(state,par) h.num.model([],state,par);
opt_ss = optimoptions(@fsolve,'algorithm',h.ssList{get(h.ssSolver,'Value')},...
                        'TolFun',h.num.TolFun,'TolX',h.num.TolX,'Jacobian','off','Display','off');
par0 = h.num.par0(:);

% Dynamical degrees of freedom at initial point
[~,~,~,~,J_ss0] = hfsolve(fun,x0,opt_ss,25,par0);

P0 = labudde(J_ss0,MassM);
nPol = length(P0)-1;
numPar = length(h.num.par0);

%Build
if numPar == 1
    meshRef = h.num.MeshRef(1);
    par = linspace(h.num.lb,h.num.ub,meshRef)';
    flag(meshRef) = 0; lambda(meshRef,nPol) = 0;
    x(meshRef,numVar) = 0;
    for i = 1:meshRef

        [x(i,:),~,flag(i),~,Jss] = fsolve(fun,x0,opt_ss,par(i));
        if flag(i) <= 0
            [x(i,:),~,flag(i),~,Jss] = hfsolve(fun,x0,opt_ss,5,par(i));
            if flag(i) <= 0
                [x(i,:),~,flag(i),~,Jss] = hfsolve(fun,x0,opt_ss,25,par(i));
            end
        end
        x0 = x(i,:);
        lambda(i,:) = roots(labudde(Jss,MassM));

        if ~mod(i,meshRef/20)
            set(h.pbar2,'Position',[0.7615 .515 .2255*i/meshRef .05],'BackgroundColor',[1 0 0]); drawnow;
            if get(h.stopBuilding,'Value')
                break
            end
        end
    end
    r = real(lambda);
    s = imag(lambda);
    
    set(h.pbar2,'Position',[0.7615 .515 .2255 .05],'BackgroundColor',[0 1 0]); drawnow;
    set(h.main,'nextplot','add','CurrentAxes',h.ax)

    plot(r,s,'LineWidth',1,'LineStyle','-','Marker','.','MarkerSize',10);
    line([0,0],h.ax.YLim,'color','k','linestyle',':')
    line(h.ax.XLim,[0,0],'color','k','linestyle',':')

    set(h.main,'nextplot','new')
    set(h.ax,'YGrid','on','Box','on','FontSize',h.fs_num+2)
    ylabel('Im(\lambda)','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
    xlabel('Re(\lambda)','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
    h.eigenMap = [par, lambda];
else
    errordlg('Cannot handle more than 1 parameter.','Error','modal')
end

    
%Enable interface
set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
    h.load, h.reset, h.loadMesh, h.buildObjective, ...
    h.exportFigure, h.saveMesh, h.getRes],'Enable','on')

set(h.stopBuilding,'Visible','off','Enable','off','Value',0)
h.isEigenMap = 1;
h.nSSEM = sum(flag<=0);
fprintf('\n Not SS evaluations = %2d\n',h.nSSEM);
guidata(ho,h)
end


%% One parameter continuation

function oneParCont(ho,~)

h = guidata(ho);

%Check for setup errors
switch h.getError
    case -1
        errordlg('Some information might be unset or wrong. Check for error','Error','modal')
        return
    case -2 
        errordlg('Invalid Model','Error','modal')
        return
    case {-3,-4}
        errordlg('Problem unsetted','Error','modal')
        return
    otherwise
        %Check model setup consistency
        try
            numVar = length(h.num.x0);
            if nargout(h.num.model) == 2
                [numF,MassM] = feval(h.num.model,[],h.num.x0,h.num.par0);
            else
                numF = feval(h.num.model,[],h.num.x0,h.num.par0);
                MassM = eye(numVar);
            end
        catch
            errordlg('Check model! Some inconsistency was identified.')
            return
        end
end

if length(numF) ~= length(h.num.x0)
    err = sprintf('Model returns a vector of length %d, but the length of initial conditions vector is %d'...
                    ,[length(numF),numVar]);
    errordlg(err,'Error','modal')
    return
end

%Check bound consistency
if length(h.num.lb) ~= length(h.num.par0(:)) || length(h.num.ub) ~= length(h.num.par0(:))
    errordlg('Lower and upper bounds dimensions do not agree','Error','modal')
    return
% elseif any(h.num.lb >= h.num.ub)
%     errordlg('Inconsistent bounds','Error','modal')
%     return
end

%Check SS solvability
testF = h.num.model([],h.num.x0,h.num.par0);
if any(isnan(testF))
    errordlg('The model return NaN with given initial guess. Steady-State cannot be found.')
    return
end
     
%Disable interface
set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
    h.save, h.load, h.reset, h.loadMesh, h.saveMesh, h.buildObjective, ...
    h.nextPlot, h.plotAll, h.exportFigure, h.getRes],'Enable','off')

set(h.resView,'Enable','Inactive')
set(h.pbar2,'Position',[0.7615 .515 .2255*eps .05],'BackgroundColor',[1 0 0])
set(h.pbar,'Position',[0.762 .515 .225*eps .05],'BackgroundColor',[1 0 0])
set(h.stopBuilding,'Visible','on','Enable','on')
drawnow;

%Setup
h.isEigenMap = 0;
x0 = h.num.x0(:);
fun = @(state,par) h.num.model([],state,par);
opt_ss = optimoptions(@fsolve,'algorithm',h.ssList{get(h.ssSolver,'Value')},...
                        'TolFun',h.num.TolFun,'TolX',h.num.TolX,'Jacobian','off','Display','off');
par0 = h.num.par0(:);

% Dynamical degrees of freedom at initial point
[~,~,~,~,J_ss0] = hfsolve(fun,x0,opt_ss,25,par0);

P0 = labudde(J_ss0,MassM);
nPol = length(P0)-1;
numPar = length(h.num.par0);

%Build
if numPar == 1
    meshRef = h.num.MeshRef(1);
    h.stateSel = 1;
    ds = (h.num.ub - h.num.lb)/meshRef;
    flag = zeros(meshRef,1); lambda = zeros(meshRef,nPol);
    x = zeros(meshRef,numVar);
    par = zeros(meshRef,1);
    [x0,fval0,flag0,~,Jss] = fsolve(fun,x0,opt_ss,par0);
    if flag0 <= 0
        [x0,~,flag0,~,Jss] = hfsolve(fun,x0,opt_ss,5,par0);
        if flag0 <= 0
            [x0,~,~,~,Jss] = hfsolve(fun,x0,opt_ss,25,par0);
        end
    end
    delta = 1e-6;
    b = [zeros(numVar,1); 1];
    for i = 1:meshRef
        parp = par0 + delta;
        dFdp = (fun(x0,parp) - fval0)/delta;
        A = [Jss dFdp; null([Jss dFdp])'];
        TgV = linsolve(A,b);
        xPred = x0 + TgV(1:end-1)*ds;
        parPred = par0 + abs(TgV(end))*ds;
%         parPred = par0 + ds;
        [x(i,:),fval0,flag(i),~,Jss] = fsolve(fun,xPred,opt_ss,parPred);
        if flag(i) <= 0
            [x(i,:),fval0,flag(i),~,Jss] = hfsolve(fun,xPred,opt_ss,5,parPred);
            if flag(i) <= 0
                [x(i,:),fval0,flag(i),~,Jss] = hfsolve(fun,xPred,opt_ss,25,parPred);
            end
        end
        x0 = x(i,:)';
        lambda(i,:) = roots(labudde(Jss,MassM));
        par(i) = parPred;
        par0 = parPred;
        if ~mod(i,meshRef/20)
            set(h.pbar2,'Position',[0.7615 .515 .2255*i/meshRef .05],'BackgroundColor',[1 0 0]); drawnow;
            if get(h.stopBuilding,'Value')
                break
            end
        end
    end
    r = real(lambda);

    h.stableC = nan(meshRef,numVar);
    h.unstableC = h.stableC;
    for i = 1:meshRef
        if any(sign(r(i,:))==1)
            h.unstableC(i,:) = x(i,:);
        else
            h.stableC(i,:) = x(i,:);
        end
    end
    
    h.par = par;
    set(h.pbar2,'Position',[0.7615 .515 .2255 .05],'BackgroundColor',[0 1 0]); drawnow;
    set(h.main,'nextplot','add','CurrentAxes',h.ax)
    
    pp = plot(par,h.stableC(:,h.stateSel),par,h.unstableC(:,h.stateSel));
    set(pp(1),'LineStyle','-','Color','k','LineWidth',1.5)
    set(pp(2),'LineStyle',':','Color','r','LineWidth',1)
    set(h.main,'nextplot','new')
    set(h.ax,'YGrid','on','Box','on','FontSize',h.fs_num+2)
    ylabel('State #1','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
    xlabel('Par #1','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
    h.isBranchDiagram = 1;
    h.numVar = numVar;
        
    h.nSSCont = sum(flag<=0);
    fprintf('\n Not SS evaluations = %2d\n',h.nSSCont);
else
    errordlg('Cannot handle more than 1 parameter.','Error','modal')
end

    
%Enable interface
set([h.edit_model, h.selModel_btn ,h.pushSwarm, h.pushHybrid, h.pushFmin, h.pushLM, h.simulate, ...
    h.load, h.reset, h.loadMesh, h.buildObjective, ...
    h.exportFigure, h.saveMesh, h.getRes, h.nextPlot],'Enable','on')

set(h.stopBuilding,'Visible','off','Enable','off','Value',0)
guidata(ho,h)



end

%% Lesser Callback functions
function editme(~,~)
open(mfilename)
end
function getResFunction(ho,~)
h = guidata(ho);

h.num.par0 = h.num.sol.par;

if ~isfield(h.num.sol,'x_T')
    h.num.x0 = h.num.sol.x;
else
    h.num.x0 = h.num.sol.x_T;
end

str = ['[' num2str(h.num.x0(:)','%.4f, ') ]; str(end) = ']';
set(h.x0,'String',str,'FontAngle','normal','ForegroundColor',[1 0 0])

str = ['[' num2str(h.num.par0(:)','%.4f, ') ]; str(end) = ']';
set(h.par0,'String',str,'FontAngle','normal','ForegroundColor',[1 0 0])

if isfield(h.num.sol,'T')
    h.num.it = [0 h.num.sol.T]';
    str = ['[' num2str(h.num.it','%.4f, ') ]; str(end) = ']';
    set(h.it,'String',str,'FontAngle','normal','ForegroundColor',[1 0 0])
end

guidata(ho,h)
end

function saveResultsFunction(ho,~)
h = guidata(ho);

defaultName = [func2str(h.num.model) '_results_' datestr(now,'yyyy-mm-dd,HH:MM:SS')];
[FileName,PathName,FilterIndex] = uiputfile('*.mat','Save mesh',defaultName);

if FilterIndex == 1
    fullName = [PathName FileName];
    oF_data = h.num; %#ok<NASGU>
    save(fullName,'oF_data')
else
    errordlg('Data not saved!','Error','modal')
    return
end

end

function loadResultsFunction(ho,~)
h = guidata(ho);
cla(h.ax)
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Load results');

if FilterIndex == 1
    vinfo = who('-file',[PathName FileName]);
    if ~ismember('oF_data',vinfo)
        errordlg('Invalid file!','Error','modal');
        return
    else
        tmp = load([PathName FileName],'oF_data');
        h.num = tmp.oF_data;
        sol = h.num.sol;
        
        if h.idScreen == 1920
            s_final{1} = sprintf('\n*** Data Loaded & Results ***            ');
            s_final{2} = sprintf('\nFile ............................................. %s',[func2str(h.num.model) '.m']);
            s_final{3} = sprintf(  'Fval ............................................. %.6e',sol.f);

            s = repmat('%.4e, ',1,length(sol.par)); s(end-1:end) = [];
            s_final{4} = sprintf([ 'Par .............................................. [' s ']'], sol.par);
            s = repmat('%.4e, ',1,length(sol.x)); s(end-1:end) = [];

            s_final{5} = sprintf([ 'Ss ............................................... [' s ']'], sol.x);
            s_final{6} = sprintf(  'Time ............................................. %.4f s',sol.time);
            s_final{7} = sprintf(  'FCount ........................................... %d',sol.nS);
            s_final{8} = sprintf(  'Derivative ....................................... %.4f',sol.derivative(1));
            s_final{9} = sprintf(  'NotSS ............................................ %.1f',sol.notSS);
            s_final{10} = sprintf( 'OptType .......................................... %s',sol.OptType);
            if isfield(sol,'T')
                s_final{11} = sprintf( 'Period ........................................... %.4f',sol.T);
                s = repmat('%.4e, ',1,length(sol.x_T)); s(end-1:end) = [];
                s_final{12} = sprintf([ 'Pp ............................................... [' s ']\n'], sol.x_T);
            end
            
        else
        
            s_final{1} = sprintf('*** Data Loaded & Results ***      ');
            s_final{2} = sprintf( '\n File ................................. %s',[func2str(h.num.model) '.m']);
            s_final{3} = sprintf( 'Fval ................................. %.6e',sol.f);

            s = repmat('%.4e, ',1,length(sol.par)); s(end-1:end) = [];
            s_final{4} = sprintf(['Par .................................. [' s ']'], sol.par);
            s = repmat('%.4e, ',1,length(sol.x)); s(end-1:end) = [];

            s_final{5} = sprintf(['Ss ................................... [' s ']'], sol.x);
            s_final{6} = sprintf( 'Time ................................. %.4f s',sol.time);
            s_final{7} = sprintf( 'FCount ............................... %d',sol.nS);
            s_final{8} = sprintf( 'Derivative ........................... %.4f',sol.derivative(1));
            s_final{9} = sprintf( 'NotSS ................................ %.1f',sol.notSS);
            s_final{10} = sprintf('OptType .............................. %s',sol.OptType);
            if isfield(sol,'T')
                s_final{11} = sprintf( 'Period .............................. %.4f',sol.T);
                s = repmat('%.4e, ',1,length(sol.x_T)); s(end-1:end) = [];
                s_final{12} = sprintf([ 'Pp .................................. [' s ']\n'], sol.x_T);
            end
        end
        
        set(h.resView,'String',s_final,'HorizontalAlignment','right','FontSize',h.fs_num,'Enable','on')
        
        
        set(h.selModel,'String',which(func2str(h.num.model)),'FontAngle','normal','ForegroundColor',[0 0 0])
        
        str = ['[' num2str(h.num.x0(:)','%.4f, ') ]; str(end) = ']';
        set(h.x0,'String',str,'FontAngle','normal','ForegroundColor',[0 0 0])
        
        str = ['[' num2str(h.num.par0(:)','%.4f, ') ]; str(end) = ']';
        set(h.par0,'String',str,'FontAngle','normal','ForegroundColor',[0 0 0])
        
        str = ['[' num2str(h.num.lb(:)','%.0f, ') ]; str(end) = ']';
        set(h.lb,'String',str,'FontAngle','normal','ForegroundColor',[0 0 0])
        
        str = ['[' num2str(h.num.ub(:)','%.0f, ') ]; str(end) = ']';
        set(h.ub,'String',str,'FontAngle','normal','ForegroundColor',[0 0 0])
        
        str = ['[' num2str(h.num.it(:)','%.0f, ') ]; str(end) = ']';
        set(h.it,'String',str,'FontAngle','normal','ForegroundColor',[0 0 0])
        
        set(h.TolFun,'String',h.num.TolFun,'FontAngle','normal','ForegroundColor',[0 0 0])
        set(h.TolX,'String',h.num.TolX,'FontAngle','normal','ForegroundColor',[0 0 0])
        set(h.ObjectiveLimit,'String',h.num.ObjectiveLimit,'FontAngle','normal','ForegroundColor',[0 0 0])
        set(h.MaxFunEvals,'String',h.num.MaxFunEvals,'FontAngle','normal','ForegroundColor',[0 0 0])
        set(h.MaxTrials,'String',h.num.MaxTrials,'FontAngle','normal','ForegroundColor',[0 0 0])
        set(h.nger,'String',h.num.nger,'FontAngle','normal','ForegroundColor',[0 0 0])
        set(h.npas,'String',h.num.npas,'FontAngle','normal','ForegroundColor',[0 0 0])
        set(h.AbsTol,'String',h.num.AbsTol,'FontAngle','normal','ForegroundColor',[0 0 0])
        set(h.RelTol,'String',h.num.RelTol,'FontAngle','normal','ForegroundColor',[0 0 0])
        
        if h.num.isSimulated
            
            %Plot first result
            set(h.main,'nextplot','add','CurrentAxes',h.ax);
            h1 = plot(h.num.t,h.num.y(:,1));
            set(h.main,'nextplot','new')
            set(h1,'LineWidth',1.5)
            %h.ax.YAxis.TickLabelFormat = '%,.4f';
            set(h.ax,'YGrid','on','Box','on','FontSize',h.fs_num+2)            

            ylabel('State #1','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)
            xlabel('Time','FontWeight','bold','FontName','Palatino','fontsize',h.fs_num+2)

            h.numVar = length(h.num.x0);
            h.plotCounterI = 2;
            h.plotCounterJ = 0;
            set([h.load, h.save, h.nextPlot, h.plotAll, h.exportFigure],'Enable','on')
           
        end
        
    end
else
    errordlg('Invalid file!','Error','modal')
    return
end
h.getError = 0;
guidata(ho,h)
end

function openExample(~,~)

if exist('orbitfinder_modelExample.m','file') == 2
    edit('orbitfinder_modelExample.m')
end

end

function editFunction(ho,~)
h = guidata(ho);

if isempty(h.num.model)
    [FileName,PathName,FilterIndex] = uiputfile({'*.m','MATALB Code files (*.m)'},'Edit a new model');
    if FilterIndex == 1
        edit([PathName FileName]);
    end
else
    edit(func2str(h.num.model))
end

end

function getPositiveNumber(ho,~)
h = guidata(ho);
t = get(ho,'String');
t(t==',') = '.';

t = str2num(t);
if any(isnan(t)) || t < 0 || isinf(t) || any(isempty(t))
    set(ho, 'String', 'error','ForegroundColor',[1 0 0],'FontAngle','italic');
    errordlg('Input must be a real positive number','Error','modal');
    h.getError = -1;
elseif strcmp(ho.Tag,'MaxTrials') && t < 1
    h.num.MaxTrials = 1;
    set(ho,'FontAngle','normal','ForegroundColor',[0 0 0],'String',1);
    h.getError = 0;
else
    h.num.(ho.Tag) = t;
    set(ho,'FontAngle','normal','ForegroundColor',[0 0 0],'String',t);
    h.getError = 0;
end
guidata(ho,h)
end

function getVector(ho,~)
h = guidata(ho);
vec =  str2num(get(ho, 'String'));
[m,n] = size(vec);

if isempty(vec) || any(isnan(vec)) ||  (m > 1 && n > 1) || any(isinf(vec))
    set(ho, 'String', 'error','ForegroundColor',[1 0 0],'FontAngle','italic');
    errordlg('Input must be a real scalar or a vector of real numbers','Error','modal')
    h.getError = -1;
else
    h.num.(ho.Tag) = vec;
    set(ho,'FontAngle','normal','ForegroundColor',[0 0 0]);
    h.getError = 0;
end
guidata(ho,h)
end

function selBtnPress(ho,~)
h = guidata(ho);
[FileName,PathName,FilterIndex] = uigetfile({'*.m','MATLAB Code files (*.m)'},'Select the Model file...');
if FilterIndex == 0 || FilterIndex > 1
    errordlg('Please select a valid *.m file','Error','modal')
    set(h.selModel,'fontangle','italic','foregroundcolor',[1 0 0],'String','error')

    h.getError = -2;
else
    addpath(PathName);
    h.num.model = str2func(FileName(1:end-2));
    set(h.selModel,'fontangle','normal','foregroundcolor',[0 0 0],'String',[PathName FileName])
    h.getError = 0;
end
guidata(ho,h)
end

function resetBtnPress(~,~)
feval(str2func(mfilename));
end

function selModelFunction(ho,~)
h = guidata(ho);

model = get(h.selModel,'String');
if exist(model,'file')
    h.num.model = str2func(model);
    set(h.selModel,'fontangle','normal','foregroundcolor',[0 0 0])
    h.getError = 0;
else
    errordlg('Invalid function! Check if it is in path.','Error','modal')
    set(h.selModel,'fontangle','italic','foregroundcolor',[1 0 0],'String','error')
    h.getError = -2; 
    return
end

guidata(ho,h)
end
%% Other Functions
function [xo,Ot,nS] = mySwarm(S,x0,Lb,Ub,nger,npas,ObjectiveLimit,MaxFunEvals,...
                                resView,stopSwarm,stopHybrid,s_initial,idScreen)

n = length(x0);
problem = -1;
c1 = 1;
c2 = 1;
viewIndex = floor(idScreen/220);

w = 0.9;
tw = (w-4e-3)/nger;
Ot = feval(S,x0)*problem;
y(1) = Ot;
p(n,npas) = 0;
v(n,npas) = 0;
x(n,npas) = 0;

del = Ub-Lb;
 
for j = 1:npas
    x(:,j) = Lb + rand(n,1).*del;
    v(:,j) = 0.1*rand(n,1).*del;
    p(:,j) = x(:,j);
end	

x(:,1) = x0;
p(:,1) = x0;
nS = 1;
ipg = 1;
 
s_initial{2} = sprintf('\n F-Count               Best value');
set(resView,'String',s_initial,'HorizontalAlignment','center'); drawnow;

y(npas) = 0;
for j = 2:npas
    y(j) = feval(S,x(:,j))*problem;
    nS = nS + 1;
    if y(j) > Ot
        Ot = y(j);
        ipg = j;
   end
end

s{1} = sprintf(' %7.0d %24.6e',[nS,abs(Ot)]);
set(resView,'String',[s_initial s{1}]); drawnow;
xo = p(:,ipg);

s{nger} = [];
for ig = 1:nger
    rnd1 = ones(n,1)*rand(1,npas);
    rnd2 = ones(n,1)*rand(1,npas);
    v = w*v + c1*rnd1.*(p-x) + c2*rnd2.*( p(:,ipg)*ones(1,npas) - x );
    x = x + v;
    
    for i = 1:n
        j = find(x(i,:) > Ub(i));
        if ~isempty(j)
            x(i,j) = Ub(i);
            v(i,j) = 0;
        end

        j = find(x(i,:) < Lb(i));
        if ~isempty(j)
            x(i,j) = Lb(i);
            v(i,j) = 0;
        end
    end

    for j = 1:npas
        val = feval(S,x(:,j))*problem;
        nS = nS + 1;

        if val > y(j)
            y(j) = val;
            p(:,j) = x(:,j);
            if val > Ot
                Ot = val;
                ipg = j;
            end
        end
    end

    s{ig+1} = sprintf(' %7.0d %24.6e',[nS,abs(Ot)]);
    if ig > viewIndex
        set(resView,'String',[s_initial s{ig-viewIndex:ig}]); drawnow;
    else
        set(resView,'String',[s_initial s{1:ig+1}]); drawnow;
    end
    
    xo = p(:,ipg);

    w = w - tw;        
    
    if get(stopSwarm,'Value') || get(stopHybrid,'Value') || Ot*problem <= ObjectiveLimit || nS > MaxFunEvals
        break
    end
end

Ot = Ot*problem;

end
