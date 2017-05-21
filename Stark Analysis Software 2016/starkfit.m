function starkfit()
%GUI for Stark spectra analysis

%% FIGURE INITIALIZATION  ===========================================================
%set figure size to fit any screen
ssize = get(0,'ScreenSize');
ssize = ssize+[9 38 -20 -68];

figMain = figure('Name','Stark Analysis','Units','pixels',...
    'Position',ssize,'Color',[1 1 1]*0.95,... 
    'MenuBar','none','ToolBar','none','NumberTitle','off',...
    'KeyReleaseFcn',{@keyShortcuts});
fsize = get(figMain,'Position');

%% TOOLBAR ==========================================================================
%Toolbar panel & buttons ------------------------------------------------------------
    pnlToolbar = uipanel('Parent',figMain,'Units','pixels',...
        'Position',[0 fsize(4)-32 650 32],'BorderType','none'); %[0 669 650 32]
% -----------------------------------------------------------------------------------
    cmdSave = uicontrol('Parent',pnlToolbar,'Position',[10 0 50 28],...
        'String','SAVE','FontWeight','bold','Callback',{@connectGUI});
    cmdLoad = uicontrol('Parent',pnlToolbar,'Position',[60 0 50 28],...
        'String','LOAD','FontWeight','bold','Callback',{@connectGUI});
% -----------------------------------------------------------------------------------
    toolFiles = uicontrol('Parent',pnlToolbar,'Position',[130 0 50 28],...
        'String','Files','Callback',{@changePanel});
    toolFitting = uicontrol('Parent',pnlToolbar,'Position',[180 0 50 28],...
        'String','Fitting','Callback',{@changePanel});
    toolSettings = uicontrol('Parent',pnlToolbar,'Position',[230 0 50 28],...
        'String','Options','Callback',{@changePanel});
    toolMonteCarlo = uicontrol('Parent',pnlToolbar,'Position',[280 0 70 28],...
        'String','MonteCarlo','Callback',{@changePanel});
% -----------------------------------------------------------------------------------
    lblWeight = uicontrol('Parent',pnlToolbar,'Position',[400 15 250 15],...
        'Style','text','String','Absorption Weighting: 1','FontWeight','bold',...
        'HorizontalAlignment','left');
    lblStatus = uicontrol('Parent',pnlToolbar,'Position',[400 0 250 15],...
        'Style','text','String','Fitting Status:','FontWeight','bold',...
        'ForegroundColor','r','HorizontalAlignment','left');

%% PANEL INITIALIZATION =============================================================
%Command panels ---------------------------------------------------------------------
psize = [10 fsize(4)-350 640 310];
pnlSpectra = uipanel('Parent',figMain,'Title','Spectra Management',...
    'Units','pixels','Position',psize,'Visible','on');
pnlParam = uipanel('Parent',figMain,'Title','Fitting Parameters',...
    'Units','pixels','Position',psize,'Visible','off');
pnlSettings = uipanel('Parent',figMain,'Title','Settings & Options',...
    'Units','pixels','Position',psize,'Visible','off');
pnlPlotting = uipanel('Parent',figMain,'Title','Plotting Options',...
    'Units','pixels','Position',psize,'Visible','off');
pnlMonteCarlo = uipanel('Parent',figMain,'Title','Monte Carlo Analysis',...
    'Units','pixels','Position',psize,'Visible','off');
pnlComments = uipanel('Parent',figMain,'Title','Comments',...
    'Units','pixels','Position',psize,'Visible','off');
%Plotting Panels --------------------------------------------------------------------
pnlPlotMain = uipanel('Parent',figMain,'Units','pixels','Position',[0 0 650 fsize(4)-350],... 
    'BorderType','none');
pnlPlotSide = uipanel('Parent',figMain,'Units','pixels','Position',[650 0 fsize(3)-650 fsize(4)],...
    'BorderType','none');

%% SPECTRA MANAGEMENT ===============================================================
%Stark spectra table ----------------------------------------------------------------
    names = {'Filename','Pol','Inc','F(Vrms)','Conc(mM)','<html>Path(&mu;m)</html>'};
    tblStkSpec = uitable('Parent',pnlSpectra,'Position',[10 85 605 200],...
        'ColumnName',names,'ColumnWidth',{200,60,60,70,75,75},...
        'ColumnEditable',logical([0 1 1 1 1 1]),'Data',[],...
        'CellEditCallback',{@connectGUI},'CellSelectionCallback',{@changeSel});
%Add/Delete Stark spectra -----------------------------------------------------------
    uicontrol('Parent',pnlSpectra,'Position',[470 55 135 20],...
        'Style','text','String','Stark Spectra =>','HorizontalAlignment','left',...
        'FontWeight','bold');
    cmdAddStark = uicontrol('Parent',pnlSpectra,'Position',[565 55 25 25],...
        'String','+','FontWeight','bold','FontSize',12,'Callback',{@connectGUI});
    cmdDelStark = uicontrol('Parent',pnlSpectra,'Position',[590 55 25 25],...
        'String','-','FontWeight','bold','FontSize',12,'Callback',{@connectGUI});
%Absorption files -------------------------------------------------------------------
    cmdAddAbs = uicontrol('Parent',pnlSpectra,'Position',[10 60 100 20],...
        'String','Load Absorption','Callback',{@connectGUI});
    lblAbs = uicontrol('Parent',pnlSpectra,'Position',[120 57 255 20],...
        'Style','text','String','','HorizontalAlignment','left');
    cmdAddRef = uicontrol('Parent',pnlSpectra,'Position',[10 40 100 20],...
        'String','Load Reference','Callback',{@connectGUI});
    lblRef = uicontrol('Parent',pnlSpectra,'Position',[120 37 255 20],...
        'Style','text','String','','HorizontalAlignment','left');
%Absorption parameters --------------------------------------------------------------
    uicontrol('Parent',pnlSpectra,'Position',[10 8 135 20],...
        'Style','text','String','Absorption Parameters:','HorizontalAlignment','left',...
        'FontWeight','bold');
    uicontrol('Parent',pnlSpectra,'Position',[150 8 100 20],...
        'Style','text','String','Concentration (mM):','HorizontalAlignment','left');
    txtAbsConc = uicontrol('Parent',pnlSpectra,'Position',[255 10 30 20],...
        'Style','edit','String','1','HorizontalAlignment','left',...
        'BackgroundColor','w','Callback',{@connectGUI});
    uicontrol('Parent',pnlSpectra,'Position',[300 8 100 20],...
        'Style','text','String','Pathlength (um):','HorizontalAlignment','left');
    txtAbsPath = uicontrol('Parent',pnlSpectra,'Position',[385 10 30 20],...
        'Style','edit','String','55','HorizontalAlignment','left',...
        'BackgroundColor','w','Callback',{@connectGUI});
    uicontrol('Parent',pnlSpectra,'Position',[430 8 100 20],...
        'Style','text','String','Incidence:','HorizontalAlignment','left');
    txtAbsInc = uicontrol('Parent',pnlSpectra,'Position',[485 10 30 20],...
        'Style','edit','String','0','HorizontalAlignment','left',...
        'BackgroundColor','w','Callback',{@connectGUI}); 

%% PARAMETERS PANEL =================================================================    
%Model Panel ------------------------------------------------------------------------
    bgpFit = uibuttongroup('Parent',pnlParam,'Units','pixels','Position',[320 250 285 40],...
        'BorderType','none','SelectionChangeFcn',{@connectGUI});
    radAbs = uicontrol('Parent',bgpFit,'Position',[5 5 100 20],...
        'Style','radiobutton','String','Absorption');
    radSim = uicontrol('Parent',bgpFit,'Position',[90 5 100 20],...
        'Style','radiobutton','String','Simultaneous');
    cmdFit = uicontrol('Parent',bgpFit,'Position',[188 2 50 28],...
        'String','FIT','FontWeight','bold','Callback',{@connectGUI});    
%Gaussian parameters ----------------------------------------------------------------
    names = {'Amplitude','Center','Width','Band'};
    tblGauParms = uitable('Parent',pnlParam,'Position',[10 120 590 120],...
        'ColumnName',names,'ColumnWidth',{100},...
        'ColumnEditable',true(1,5),'Data',[ones(1,4)],...
        'CellEditCallback',{@connectGUI},'CellSelectionCallback',{@changeSel});
%Stark parameters -------------------------------------------------------------------
    names = {'<html>Tr&Delta;&alpha;</html>',...
        '<html>m&#8226;&Delta;&alpha;&#8226;m</html>',...
        '<html>&Delta;&mu;</html>','<html>&zeta;</html>'};
    tblStarkParam = uitable('Parent',pnlParam,'Position',[10 10 590 100],...
        'ColumnName',names,'ColumnWidth',{100},...
        'ColumnEditable',true(1,4),'Data',[ones(1,4)],...
        'CellEditCallback',{@connectGUI},'CellSelectionCallback',{@changeSel});
    cmdAddTransition = uicontrol('Parent',pnlParam,'Position',[605 85 25 25],...
        'String','+','FontWeight','bold','FontSize',12,'Callback',{@connectGUI});
    cmdDelTransition = uicontrol('Parent',pnlParam,'Position',[605 60 25 25],...
        'String','-','FontWeight','bold','FontSize',12,'Callback',{@connectGUI});
    cmdError = uicontrol('Parent',pnlParam,'Position',[605 10 25 25],...
        'String','<html>&sigma;</html>','FontSize',14,'FontWeight','bold',...
        'Callback',{@connectGUI});
%Achi -------------------------------------------------------------------------------
    cmdAchi = uicontrol('Parent',pnlParam,'Position',[572 251 28 28],...
        'Style','togglebutton','String','<html>A&chi;</html>','FontWeight','bold',...
        'Callback',{@connectGUI});
    cmdResetAchi = uicontrol('Parent',pnlParam,'Position',[552 220 50 25],...
        'String','Reset','Visible','off','Callback',{@connectGUI});
    tblAchi = uitable('Parent',pnlParam,'Position',[10 10 590 200],...
        'Data',[],'Visible','off',...
        'CellEditCallback',{@connectGUI},'CellSelectionCallback',{@changeSel});

%% SETTINGS PANEL ===================================================================
%Options ----------------------------------------------------------------------------
    names = {'Weight','Gau/Prog','MaxEval','Abs Tol','Sim Tol','TolX'};
    tblOpts = uitable('Parent',pnlSettings,'Position',[10 150 182 130],...
        'RowName',names,'ColumnName','Options','ColumnEditable',true(1,1),...
        'ColumnWidth',{92},'Data',[],...
        'CellEditCallback',{@connectGUI});
%baseline correction & solvent ------------------------------------------------------
    names = {'8M LiCl','Buffer/Glycerol','Ethanol','Butanol','2Me-THF','Toluene'};
    uicontrol('Parent',pnlSettings,'Position',[200 255 50 25],...
        'Style','text','String','Solvent Ref. Ind.:',...
        'HorizontalAlignment','left');
    cboSolvent = uicontrol('Parent',pnlSettings,'Position',[260 255 70 25],...
        'Style','edit','String','1.4','BackgroundColor','w',...
        'Callback',{@connectGUI});
    cmdBaseline = uicontrol('Parent',pnlSettings,'Position',[200 220 55 28],...
        'String','Abs. BL','Callback',{@connectGUI});
    txtBaseline = uicontrol('Parent',pnlSettings,'Position',[260 222 70 25],...
        'Style','edit','String',num2str([600 650]),'BackgroundColor','w',...
        'HorizontalAlignment','left','Callback',{@connectGUI});
    cmdLimits = uicontrol('Parent',pnlSettings,'Position',[200 190 55 28],...
        'String','Fit Range','Callback',{@connectGUI});
    txtLimits = uicontrol('Parent',pnlSettings,'Position',[260 192 70 25],...
        'Style','edit','String',num2str([600 650]),'BackgroundColor','w',...
        'HorizontalAlignment','left','Callback',{@connectGUI});  
%Stark baseline ---------------------------------------------------------------------
    uicontrol('Parent',pnlSettings,'Position',[380 250 200 20],...
        'Style','text','String','Stark Spectra Baseline Parameters:',...
        'HorizontalAlignment','left','FontWeight','bold');
    tblStkBaseline = uitable('Parent',pnlSettings,'Position',[380 10 240 240],...
        'Data',[],'ColumnName',{'Min' 'Max'},'ColumnWidth',{80},...
        'ColumnEditable',true(1,2),'CellEditCallback',{@connectGUI});
    
%% MONTE CARLO PANEL ================================================================
    cmdMonteCarlo = uicontrol('Parent',pnlMonteCarlo,'Position',[10 250 105 25],...
        'String','Run Monte Carlo','Callback',{@connectGUI});
    lblMonteCarlo = uicontrol('Parent',pnlMonteCarlo,'Position',[275 10 350 265],...
        'Style','text','String','Iteration','BackgroundColor','w','HorizontalAlignment','left');
    uicontrol('Parent',pnlMonteCarlo,'Position',[10 215 100 20],...
        'Style','text','String','Iterations:','HorizontalAlignment','left');
    txtMCiter = uicontrol('Parent',pnlMonteCarlo,'Position',[65 220 50 20],...
        'Style','edit','String','50','BackgroundColor','w','Callback',{@connectGUI});
    cmdPlotMC = uicontrol('Parent',pnlMonteCarlo,'Position',[10 180 105 25],...
        'String','Plot Results','Callback',{@connectGUI});
    cmdExportResults = uicontrol('Parent',pnlMonteCarlo,'Position',[10 10 105 25],...
        'String','Export Results','FontWeight','bold','Callback',{@connectGUI});
      
%% PLOTTING PANELS ==================================================================
axMain = axes('Parent',pnlPlotMain,'Position',[0.13,0.11,0.86,0.815],...
    'NextPlot','add','Box','on');
axStk = subplot(6,1,1:3,'parent',pnlPlotSide,...
    'NextPlot','add','Box','on','XTickLabel',[]);
axAbs = subplot(6,1,4:5,'parent',pnlPlotSide,...
    'NextPlot','add','Box','on');%,'XTickLabel',[]);
axRes = subplot(6,1,6,'parent',pnlPlotSide,...
    'NextPlot','add','Box','on');

%% CALLBACKS ========================================================================
    function keyShortcuts(~,evdat)
        if strcmp(evdat.Key,'escape') == 1
            close(figMain);
        end
    end %press ESC to close program

    function changePanel(src,~) %navigation between comman panels
        set([pnlSpectra,pnlParam,pnlSettings,pnlPlotting,pnlMonteCarlo,pnlComments],'Visible','off');
        switch src
            case toolFiles
                set(pnlSpectra,'Visible','on');
            case toolFitting
                set(pnlParam,'Visible','on');
            case toolSettings
                set(pnlSettings,'Visible','on');
            case toolMonteCarlo
                set(pnlMonteCarlo,'Visible','on');
        end
    end

    function changeSel(src,evdat) %selected table cells
        if numel(evdat.Indices) == 0, return; end
        switch src
            case tblGauParms
                vibIndex = evdat.Indices(1);
            case tblStarkParam
                parmIndex = evdat.Indices(1);
            case tblStkSpec
                specIndex = evdat.Indices(1);
        end
    end

    function connectGUI(src,evdat) %all button functionality!!!
        switch src
            case cmdLoad
                [fname,pname] = uigetfile('.mat','Select Stark Project');
                if fname == 0, return; end
                D = load([pname,fname]);
                obj = D.obj;
                obj.ax = [axStk,axAbs,axRes,axMain];
                [p,f,~] = fileparts([pname,fname]);
                projname = [p,filesep,f];
                populateGUI;
            case cmdSave
                [sname,pname] = uiputfile('.mat','Save project as...');
                if sname == 0, return; end
                save([pname,sname],'obj');
                [p,f,~] = fileparts([pname,sname]);
                projname = [p,filesep,f];
            case tblStkSpec
                dat = get(src,'Data'); dat(:,1) = [];
                for j = 1:size(dat,1)
                    obj.sobj(j).setPol(dat{j,1});
                    obj.sobj(j).setInc(dat{j,2});
                    obj.sobj(j).setField(dat{j,3});
                    obj.sobj(j).setConc(dat{j,4});
                    obj.sobj(j).setPath(dat{j,5});
                end
                populateGUI;
            case cmdAddStark
                sobj = stkSpecObj;
                sobj.loadSpec
                obj.sobj = [obj.sobj,sobj];
                populateGUI;
            case cmdDelStark
                if specIndex == [], return; end;
                obj.sobj(specIndex) = [];
                if specIndex > max(length(obj.sobj))
                    specIndex = length(obj.sobj);
                end
                populateGUI;
            case cmdAddAbs
                obj.aobj.loadSpec('sig');
                populateGUI;
            case cmdAddRef
                obj.aobj.loadSpec('ref');
                populateGUI;
            case txtAbsConc
                obj.aobj.setConc(str2num(get(src,'String')));
                populateGUI;
            case txtAbsPath
                obj.aobj.setPath(str2num(get(src,'String')));
                populateGUI;
            case txtAbsInc
                obj.aobj.setInc(str2num(get(src,'String')));
                populateGUI;
            case cmdFit
                set(lblStatus,'String','Fitting Status: working...'); drawnow;
                if get(bgpFit,'SelectedObject') == radAbs
                    obj.fitAbs;
                    populateGUI;
                elseif get(bgpFit,'SelectedObject') == radSim
                    obj.fitSim;
                    populateGUI;
                end
                set(lblStatus,'String','Fitting Status: COMPLETE'); drawnow;
            case tblGauParms
                dat = get(src,'Data');
                obj.gauparms = dat(:,1:3); 
                %make sure max band # is equal to # of transitions
                dat = dat(:,4); dat(dat > size(obj.results,1)) = size(obj.results,1);
                obj.bands = dat;
                obj.tryFit;
                populateGUI;
            case tblStarkParam
                dat = get(src,'Data');
                obj.results = dat(:,1:4);
                populateGUI;
            case cmdAddTransition
                obj.results = [obj.results;[100 50 2 45]];
                axes(axAbs);
                [x,y] = ginput;
                y = y.*x.*obj.opts.weightval;
                w = ones(size(x)).*400;
                if numel(obj.bands) ~= 0
                    b = ones(size(x)).*max(obj.bands)+1;
                    obj.bands = [obj.bands;b];
                    obj.gauparms = [obj.gauparms;[y,x,w]];
                    obj.Achi = [obj.Achi;ones(1,length(obj.sobj)).*1e-21;];
                else
                    b = ones(size(x));
                    obj.bands = b;
                    obj.gauparms = [y,x,w];
                    obj.Achi = ones(1,length(obj.sobj)).*1e-21;
                end
                populateGUI;
            case cmdDelTransition
                obj.results(parmIndex,:) = [];
                obj.gauparms(obj.bands == parmIndex,:) = [];
                obj.bands(obj.bands == parmIndex) = [];
                populateGUI;
            case cmdAchi
                if get(src,'value') == 1
                    set([tblGauParms,...
                        tblStarkParam,cmdAddTransition,cmdDelTransition],'Visible','off');
                    set([tblAchi,cmdResetAchi],'Visible','on');
                else
                    set([tblGauParms,...
                        tblStarkParam,cmdAddTransition,cmdDelTransition],'Visible','on');
                    set([tblAchi,cmdResetAchi],'Visible','off');
                end
            case cmdResetAchi
                obj.Achi = ones(max(obj.bands),length(obj.sobj)).*1e-21;
                populateGUI;
            case tblAchi
                obj.Achi = get(src,'Data');
            case tblOpts
                dat = get(src,'Data');
                obj.opts.weight = dat(1);
                obj.opts.v = dat(2);
                obj.opts.MaxEval = dat(3);
                obj.opts.TolFunAbs = dat(4);
                obj.opts.TolFunSim = dat(5);
                obj.opts.TolX = dat(6);
                populateGUI;
            case cboSolvent
                v = str2num(get(src,'String'));
                obj.aobj.setRefInd(v);
                for j = 1:length(obj.sobj)
                    obj.sobj(j).setRefInd(v);
                end
            case cmdBaseline
                axes(axAbs);
                val = ginput(2); val = val(:,1)';
                set(txtBaseline,'String',num2str(1e7./val));
                obj.aobj.setBaseline(1e7./val);
                populateGUI;
            case txtBaseline
                obj.aobj.setBaseline(str2num(get(src,'String')));
                populateGUI;
            case tblStkBaseline
                bb = get(src,'Data');
                for j = 1:length(bb)
                    obj.sobj(j).setBaseline(bb(j,:));
                end
                populateGUI;
            case cmdLimits
                axes(axStk);
                val = ginput(2); val = val(:,1)';
                set(txtLimits,'String',num2str(1e7./val));
                obj.lims = floor(1e7./val);
                populateGUI;
            case txtLimits
                obj.lims = str2num(get(src,'String'));
                populateGUI;
            case cmdMonteCarlo
                obj.fitMonteCarlo(lblMonteCarlo);
            case txtMCiter
                obj.opts.MCiter = str2num(get(src,'String'));
            case cmdPlotMC
                obj.plotMCresults;
            case cmdError
                set(lblStatus,'String','Fitting Status: working...'); drawnow;
                obj.fitStark;
                populateGUI;
                set(lblStatus,'String','Fitting Status: COMPLETE'); drawnow;
            case cmdExportResults
                if strcmp(projname,'') == 1 %save project before exporting
                    connectGUI(cmdSave,[]);
                end
                obj.export(projname);
                disp('hello')
        end
    end

    function populateGUI() %populates fields in GUI with data from object
    %options
        dat = [obj.opts.weight;obj.opts.v;obj.opts.MaxEval;obj.opts.TolFunAbs; ...
            obj.opts.TolFunSim;obj.opts.TolX];
        set(tblOpts,'Data',dat);
    %Stark spectra
        if numel(obj.sobj) ~= 0
            fname = [obj.sobj.name]'; fname = fname(2:2:length(fname)); %remove path
            pol = num2cell([obj.sobj.pol]');
            inc = num2cell([obj.sobj.inc]');
            F = num2cell([obj.sobj.field]');
            conc = num2cell([obj.sobj.conc]');
            path = num2cell([obj.sobj.path]');
            dat = [fname,pol,inc,F,conc,path];
            set(tblStkSpec,'Data',dat);
        end
    %absorption spectra 
        if numel(obj.aobj.name) ~= 0
            set(lblAbs,'String',obj.aobj.name{1,2});
            if size(obj.aobj.name,1) == 2
                set(lblRef,'String',obj.aobj.name{2,2});
            else
                set(lblRef,'String','');
            end
            set(txtAbsConc,'String',num2str(obj.aobj.conc));
            set(txtAbsPath,'String',num2str(obj.aobj.path));
            set(txtAbsInc,'String',num2str(obj.aobj.inc));
        else
            set(lblAbs,'String','');
        end
    %vibronic parameters
        if size(obj.bands,2) > size(obj.bands,1), obj.bands = obj.bands'; end
        dat = [obj.gauparms,obj.bands];
        set(tblGauParms,'Data',dat);
    %Stark parameters
        dat = obj.results;
        set(tblStarkParam,'Data',dat);
    %solvent
        v = num2str(obj.aobj.n2);
        set(cboSolvent,'String',v);
    %baseline & limits
        set(txtBaseline,'String',num2str(obj.aobj.bsln));
        set(txtLimits,'String',num2str(obj.lims));
        if numel(obj.sobj) ~= 0
            bb = [];
            for j = 1:length(obj.sobj)
                tmpb = obj.sobj(j).bsln;
                if numel(tmpb) == 0
                    tmpb = [0 0];
                end
                if max(tmpb == 0) == 0, tmpb = sort(tmpb); end
                bb(j,:) = tmpb;
            end
            set(tblStkBaseline,'Data',bb);
        end
    %achi
        set(tblAchi,'Data',obj.Achi);
        set(tblAchi,'ColumnEditable',true(1,size(obj.Achi,2)));
        obj.tryFit;
        set(lblWeight,'String',['Absorption Weighting: ',num2str(floor(obj.opts.weightval))]);
    %MonteCarlo options
        set(txtMCiter,'String',num2str(obj.opts.MCiter));
        set(lblMonteCarlo,'String',obj.mc.results);
    end

%% MAIN LOOP ========================================================================
specIndex = []; vibIndex = []; parmIndex = [];
obj = stkFitObj; projname = '';
obj.ax = [axStk,axAbs,axRes,axMain];
populateGUI;

end