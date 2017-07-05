classdef segmentgui < maingui
%Contains methods and properties that are specific to (original) Segment.
 
properties
    AxesTables = [];
    SegmentMeasurements=[];
    LVNO = [];
    RVNO = [];
    FlowNO = [];
    FlowROI = [];
  end
  
  methods
  
    %--------------------------------------
    function g = segmentgui(programversion)
    %--------------------------------------
    %Constructor
    
    %Call superclass constructor
    g = g@maingui(programversion);
    
%     %redo some Gui settings.
%     g.GUISettings.TopGapHeight=0.1350;
%     g.GUISettings.BottomGapHeight= 0.1330;
%     g.GUISettings.LeftGapWidth= 0.1200;
%     g.GUISettings.RightGapWidth=0.2100;
    
    end
    

    %--------------------------
    function initmaintoolbar(g)
    %--------------------------
    %Intialization of main toolbar
    %If new handles are added, add those to list in killhandles.m as well!

    %Create
    g.Handles.maintoolbar = uitoolbar(g.imagefig);
    set(g.fig,'color',g.GUISettings.BackgroundColor); %Added EH:

    %fileopen
    g.addmainicon_helper('segment(''fileopen_Callback'')', ...
      'Open images',...
      g.Icons.opendoc,...
      'fileopenicon');

    %filesave
    g.addmainicon_helper('segment(''filesaveallas_Callback'')',...
      'Save All Image Stacks and Segmentation to Disc',...
      g.Icons.save,...
      'filesaveicon');

    %database
    g.addmainicon_helper('patientdatabase',...
      'Open Segment Patient Database',...
      g.Icons.database,...
      'databaseicon','on');

    %database add patient to
    g.addmainicon_helper('filemenu(''savetopatientdatabase_Callback'')',...
      'Save Image Stacks and Segmentation to Patient Database',...
      g.Icons.databaseadd,...
      'databaseaddicon');

    %pacs
    g.addmainicon_helper('pacs(''init_Callback'')',...
      'Retrieve Images from PACS',...
      g.Icons.connect,...
      'pacsicon','on');

    %pacs add to
    g.addmainicon_helper('filemenu(''savetopacs_Callback'')',...
      'Save Image Stacks and Segmentation to PACS',...
      g.Icons.connectadd,...
      'pacsaddicon');

    %import from CD
    g.addmainicon_helper('utility(''utilitycopyandsortfromcd_Callback'')',...
      'Import Files from CD',...
      g.Icons.importfromcd,...
      'importfromcdicon');

    %View 1 panel
    g.addmainicon_helper('drawfunctions(''drawall'',1)',...
      'View one image panel',...
      g.Icons.panel1,...
      'view1panelicon',...
      'on');

    %View 2 panel
    g.addmainicon_helper('drawfunctions(''drawall'',1,2)',...
      'View two image panels',...
      g.Icons.panel2,...
      'view2panelicon');

    %View 2x1 panel
    g.addmainicon_helper('drawfunctions(''drawall'',2,1)',...
      'View two image panels',...
      g.Icons.panel2x1,...
      'view2x1panelicon');

    %View three*1 panel
    g.addmainicon_helper('drawfunctions(''drawall'',3,1)',...
      'View three image panels',...
      g.Icons.panel3,...
      'view3panelicon');

    %View 1*three panel
    g.addmainicon_helper('drawfunctions(''drawall'',1,3)',...
      'View three image panels',...
      g.Icons.panel1x3,...
      'view1x3panelicon');

    %View four panel
    g.addmainicon_helper('drawfunctions(''drawall'',4)',...
      'View four image panels',...
      g.Icons.panel4,...
      'view4panelicon');

    %View six panel
    g.addmainicon_helper('drawfunctions(''drawall'',6)',...
      'View six image panels',...
      g.Icons.panel6,...
      'view6panelicon');

    %View 9 panel
    g.addmainicon_helper('drawfunctions(''drawall'',9)',...
      'View nine image panels',...
      g.Icons.panel9,...
      'view9panelicon');
    
    %Open view loader
    g.addmainicon_helper('segmentview',...
      'Launch GUI for saving custom views',...
      g.Icons.saveview,...
      'saveviewicon');

    %Show/hide pins
    g.addmainicon_helper('segment(''viewhidepins_Callback'')',...
      'Hide/show pins',...
      g.Icons.hidepins,...
      'hidepinsicon',...
      'on');

    %Show/hide contour from other image stacks.
    g.addmainicon_helper('segment(''viewhideothercontour_Callback'')',...
      'Hide/show LV/RV contour from other image stacks',...
      g.Icons.hideothercontour,...
      'hideothercontouricon',...
      'off');

    %Show/hide interpolation points
    g.addmainicon_helper('segment(''viewhideinterp_Callback'')',...
      'Hide/show interpolation points',...
      g.Icons.hideinterp,...
      'hideinterpicon',...
      'off');

    %Show/hide lv segmentation
    g.addmainicon_helper('segment(''viewhidelv_Callback'')',...
      'Hide/show left ventricle segmentation',...
      g.Icons.hidelv,...
      'hidelvicon',...
      'off');

    %Show/hide rv segmentation
    g.addmainicon_helper('segment(''viewhiderv_Callback'')',...
      'Hide/show right ventricle segmentation',...
      g.Icons.hiderv,...
      'hidervicon',...
      'off');

    %Show/hide scar
    g.addmainicon_helper('segment(''viewhidescar_Callback'')',...
      'Hide/show scar',...
      g.Icons.hidescar,...
      'hidescaricon',...
      'off');

    %Show/hide mar
    g.addmainicon_helper('segment(''viewhidemar_Callback'')',...
      'Hide/show MaR',...
      g.Icons.hidemar,...
      'hidemaricon',...
      'off');

    %Show/hide roi
    g.addmainicon_helper('segment(''viewhideroi_Callback'')',...
      'Hide/show roi''s',...
      g.Icons.hideroi,...
      'hideroiicon',...
      'off');

    %Show/hide measures
    g.addmainicon_helper('segment(''viewhidemeasures_Callback'')',...
      'Hide/show measurements',...
      g.Icons.hidemeasure,...
      'hidemeasuresicon',...
      'off');
    
    %Show/hide annotation points
    g.addmainicon_helper('segment(''viewhidepoints_Callback'')',...
      'Hide/show annotation points',...
      g.Icons.hidepoint,...
      'hidepointsicon',...
      'off');

    %Show/hide plus/center point
    g.addmainicon_helper('segment(''viewhideplus_Callback'')',...
      'Hide/show image center point',...
      g.Icons.hideplus,...
      'hideplusicon',...
      'off');

    %Show/hide intersections
    g.addmainicon_helper('segment(''viewhideintersections_Callback'')',...
      'Hide/show image plane intersections',...
      g.Icons.hideintersections,...
      'hideintersectionsicon',...
      'off');

    %Show/hide text
    g.addmainicon_helper('segment(''viewhidetext_Callback'')',...
      'Hide/show text',...
      g.Icons.hidetext,...
      'hidetexticon',...
      'off');

    if isequal(g.GUISettings.PapilarColor,1)
      g.addmainicon_helper('segment(''viewhidepap_Callback'')',...
        'Hide/show papillary overlay',...
        g.Icons.hidepapwhite,...
        'hidepapicon','off');
    else
      g.addmainicon_helper('segment(''viewhidepap_Callback'')',...
        'Hide/show papillary overlay',...
        g.Icons.hidepapblack,...
        'hidepapicon','off');
    end
        %show/hide manual interaction overlay
    g.addmainicon_helper('segment(''viewhidemanualinteraction_Callback'')',...
      'Hide/show manual interaction for Scar/MaR',...
      g.Icons.hideoverlaynew,...
      'hideoverlayicon','off');
        
    %Show colorbar
    g.addmainicon_helper('segment(''viewhidecolorbar_Callback'')',...
      'Show colorbar in current image stack',...
      g.Icons.colorbar,...
      'colorbaricon',...
      'off');
    
    %View pixely/interpolated
    g.addmainicon_helper('segment(''viewinterp_Callback'')',...
      'View image pixels',...
      g.Icons.viewpixely,...
      'viewpixelyicon',...
      'off');

    %Help icon
    props.ClickedCallback = 'segmenthelp(''about_Callback'')';
    props.ToolTip = 'About';
    props.CData = g.Icons.help;
    props.Tag = 'about';
    props.separator = 'on';
    g.Handles.abouticon = uitoggletool(g.Handles.maintoolbar,props);
    end
    
        
    %--------------------------
    function measurementinit(g)
    %--------------------------
    %init objects for measurement panel in AxesTables
    
%     set(g.Handles.lvuipanel,'BackgroundColor','black');
%     set(g.Handles.flowuipanel,'BackgroundColor','black');
%     set(g.Handles.reportpanel,'BackgroundColor','black');
%     set(g.Handles.barpanel,'BackgroundColor','black');
    
    bl = '---'; %no values on Init
    bl2 = {bl,bl};
    
    g.AxesTables.volume = axestable(g.Handles.volumeresultaxes);
    g.AxesTables.flow = axestable(g.Handles.flowresultaxes);
    g.AxesTables.measurement = axestable(g.Handles.measurementresultaxes);
    
    m=g.getmeasurementstruct();
    
    %LV and RV report table  
    g.AxesTables.volume.backgroundcolor = [0.94 0.94 0.94]; %[0 0 0];%
    g.AxesTables.volume.fontcolor = [0 0 0]; % [1 1 1];%
    g.AxesTables.volume.fontsize = 10;
    g.AxesTables.volume.ystep = 15;
%     g.AxesTables.volume.FigWith = 300;
    g.AxesTables.volume.addTable('LV',2,1,[0.6 0.4]);
    g.AxesTables.volume.addKey('LVM','ED/ES LVM','g',bl);
    g.AxesTables.volume.addKey('LVEDV','EDV','ml',bl);
    g.AxesTables.volume.addKey('LVESV','ESV','ml',bl);
    g.AxesTables.volume.addKey('LVSV','SV','ml',bl);
    g.AxesTables.volume.addKey('LVEF','EF','%',bl);
    g.AxesTables.volume.addKey('LVHR','HR','bpm',bl);
    
%     g.AxesTables.volume.addSpace();
    g.AxesTables.volume.addTable('RV',2,1,[0.6 0.4]);
    g.AxesTables.volume.addKey('RVM','ED/ES RVM','g',bl);
    g.AxesTables.volume.addKey('RVEDV','EDV','ml',bl);
    g.AxesTables.volume.addKey('RVESV','ESV','ml',bl);
    g.AxesTables.volume.addKey('RVSV','SV','ml',bl);
    g.AxesTables.volume.addKey('RVEF','EF','%',bl);
    g.AxesTables.volume.addKey('RVHR','HR','bpm',bl);
    
    %Flow report table
    g.AxesTables.flow.backgroundcolor = [0.94 0.94 0.94]; %[0 0 0];%
    g.AxesTables.flow.fontcolor = [0 0 0]; %[1 1 1];%
    g.AxesTables.flow.ystep = 20;%16;
    g.AxesTables.flow.fontsize = 10;
    g.AxesTables.flow.addTable('Flow',2,1,[0.65 0.35]);
    g.AxesTables.flow.addKey('ROI','ROI','','');
    g.AxesTables.flow.addKey('Netvol','Net vol','ml',bl);
    g.AxesTables.flow.addKey('Forward','Forward','ml',bl);
    g.AxesTables.flow.addKey('Backward','Backward','ml',bl);
    g.AxesTables.flow.addKey('Regfrac','Regurg. frac.','%',bl);
    g.AxesTables.flow.addKey('FlowHR','HR','bpm',bl);
   
    %measurement report table
    g.AxesTables.measurement.backgroundcolor = [0.94 0.94 0.94]; %[0 0 0];%
    g.AxesTables.measurement.fontcolor = [0 0 0]; %[1 1 1];%
    g.AxesTables.measurement.ystep = 20;%16;
    g.AxesTables.measurement.fontsize = 10;
    g.AxesTables.measurement.addTable('Measurement',2,1,[0.6 0.4]);
    g.AxesTables.measurement.addKey('m1',bl,'',bl);
    g.AxesTables.measurement.addKey('m2',bl,'',bl);
            
    end
    
    %----------------------
    function inittoolbar(g)
    %----------------------
    g.GUISettings.VolumeAxesColor = [0 0 0];%[1 1 1];%
    g.GUISettings.VolumeColorGraph = [1 1 1];%[0 0 0];%
    g.GUISettings.TopGapHeight = 0.135;
    g.GUISettings.BottomGapHeight = 0.133;
    g.GUISettings.RightGapWidth = 0.21;
    g.GUISettings.ShowColorbar =  false;
    
    %Initiate AxesTables objects
    g.measurementinit;
    
    %Clear all measurements
    g.updateaxestables('flowclearall');
    g.updateaxestables('volumeclearall');
    
    %initateplugins
    segment('initmenu');
    
%     if isequal(36,getmodule(35,'KI',[],true));
%       uimenu(g.Handles.analysismenu,'Label','RAMP', ...
%         'Callback','plugin_RAMP(''init'')');
%     end
    end
    
    
    %---------------------------------
    function m=getmeasurementstruct(g)
    %---------------------------------
    %Set up the structure of the result tables in main interface
    %For now mostly copied from RVQ
    
    %category   type  code  shortname longname                        showinpanel showinreport
    %LV         mm    ASW   Asw       Anterior Septal Wall Thickness  1
    
    if isempty([]) %Replace property CVQMeasurements
      try
        all=textread('SegmentMeasurements.txt','%s','delimiter',';');
      catch
        myfailed('Could not read SegmentMeasurements.txt file.');
        return;
      end;

      write_pos=1;
      read_pos=1;
      m=cell(1,1);

      while read_pos<length(all)
        comment = false;
        if isequal(all{read_pos},'#')
          comment = true;
        end;
        if not(comment)
          switch char(all(read_pos+1))
            case 'ta'%parse a table; 10 args
              nbr_args=10;
              values_loop=6:10;
            case 'sp'%parse a special; 9 args
              nbr_args=9;
              values_loop=6:7;
            case {'mm','roi'}%parse a other; 7 args
              nbr_args=7;
              values_loop=6:7;
            case 'em'
              nbr_args=2;
              values_loop=[];
            otherwise
              mywarning('Parser error of SegmentMeasurements.txt panel will probaly not look good',g.GUI.Segment);
              return;
          end
          m{write_pos}=all(read_pos:read_pos+nbr_args-1);
          for i=values_loop;
            m{write_pos}{i}=eval(m{write_pos}{i});
          end
          read_pos=read_pos+nbr_args;
          write_pos=write_pos+1;
        else
          %is comment
          read_pos=read_pos+2;
        end;
      end
      g.SegmentMeasurements=m;
    else
      m=g.SegmentMeasurements;
    end

    end
    
    
    %------------------------------------------
    function updateaxestables(g, arg, varargin)
    %------------------------------------------
    %Method to update AxesTables, in GUIs where present
    global NO SET
    
    g.updatetimebaraxes;

    %In segment we want to show current stack values. 
    g.LVNO=NO;
    g.RVNO=NO;
    
    if ~isempty(SET) && ~isempty(SET(NO).Flow)
      g.FlowNO=NO;
    else
      g.FlowNO=[];
    end
    
    switch arg
      case 'volume'
        g.updatevolumeaxes;
        g.volumereportupdate(varargin{:});
      case 'volumeclearall'
        g.updatevolumeaxes;
        g.volumereportupdate(varargin{:});
      case {'area','flow'}
        g.updateflowaxes;
        g.flowreportupdate;
      case {'areaclearall','flowclearall'}
        g.updateflowaxes;
        g.flowreportupdate;
      case {'measure'}
        g.measurementreportupdate;
      case {'measureclearall'}
        g.measurementreportupdate;
    end
    
    end
    
    
    %-----------------------------
    function volumereportupdate(g)
    %-----------------------------
    %Update LV and RV report in axestables
    global SET DATA
    
      m=g.getmeasurementstruct();
      bla = {'---'};
      
      % -- LV --
      
      s = [];
      s.LVM = bla;
      s.LVEDV = bla;
      s.LVESV = bla;
      s.LVSV = bla;
      s.LVEF = bla;
      s.LVHR = bla;
      updatestruct = s;
      %       title{1} = 'LV';
        
      if ~isempty(DATA.LVNO) % isfield(DATA,'LVNO') && 
        no = DATA.LVNO;
        haslv = ~isempty(SET(no).EpiX) || ~isempty(SET(no).EndoX);
        
        if haslv          
          %Determine last character of key
          prepost = 'LV';
          if isnan(SET(no).LVM(SET(no).EDT))
            lvmedt = bla{1};
          else
            lvmedt = num2str(round(SET(no).LVM(SET(no).EDT)*1.05));
          end
          if isnan(SET(no).LVM(SET(no).EST))
            lvmest = bla{1};
          else
            lvmest = num2str(round(SET(no).LVM(SET(no).EST)*1.05));
          end
          updatestruct.([prepost 'M']) = sprintf('%s / %s',lvmedt,lvmest);
          updatestruct.([prepost 'EDV']) = round(SET(no).EDV);
          updatestruct.([prepost 'ESV']) = round(SET(no).ESV);
          updatestruct.([prepost 'SV']) = round(SET(no).SV);
          updatestruct.([prepost 'EF']) = round(100*SET(no).EF);
%           title{1} = sprintf('LV  [Image Stack %d]',no);
        end
        updatestruct.LVHR = round(SET(no).HeartRate);
        %set(DATA.Handles.lvstackpushbutton,'String',sprintf('Stack #%d',no));
      else        
        %set(DATA.Handles.lvstackpushbutton,'String','Set stack');
      end
          
      fnames = fieldnames(updatestruct);
      for floop = 1:numel(fnames)
        fname = fnames{floop};
        g.AxesTables.volume.updateKey(fname,updatestruct.(fname),true);
      end
    
      % -- RV --
      
      s = [];
      s.RVM = bla;
      s.RVEDV = bla;
      s.RVESV = bla;
      s.RVSV = bla;
      s.RVEF = bla;
      s.RVHR = bla;
      updatestruct = s;
      %       title{2} = 'RV';
        
      if ~isempty(DATA.RVNO) % isfield(DATA,'RVNO') && 
        no = DATA.RVNO;
        hasrv = ~isempty(SET(no).RVEpiX) || ~isempty(SET(no).RVEndoX);
        
        if hasrv          
          %Determine last character of key
          prepost = 'RV';   
          if isnan(SET(no).RVM(SET(no).EDT))
            rvmedt = bla{1};
          else
            rvmedt = num2str(round(SET(no).RVM(SET(no).EDT)*1.05));
          end
          if isnan(SET(no).RVM(SET(no).EST))
            rvmest = bla{1};
          else
            rvmest = num2str(round(SET(no).RVM(SET(no).EST)*1.05));
          end
          updatestruct.([prepost 'M']) = sprintf('%s / %s',rvmedt,rvmest);
          updatestruct.([prepost 'EDV']) = round(SET(no).RVEDV);
          updatestruct.([prepost 'ESV']) = round(SET(no).RVESV);
          updatestruct.([prepost 'SV']) = round(SET(no).RVSV);
          updatestruct.([prepost 'EF']) = round(100*SET(no).RVEF);
%           title{2} = sprintf('RV  [Image Stack %d]',no);
        end
        updatestruct.RVHR = round(SET(no).HeartRate);
        %set(DATA.Handles.rvstackpushbutton,'String',sprintf('Stack #%d',no));
      else        
        %set(DATA.Handles.rvstackpushbutton,'String','Set stack');
      end
      
      fnames = fieldnames(updatestruct);
      for floop = 1:numel(fnames)
        fname = fnames{floop};
        g.AxesTables.volume.updateKey(fname,updatestruct.(fname),true);
      end
%       g.AxesTables.volume.updateTitle(title);      
      g.AxesTables.volume.draw();
    end
    
    
    %---------------------------
    function flowreportupdate(g)
    %---------------------------
    %Update flow report in axestables
    
    global SET DATA
    
    updatestruct = struct(...
      'ROI', '',...
      'Netvol', '---',...
      'Forward', '---',...
      'Backward', '---',...
      'Regfrac', '---',...
      'FlowHR', '---');
    name.ROI = 'ROI';
    %set(DATA.Handles.flowstackpushbutton,'String','Set stack');
%     title{1} = 'Flow';
    
    if ~isempty(DATA.FlowNO) && ~isempty(DATA.FlowROI)
      no = DATA.FlowNO;
      roinbr = DATA.FlowROI;
      if ~isempty(SET) && ~isempty(no) && SET(no).RoiN > 0 && ~isempty(SET(no).Flow.Result) && not(isempty(roinbr)) && length(SET(no).Flow.Result)>=roinbr && isfield(SET(no).Flow.Result(roinbr),'nettotvol')
        
        updatestruct.ROI = '';
        updatestruct.Netvol = SET(no).Flow.Result(roinbr).nettotvol;
        updatestruct.Forward = SET(no).Flow.Result(roinbr).netforwardvol;
        updatestruct.Backward = SET(no).Flow.Result(roinbr).netbackwardvol;
        updatestruct.Regfrac = SET(no).Flow.Result(roinbr).regfrac;
        updatestruct.FlowHR = SET(no).HeartRate;
        nop = SET(no).Flow.PhaseNo;
        if ~isempty(SET(nop).Flow.PhaseCorr)
          if ~isfield(SET(nop).Flow,'PhaseCorrTimeResolved')
            mywarning('Incompatible eddy current correction. Correction reset.',DATA.GUI.Segment);
            SET(nop).Flow.PhaseCorr = [];
            name.ROI = SET(no).Roi(roinbr).Name;
          else
            name.ROI = sprintf('%s [Eddy]',SET(no).Roi(roinbr).Name);
          end
        else
          name.ROI = SET(no).Roi(roinbr).Name;
        end
        %set(DATA.Handles.flowstackpushbutton,'String',sprintf('Stack #%d',no));
      end
    end      
    DATA.updateflowaxes;
    
    g.AxesTables.flow.updateName('ROI',name.ROI,true);
    fnames = fieldnames(updatestruct);
    for floop = 1:numel(fnames)
      fname = fnames{floop};
      g.AxesTables.flow.updateKey(fname,updatestruct.(fname),true);
    end
%     g.AxesTables.flow.updateTitle(title);
    g.AxesTables.flow.draw();
            
    end
    
    
    %---------------------------------
    function measurementreportupdate(g)
    %---------------------------------
    %Update scar report in axestables
    
    global SET NO
                
      no = NO;
      hasscar = ~isempty(SET(no).Scar) && ~isempty(SET(no).Scar.Result);
      hasmar = ~isempty(SET(no).MaR) && ~isempty(SET(no).MaR.Result);
      hasflow = ~isempty(SET) && SET(no).RoiN > 0 && isfield(SET(no).Flow,'nettotvol');
      hasmmode = strcmp(g.ViewPanelsType{g.CurrentPanel},'mmodespatial');
      updatestruct = [];
      if hasmar && hasscar
        updatestruct.m1 = round(SET(no).Scar.Percentage);%([SET(no).EDT, SET(no).EST]));
        updatestruct.m2 = round(SET(no).MaR.Percentage([SET(no).EDT,SET(no).EST]));
        name.m1 = 'Scar';
        unit.m1 = '%';
        name.m2 = 'MaR ED ES';
        unit.m2 = '%';
        title{1} = sprintf('Scar & MaR');
           
      elseif hasscar
        %updatestruct.m1 = round(SET(no).LVM(SET(no).EDT)*1.05);
        updatestruct.m1 = round(SET(no).Scar.Percentage);%([SET(no).EDT, SET(no).EST]));
        %name.m1 = 'LVM';
        name.m1 = 'Scar';
        %unit.m1 = 'g';
        unit.m1 = '%';
        title{1} = sprintf('Scar');
      elseif hasmar
        %updatestruct.m1 = round(SET(no).LVM(SET(no).EDT)*1.05);
        updatestruct.m1 = round(SET(no).MaR.Percentage([SET(no).EDT,SET(no).EST]));
        %name.m1 = 'LVM';
        name.m1 = 'MaR ED ES';
        %unit.m1 = 'g';
        unit.m1 = '%';
        title{1} = sprintf('MaR');  

      elseif hasflow
        rloop = 1;
        updatestruct.m1 = SET(no).Flow.nettotvol(rloop);
        updatestruct.m2 = SET(no).Flow.regfrac(rloop);
        name.m1 = 'Net vol';
        name.m2 = 'Regurg. frac.';
        unit.m1 = 'ml';
        unit.m2 = '%';
        title{1} = sprintf('Flow');
      elseif hasmmode
        [dist,timedist] = calcfunctions('calcmmodedists',no);
        updatestruct.m1 = dist;
        updatestruct.m2 = timedist;
        name.m1 = 'Distance';
        name.m2 = 'Time';
        unit.m1 = 'mm';
        unit.m2 = 'ms';
        title{1} = sprintf('Mmode: Distance between lines');
      else
        bla = {'---'};
        updatestruct.m1 = bla;
        updatestruct.m2 = bla;
        name.m1 = bla;
        name.m2 = bla;
        unit.m1 = '';
        unit.m2 = '';
        title{1} = '';
      end
          
      if hasscar || hasmar || hasflow || hasmmode
        set(g.Handles.measurementuipanel,'Visible','on');
        g.AxesTables.measurement.show;
      end
      
      fnames = fieldnames(updatestruct);
      for floop = 1:numel(fnames)
        fname = fnames{floop};
        g.AxesTables.measurement.updateKey(fname,updatestruct.(fname),true);
        g.AxesTables.measurement.updateName(fname,name.(fname),true);
        g.AxesTables.measurement.updateUnit(fname,unit.(fname),true);
      end
      g.AxesTables.measurement.updateTitle(title);
      g.AxesTables.measurement.draw();
      
      if not(hasscar) && not(hasmar) && not(hasflow) && not(hasmmode)
        g.AxesTables.measurement.hide;
        set(g.Handles.measurementuipanel,'Visible','off');
%       else        
%         xval = get(g.Handles.measurementresultaxes,'xlim');
%         yval = get(g.Handles.measurementresultaxes,'ylim');
%         hold(g.Handles.measurementresultaxes,'on');
%         plot(g.Handles.measurementresultaxes,[xval(1) xval(2) xval(2) xval(1) xval(1)],[yval(1) yval(1) yval(2) yval(2) yval(1)],'k-');
%         hold(g.Handles.measurementresultaxes,'off');
      end
      
    end
    
      
    %------------------------
    function failedaborted(g)
    %------------------------
    %Message when a function is aborted.
    myfailed('Aborted by user.',g.GUI.Segment);
    end
      
    %------------------------------------
    function updateprefhandlesadvanced(g)
    %------------------------------------
    %Make update of handles in advanced preferences GUI
    updateprefhandlesadvanced@maingui(g);
    g.normalizephaseupdate;
    end
    
    %-----------------------
    function checkversion(g)
    %-----------------------
    %Check if new version is available.
    if g.Pref.CheckVersion
      checkversion();
    end
    end
    
    %---------------------------
    function updateicons(g,mode)
    %---------------------------
    %Update icons when new mode is selected
    %Overloads main GUI method.
    set([g.Handles.t2starpushbutton g.Handles.perfusionpushbutton ], 'visible', 'off', 'enable','off');
    if nargin<2
      mode = g.CurrentTheme;
    end;
    updateicons@maingui(g,mode);
    end
    
  end
  
  methods(Access = 'protected')
    
    %---------------------------------------
    function switchtoimagestack_part2(g, no)
    %---------------------------------------
    %switchtoimagestack continued. Overloads method in maingui.
    if ~isempty(g.BlockingFigs)
      deleteguis=yesno('Need to close all reviewing windows before changing current image stack. Do you want to do that now?',[],g.GUI.Segment);
      if deleteguis
        try
          delete(g.BlockingFigs);
        catch %#ok<CTCH>
        end;
      end
    end;
    switchtoimagestack_part2@maingui(g,no);
    end
    
  end
  
  methods(Static)
    
    %---------------------
    function licensemsg(m)
    %---------------------
    %Display license message
    
    %h = msgbox(m,'License');
    %set(h,'windowstyle','modal');
    %uiwait(h);
    
    %mymsgbox(m,'License')
    versionhello('init',m,'License');
    end
    
    %----------------------------
    function tip = singleframetip
    %----------------------------
    %Return tip for single frame mode
    tip = 'single frame mode should be used';
    end
    
    %------------------------------------
    function out = viabilityallowweighted
    %------------------------------------
    %Return whether to allow weighted viability mode
    out = true;
    end
    
    %-------------------------
    function copyrvendoforward
    %-------------------------
    %Overloads empty method in maingui
    global SET NO
    tf = SET(NO).CurrentTimeFrame;
    tf2 = mod(tf,SET(NO).TSize) + 1;
    slc = SET(NO).CurrentSlice;
    SET(NO).RVEndoX(:,tf2,slc) = SET(NO).RVEndoX(:,tf,slc);
    SET(NO).RVEndoY(:,tf2,slc) = SET(NO).RVEndoY(:,tf,slc);
    segment('setcurrenttimeframe',tf2);
    segment('updatevolume');
    end
  end
  
end