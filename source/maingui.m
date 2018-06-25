classdef maingui < handle %Handle class
  
  properties
    LogFile = [];
    UserFile = [];
    ProgramVersion = [];
    ProgramName = 'Segment';
    ProgramFolderName = 'Segment';
    fig = [];
    BlockingFigs = [];
    imagefig = [];
    DataLoaded = 0;
    Locked = 0;
    Silent = 0;
    Handles = [];

    Interactionlock = false;
    Undo = []; %length(SET)*g.Pref.UndoHistory
    UndoN = [];
    LastSaved = now;
    Volrend = false;
    ViewPanels = [];
    ViewPanelsType = {};
    ViewPanelsMatrix = {};
    ViewIM = {};
    ViewMatrix = [];
    LastView = [];
    
    CurrentPanel = 1;
    CurrentTool = 'select';
    CurrentTheme = 'lv';
    Tools = [];

    SegmentFolder = pwd;
    
    imagedescriptionfile = 'imagedescription.txt';
    manualfile = 'manual.pdf';
    IFUfile = '';
    
    GUI = [];
    GUISettings = [];
    GUIPositions = [];
    
    AllowInt16 = false;
    Pref = [];
    PrefHandles = [];
    PrefHandlesAdvanced = [];
    PrefHandlesPacs = [];
    Synchonize = 0;

    Preview = [];
    
    %Video tools
    StartTime = [];
    Run = false;
    Record = false;
    LastPointer = 'arrow';
    LastPointerShapeCData = [];
    VisibleThumbnails = [];
    
    Icons = [];
    
    ImageTypes = [];
    ImageViewPlanes = [];
    ImagingTechniques = [];
    ImagingTechniquesFullNames = [];
    
    ThisFrameOnly = 1;
    Synchronize=0;
    StartFrame = 1;
        
    SegIntersection = [];
    IntersectionIM = [];
    ShowIntersectionMask = false;
    
    BalloonLevel = 0;
    BALLOON = [];
		EndoBalloonForce = [];
		EpiBalloonForce = [];
    LevelSet = [];
    DATASETPREVIEW = [];
    EndoEDGE0 = [];
    EndoEDGE1 = [];
    EndoEDGE2 = [];
    EndoEDGE3 = [];
    EpiEDGE0 = [];
    EpiEDGE1 = [];
    EpiEDGE2 = [];
    EpiEDGE3 = [];
    EndoEdgeDetected = 0;
    EpiEdgeDetected = 0;
    MovieFrame = [];
    NumPoints = 80;
    
    Overlay = [];
    
    BpInt = 0;
    MInt = 0;
    TInt = 0;
    
    Pin = [];
    
    MeasureX = [];
    MeasureY = [];
    MeasureT = [];
    MeasureZ = [];
    MeasureN = [];
    MeasureName = [];
    MeasureOffsetX = [];
    MeasureOffsetY = [];
    
    CursorX = [];
    CursorY = [];
    CursorN = [];
    CursorXOfs = [];
    CursorYOfs = [];
    
    NeedToSave = 0;
    
    Testing = false; %This property is used by maketest only.
    RecordMacro = false; %Used by macro_helper;
    Macro = []; %Used by macro_helper
    MacroN = []; %Used by macro_helper
    Buffer = []; %Used for macros and testing
    
    CineTimer = [];
    SegmentServerMonitorTimer = [];
    
    contourp = [];
    contourbimage = []; %used for dragging contours in lv.m
    
    MainMotionFcn = '';
    
    DynamicPACS = [];
    
    LicenseNumber = '';
    PayPerCaseS = [];
    
        
  end
  
  methods
    
    %--------------------------------------
    function g = maingui(programversion)
    %--------------------------------------
    %Constructor
    %Initialization of Segment GUI
    
    %System settings
    
    %Removed this when we switched to R2014b since 
    %opengl software %This avoids the bullseye plot to get "partly upside down" on some computers.

    %Removed this when we switched to R2011a since this the bug with race
    %condition in logical operands is elimenated since R2008b according to
    %Mathworks. EH:
    %maxNumCompThreads(1); %This (hopefully!) eliminates a loader bug problem, see ticket #598.

    pathname = getpreferencespath;
    
    g.ProgramVersion = programversion;
    g.Buffer = clearbuffer;  %Clear and initializes buffers used for macro purposes    

    %Remove old logfiles
    f = dir([pathname filesep '*.log']);
    for loop=1:length(f)
      if (now-f(loop).datenum)>7
        delete([pathname filesep f(loop).name]);
      end;
    end;

    %Remove old inbox files
    f = dir([pathname filesep '*.inbox']);
    for loop=1:length(f)
      if (now-f(loop).datenum)>7
        delete([pathname filesep f(loop).name]);
      end;
    end;

    %Remove old outbox files
    f = dir([pathname filesep '*.outbox']);
    for loop=1:length(f)
      if (now-f(loop).datenum)>7
        delete([pathname filesep f(loop).name]);
      end;
    end;

    g.initLogFile;
    g.dispwelcometext;    

    % Register where segment.m is located if we are running from source
    % from another directory (JT 090402)
    if isdeployed()
      g.SegmentFolder = pwd;
    else
      [segmentfolder_, name_, ext_] = fileparts(which('segment')); %#ok<NASGU>
      g.SegmentFolder = segmentfolder_;
      clear segmentfolder_ name_ ext_ 
    end

    checkpath(g.SegmentFolder); %Ensure that segment directory is on the path
    %This is important if the user switches directory, but even more important
    %is that the package technology will not work. Matlab will crash badly when
    %exiting on packaged functions when one change directory, and some other
    %unknown places.

    %Here all that uses mygui should be added.
    g.GUI.Segment = [];
    g.GUI.AVPlane = [];
    g.GUI.BrukerReader = [];
    g.GUI.Bullseye = [];
    g.GUI.CESSFPMaR = [];
    g.GUI.Communicate = [];
    g.GUI.CropLVStacks = [];
    g.GUI.CTPerfusion = [];
    g.GUI.DicomSorter = [];
    g.GUI.ECVRegistration = [];
    g.GUI.Flow = [];
    g.GUI.FlowUnwrap = [];
    g.GUI.Fusion = [];
    g.GUI.GreyZoneHist = [];
    g.GUI.GreyZoneHistAlt = [];
    g.GUI.GObject = [];
    g.GUI.GObjectSpeedIM = [];
    g.GUI.Hotkeys = [];
    g.GUI.ImageStacksCT = [];
    g.GUI.InterchangeDelineation =[];
    g.GUI.LCS2D = [];
    g.GUI.LCS2DFlowCurves = [];
    g.GUI.LevelSet = [];
    g.GUI.LevelSetSpeedIM = [];
    g.GUI.LongAxis = []; %GUI for longaxis plots
    g.GUI.LVLongaxismotion = [];
    g.GUI.LVSegmentation = [];
    g.GUI.MOSlider = [];
    g.GUI.OpenFile = [];
    g.GUI.PacsCon = [];
    g.GUI.PacsTransfer = [];
    g.GUI.PapillarySlider = [];
    g.GUI.PatientDataBase = [];
    g.GUI.PatientDataBaseSelectStudies = [];
    g.GUI.PatientDataBaseSendFiles = [];
    g.GUI.PatientDataBaseSelectFiletype = [];
    g.GUI.PatientInfo = [];
    g.GUI.Perfusion = [];
    g.GUI.PerfusionGraphs = [];
    g.GUI.PerfusionMR = [];
    g.GUI.PerfusionScoring = [];
    g.GUI.Pressure = [];
    g.GUI.PulseWaveVelocity = [];
    g.GUI.Reformat = [];
    g.GUI.ReportSheetGenerator = []; %GUI for generating reports
    g.GUI.Report3DModel = [];
    g.GUI.ROI = [];
    g.GUI.RVSegmentation = [];
    g.GUI.SegmentView = [];
    g.GUI.SeriesSelector=[];
    g.GUI.SetImageDescription =[];
    g.GUI.SetImageInfo = [];
    g.GUI.SliceReport = [];
    g.GUI.SpectPlot3d = [];
    g.GUI.SpectPlot2d = [];
    g.GUI.SpectPerfusion = [];
    g.GUI.SpectRegistration = [];
    g.GUI.Spectconfirmcenter = [];
    g.GUI.Spectdefectpreferences = [];
    g.GUI.Strain = [];
    g.GUI.StrainTagging = [];
%     g.GUI.StrainShortAxis = [];
%     g.GUI.StrainPolarPlot = [];
    g.GUI.T2Star = [];
    g.GUI.T2wMaR = [];
    g.GUI.TransitTimeTool = [];
    g.GUI.VolumePlot = [];
    
    g.GUISettings = [];
    g.GUISettings.TopGapHeight=0.01;
    g.GUISettings.BottomGapHeight = 0.013;
    g.GUISettings.LeftGapWidth =0.12;
    g.GUISettings.RightGapWidth=0.21;
    g.GUISettings.ReportPanelPixelMax=220;
    g.GUISettings.ThumbnailPanelWidth=100;
    g.GUISettings.BackgroundColor=[0.94 0.94 0.94];%Will later be overwritten by preferences.
    g.GUISettings.ButtonColor = [0.94 0.94 0.94];%[0.92 0.91 0.84];
    g.GUISettings.ButtonSelectedColor = [0.5 0.5 1]; %sometimes [0.6 0.6 0.6]
    g.GUISettings.AxesColorDefault=[1 0.6 0.2];
    g.GUISettings.AxesColor = [1,1,1];%g.GUISettings.AxesColorDefault;
    g.GUISettings.TimebarAxesColor = [0,0,0];
    g.GUISettings.BarColorDefault=[55 119 106]/256;
    g.GUISettings.BarColor=[55 119 106]/256;

    g.GUISettings.ThumbLineColor = [1 0.6 0.2];
    g.GUISettings.ThumbFlowLineColor = [0.9 0.9 0.9];
    g.GUISettings.MontageBorder=true;

    g.GUISettings.IntersectColor=[1 0.4 0.2];
    g.GUISettings.BoxAxesLineSpec='-';
    g.GUISettings.SliceLineSpec='y-';
    g.GUISettings.SliceLineSpecMultiSlice='y-';
    g.GUISettings.SliceLineSpecOneSlice='y-';
    g.GUISettings.MeasureLineSpec='w+-';
    g.GUISettings.MeasureLineMarkerSize=6;
    g.GUISettings.BoxAxesColor=[0.1 0.1 0.1];%[0 0 0]; %[0.7 0.2 0.2];
    g.GUISettings.VolumeAxesColor = [0 0 0];
    g.GUISettings.ViewPanelsTypeDefault='one';
    g.GUISettings.AskToExitProgram = true;
    g.GUISettings.PointsEnabled = true;
    g.GUISettings.ColorMapSize = 256;
    g.GUISettings.ThumbnailSize = 256;

    %Added by eriks
    g.AllowInt16=false;
    g.GUISettings.DrawPointer='cross';
    g.GUISettings.MeasurePointer='cross';
    g.GUISettings.PapilarColor=0;%this is casted to uint8 and multiplied by 255
    g.GUISettings.VolumeColorGraph=[1 1 1];
    g.GUISettings.DefaultROISpec='b-';
    g.GUISettings.DefaultROIDrawColor='b';
    
    %GUI positions default values
    g.GUIPositions(1).FileName='monitor';
    g.GUIPositions(1).Position=get(0,'MonitorPosition');%normalized position
    
    g.GUIPositions(2).FileName='segment.fig';
    g.GUIPositions(2).Position=[0.05 0.05 0.9 0.85];%normalized position
    
    

    %load saved GUIpositions
    try
      g.loadguipositions;
    catch
    end

    end
    
    
    %-------------------------------------------------------
    function hidetimebar_Callback(varargin)
    %-------------------------------------------------------
    end
    
    %-------------------------------------------------------
    function hideresult_Callback(varargin)
    %-------------------------------------------------------
    
    end
    
    %-------------------------------------------------------
    function hideicons_Callback(varargin)
    %-------------------------------------------------------
    end
    
    %--------------------------------------------------------------
    function inithideplaceholder(varargin)
    %-----------------------------------------------------------
    g=varargin{1};
    
    hideiconcell={};
    hideiconcell{1,end+1}=myicon('synchronize',g.Handles.hideiconholder,g.Icons.config.syncronize,'Synchronize all image stacks',@() segment('synchronize'),2);
    hideiconcell{1,end+1}=myicon('singleframemode',g.Handles.hideiconholder,g.Icons.config.singleframemode,'Single frame mode',@() segment('thisframeonly_Callback',1,0) ,1,1);
    hideiconcell{1,end+1}=myicon('allframesmode',g.Handles.hideiconholder,g.Icons.config.allframesmode,'All frames mode',@()  segment('thisframeonly_Callback',0,0),1,1);
    hideiconcell{1,end+1}=myicon('hideicons',g.Handles.hideiconholder,g.Icons.config.hideiconpanel,'Hide icons',@() segment('hidepanels'),2);
    hideiconcell{1,end+1}=myicon('hideresult',g.Handles.hideiconholder,g.Icons.config.hideresultpanel,'Hide result panel',@() segment('hidepanels'),2);
    g.Handles.hideiconholder.add(hideiconcell);
    pos=plotboxpos(g.Handles.hideiconholder.axeshandle);
    currentpos=get(g.Handles.hideiconholder.axeshandle,'position');
    currentpos(1)=currentpos(1)+1-(pos(1)+pos(3));
    set(g.Handles.hideiconholder.axeshandle,'position',currentpos);%-[pos(3),0,0,0]);
    g.Icons.hidecell=hideiconcell;
    end
    
    %--------------------------------------------
    function hidepanels(varargin)
    %----------------------------------------------
    %The function hide all panels hides all panels and displays only image
    %you can use the current tool

    g=varargin{1};
    %First check which state we want to toggle to
    stateandicon = segment('iconson',{'hideresult','hideicons'});
    
%     if get(g.Handles.hideallpanelscheckbox,'value')==1 
%       state2toggle2  = 'off';
%     else
%       state2toggle2  = 'on';
%     end
    
    if stateandicon{1,1}
      %Start setting the visible command
      set(g.Handles.barpanel,'Visible','off')
      set(g.Handles.flowuipanel,'Visible','off')
      %datasetaxeschildren=get(g.Handles.datasetaxes,'children');
      %set(datasetaxeschildren,'Visible','off')
      %set(g.Handles.datasetaxes,'Visible','off')
      set(g.Handles.lvuipanel,'Visible','off')
      set(g.Handles.reportpanel,'Visible','off')
      g.GUISettings.BottomGapHeight = 0.02;
      g.GUISettings.RightGapWidth = 0.02;
    else
      set(g.Handles.barpanel,'Visible','on')
      set(g.Handles.flowuipanel,'Visible','on')
      %datasetaxeschildren=get(g.Handles.datasetaxes,'children');
      %set(datasetaxeschildren,'Visible','on')
      %set(g.Handles.datasetaxes,'Visible','on')
      set(g.Handles.lvuipanel,'Visible','on')
      set(g.Handles.reportpanel,'Visible','on')
      g.GUISettings.BottomGapHeight = 0.133;
      g.GUISettings.RightGapWidth = 0.21;
    end
    
    if stateandicon{2,1}
      currentconfigpos = get(g.Handles.configaxes,'position');
      currenthidepos = get(g.Handles.hideaxes,'position');
      currentdatasetaxespos = get(g.Handles.datasetaxes,'position');
      currentreportpanelpos = get(g.Handles.reportpanel,'position');
      currentthumbnailsliderpos = get(g.Handles.thumbnailslider,'position');
     
      currentreportpanelpos(4)=1-currentconfigpos(4);
      currentconfigpos(2)=1-currentconfigpos(4);
      currenthidepos(2)=1-currentconfigpos(4);
      currentdatasetaxespos(4)=1-currentconfigpos(4);     
      currentthumbnailsliderpos(4) = 1-currentconfigpos(4);
      
      set(g.Handles.configaxes,'position',currentconfigpos)
      set(g.Handles.hideaxes,'position',currenthidepos)
      set(g.Handles.datasetaxes,'position',currentdatasetaxespos)
      set(g.Handles.reportpanel,'position',currentreportpanelpos);
      set(g.Handles.thumbnailslider,'position',currentthumbnailsliderpos);
      g.GUISettings.TopGapHeight = 0.06;  
    set(g.Handles.iconuipanel,'Visible','off')
      
    else
      currentpanelpos = get(g.Handles.iconuipanel,'position');
      currentconfigpos = get(g.Handles.configaxes,'position');
      currenthidepos = get(g.Handles.hideaxes,'position');
      currentdatasetaxespos = get(g.Handles.datasetaxes,'position');
      currentreportpanelpos = get(g.Handles.reportpanel,'position');
      currentthumbnailsliderpos = get(g.Handles.thumbnailslider,'position');
      
      currentdatasetaxespos(4)=currentpanelpos(2);     
      currentreportpanelpos(4)=currentpanelpos(2);
      currentthumbnailsliderpos(4)=currentpanelpos(2);
      currentconfigpos(2)=currentpanelpos(2);
      currenthidepos(2)=currentpanelpos(2);
      set(g.Handles.configaxes,'position',currentconfigpos)
      set(g.Handles.hideaxes,'position',currenthidepos)
      set(g.Handles.datasetaxes,'position',currentdatasetaxespos)
      set(g.Handles.reportpanel,'position',currentreportpanelpos);
      set(g.Handles.thumbnailslider,'position',currentthumbnailsliderpos);
      g.GUISettings.TopGapHeight = 0.135;
      set(g.Handles.iconuipanel,'Visible','on')
    end
    drawfunctions('drawall',g.ViewMatrix(1),g.ViewMatrix(2))
    %updatetool
    
    end
    
    %-----------------------------------
    function togglebuttonLV_Callback(varargin)
    %-----------------------------------
    %Get clicked toggle button set and set configaxes
    g=varargin{1};
    g.Handles.configiconholder.add(g.Icons.lviconcell);
    pos=plotboxpos(g.Handles.configiconholder.axeshandle);
    currentpos=get(g.Handles.configiconholder.axeshandle,'position');
    set(g.Handles.configiconholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);
    g.CurrentTheme='lv';
    
    %If no button is indented in config place holder choose the initial
    %buttons according to discussion with Helen
    %for LV this is the endo pen
    iconcell=g.Icons.lviconcell;
        clickedbutton=0;
    
    for i= 1:g.Handles.configiconholder.numberoficons
      if iconcell{i}.isindented && iconcell{i}.type==1
        clickedbutton=i;
      end
    end
    
    if clickedbutton==0
      indent(g.Handles.configiconholder,'endopen',1)
    else
      indent(g.Handles.configiconholder,iconcell{clickedbutton}.name,1)
    end
    end

    
    %-----------------------------------
    function togglebuttonRV_Callback(varargin)
    %-----------------------------------
    %Get clicked toggle button set and set configaxes
    g=varargin{1};
    
    g.Handles.configiconholder.add(g.Icons.rviconcell); 
    pos=plotboxpos(g.Handles.configiconholder.axeshandle);
    currentpos=get(g.Handles.configiconholder.axeshandle,'position');
    set(g.Handles.configiconholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);     
   
    g.CurrentTheme='rv';
    
    %If no button is indented in config place holder choose the initial
    %buttons according to discussion with Helen
    %for RV this is the endo pen
    iconcell=g.Icons.rviconcell;
        clickedbutton=0;
    
    for i= 1:g.Handles.configiconholder.numberoficons
      if iconcell{i}.isindented && iconcell{i}.type==1
        clickedbutton=i;
      end
    end
    
    if clickedbutton==0
      indent(g.Handles.configiconholder,'rvendopen',1)
    else
      indent(g.Handles.configiconholder,iconcell{clickedbutton}.name,1)
    end
    end
    
    %-----------------------------------
    function iconcell = togglebuttonAnalysis_Callback(varargin)
    %-----------------------------------
    %Get clicked toggle button set configuration
    g=varargin{1};
    
    g.Handles.configiconholder.add(g.Icons.analysisiconcell);

    pos=plotboxpos(g.Handles.configiconholder.axeshandle);
    currentpos=get(g.Handles.configiconholder.axeshandle,'position');
    set(g.Handles.configiconholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);   
    
    g.CurrentTheme='misc';
    
    %If no button is indented in config place holder choose the initial
    %buttons according to discussion with Helen
    %for analysis this is the select tool
    iconcell=g.Icons.analysisiconcell;
    clickedbutton=0;
    
    for i= 1:g.Handles.configiconholder.numberoficons
      if iconcell{i}.isindented && iconcell{i}.type==1
        clickedbutton=i;
      end
    end
    
    if clickedbutton==0
      indent(g.Handles.configiconholder,'select',1)
    else
      indent(g.Handles.configiconholder,iconcell{clickedbutton}.name,1)
    end
    end
    %-----------------------------------
    function iconcell = togglebuttonROIFLOW_Callback(varargin)
    %-----------------------------------
    %Get clicked toggle button set configuration
    g=varargin{1};
    g.Handles.configiconholder.add(g.Icons.roiflowiconcell);
    pos=plotboxpos(g.Handles.configiconholder.axeshandle);
    currentpos=get(g.Handles.configiconholder.axeshandle,'position');
    set(g.Handles.configiconholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);   
    g.CurrentTheme='roi';
    
    
    %If no button is indented in configplaceholder choose the initial
    %buttons according to discussion with Helen
    %for ROIFLOW this is the putroi tool
    iconcell=g.Icons.roiflowiconcell;
    clickedbutton=0;
    
    for i= 1:g.Handles.configiconholder.numberoficons
      if iconcell{i}.isindented && iconcell{i}.type==1
        clickedbutton=i;
      end
    end
    
    if clickedbutton==0
      indent(g.Handles.configiconholder,'putroi',1)
    else
      indent(g.Handles.configiconholder,iconcell{clickedbutton}.name,1)
    end
    end

    %-----------------------------------
    function iconcell = togglebuttonImage_Callback(varargin)
    %-----------------------------------
    %Get clicked toggle button icon configuration
    g=varargin{1};
    g.Handles.configiconholder.add(g.Icons.imageiconcell);
    pos=plotboxpos(g.Handles.configiconholder.axeshandle);
    currentpos=get(g.Handles.configiconholder.axeshandle,'position');
    set(g.Handles.configiconholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);    
   
    g.CurrentTheme='misc';
    %If no button is indented in config place holder choose the initial
    %buttons according to discussion with Helen
    %for image this is the select tool
    iconcell=g.Icons.imageiconcell;
        clickedbutton=0;
    
    for i= 1:g.Handles.configiconholder.numberoficons
      if iconcell{i}.isindented && iconcell{i}.type==1
        clickedbutton=i;
      end
    end
    
    if clickedbutton==0
      indent(g.Handles.configiconholder,'select',1)
    else
      indent(g.Handles.configiconholder,iconcell{clickedbutton}.name,1)
    end
    end

    
    %-----------------------------------
    function iconcell = togglebuttontxmap_Callback(varargin)
    %-----------------------------------
    %Get clicked toggle button icon configuration
    g=varargin{1};
    g.Handles.configiconholder.add(g.Icons.txmapiconcell);
    pos=plotboxpos(g.Handles.configiconholder.axeshandle);
    currentpos=get(g.Handles.configiconholder.axeshandle,'position');
    set(g.Handles.configiconholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);    
   
    g.CurrentTheme='misc';
    
    %If no button is indented in config place holder choose the initial
    %buttons according to discussion with Helen
    %for image this is the select tool
    iconcell=g.Icons.txmapiconcell;
    clickedbutton=0;
    
    for i= 1:g.Handles.configiconholder.numberoficons
      if iconcell{i}.isindented && iconcell{i}.type==1
        clickedbutton=i;
      end
    end
    
    if clickedbutton==0
      indent(g.Handles.configiconholder,'select',1)
    else
      indent(g.Handles.configiconholder,iconcell{clickedbutton}.name,1)
    end
    end
    
    %-----------------------------------
    function iconcell = togglebutton3dp_Callback(varargin)
    %-----------------------------------
    %Get clicked toggle button icon configuration
    g=varargin{1};
    g.Handles.configiconholder.add(g.Icons.printiconcell);
    pos=plotboxpos(g.Handles.configiconholder.axeshandle);
    currentpos=get(g.Handles.configiconholder.axeshandle,'position');
    set(g.Handles.configiconholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);    
   
    g.CurrentTheme='misc';
    %If no button is indented in config place holder choose the initial
    %buttons according to discussion with Helen
    %for image this is the select tool
    iconcell=g.Icons.printiconcell;
    clickedbutton=0;
    
    for i= 1:g.Handles.configiconholder.numberoficons
      if iconcell{i}.isindented && iconcell{i}.type==1
        clickedbutton=i;
      end
    end
    
    if clickedbutton==0
      indent(g.Handles.configiconholder,'select',1)
    else
      indent(g.Handles.configiconholder,iconcell{clickedbutton}.name,1)
    end
    end
    
    %-----------------------------------
    function iconcell = togglebuttonVia_Callback(varargin)
    %-----------------------------------
    %Get clicked toggle button icon configuration
    g=varargin{1};
    g.Handles.configiconholder.add(g.Icons.viabilityiconcell);
    pos=plotboxpos(g.Handles.configiconholder.axeshandle);
    currentpos=get(g.Handles.configiconholder.axeshandle,'position');
    set(g.Handles.configiconholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);
    g.CurrentTheme='scar';
    
    %If no button is indented in config place holder choose the initial
    %buttons according to discussion with Helen
    %for viability this is the select tool
    iconcell=g.Icons.viabilityiconcell;
    clickedbutton=0;
    for i= 1:g.Handles.configiconholder.numberoficons
      if iconcell{i}.isindented && iconcell{i}.type==1
        clickedbutton=i;
      end
      
      if strcmp(iconcell{i}.name,'select')
        ind=i;
      end
    end
    
    if clickedbutton==0
      iconcell{ind}.cdataDisplay=iconcell{ind}.cdataIndent;
      iconcell{ind}.isindented=1;
      feval(iconcell{ind}.execute);
    else
      iconcell{clickedbutton}.cdataDisplay=iconcell{clickedbutton}.cdataIndent;
      iconcell{clickedbutton}.isindented=1;
      feval(iconcell{clickedbutton}.execute);
    end
      g.Handles.configiconholder.render;
    end
    
    %-------------------------------------------------
      function toggleplaceholdermotion(varargin)
        %-------------------------------------------------
        g=varargin{1};
        
        handleAddress=hittest(g.fig);
        if ~isempty(get(g.Handles.configaxes,'Children'))
          if isequal(handleAddress,g.Handles.configiconholder.imagehandle)
            g.Handles.configiconholder.motion
          else
            g.Handles.configiconholder.notover
          end
        end
        
        if isequal(handleAddress,g.Handles.permanenticonholder.imagehandle)
          g.Handles.permanenticonholder.motion
        else
          g.Handles.permanenticonholder.notover
        end
        
        if isequal(handleAddress,g.Handles.toggleiconholder.imagehandle)
          g.Handles.toggleiconholder.motion
        else
          g.Handles.toggleiconholder.notover
        end
        
        if ~isempty(get(g.Handles.hideaxes,'Children'))
        if isequal(handleAddress,g.Handles.hideiconholder.imagehandle)
          g.Handles.hideiconholder.motion
        else
          g.Handles.hideiconholder.notover
        end
        end
      end
      
%     %----------------------------------------
%     function initpermanentplaceholder(varargin)
%       %--------------------------------------
%     g=varargin{1};
%       
%     iconcell=cell(1,1);
%     iconcell{1,1}=myicon('database',g.Handles.permanenticonholder,g.Icons.config.database,'Open patient database',@() segment('fileopen_Callback'),0);
%     iconcell{1,end+1}=myicon('databaseadd',g.Handles.permanenticonholder,g.Icons.config.databaseadd,'Save to patient database',... 
%       @() filemenu('savetopatientdatabase_Callback'),0);
% %     iconcell{1,end+1}=myicon(g.Handles.permanenticonholder,g.Icons.config.connect,'Open PACS connection','pacs(''init_Callback'')',0);
% %     iconcell{1,end+1}=myicon(g.Handles.permanenticonholder,g.Icons.config.connectadd,'Save to PACS','filemenu(''savetopacs_Callback'')',0);
%     iconcell{1,end+1}=myicon('closeall',g.Handles.permanenticonholder,g.Icons.config.closeall,'Close all image stacks',@() segment('filecloseall_Callback'),0);
%     
%     iconcell{1,end+1}=myicon('panel1',g.Handles.permanenticonholder,g.Icons.config.panel1,'View one image panel',@() drawfunctions('drawall',1),1,1);
%     iconcell{1,end+1}=myicon('panel2',g.Handles.permanenticonholder,g.Icons.config.panel2,'View two image panels',@() drawfunctions('drawall',1,2),1,1);
%     iconcell{1,end+1}=myicon('panel2x1',g.Handles.permanenticonholder,g.Icons.config.panel2x1,'View two image panels',@() drawfunctions('drawall',2,1),1,1);
%     iconcell{1,end+1}=myicon('panel3x1',g.Handles.permanenticonholder,g.Icons.config.panel3x1,'View three image panels',@() drawfunctions('drawall',3,1),1,1);
%     iconcell{1,end+1}=myicon('panel3',g.Handles.permanenticonholder,g.Icons.config.panel3,'View three image panels',@() drawfunctions('drawall',1,3),1,1);
%     iconcell{1,end+1}=myicon('panel4',g.Handles.permanenticonholder,g.Icons.config.panel4,'View four image panels',@() drawfunctions('drawall',2,2),1,1);
%     iconcell{1,end+1}=myicon('panel6',g.Handles.permanenticonholder,g.Icons.config.panel6,'View six image panels',@() drawfunctions('drawall',6),1,1);
%     iconcell{1,end+1}=myicon('saveview',g.Handles.permanenticonholder,g.Icons.config.saveview,'Save view',@() segmentview,0);
%     
%     iconcell{1,end+1}=myicon('viewone',g.Handles.permanenticonholder,g.Icons.config.viewone,'View one slice',@() segment('viewimage_Callback','one'),1,2);
%     iconcell{1,end+1}=myicon('viewall',g.Handles.permanenticonholder,g.Icons.config.viewall,'View all slices',@() segment('viewimage_Callback','montage'),1,2);
%     iconcell{1,end+1}=myicon('viewrow',g.Handles.permanenticonholder,g.Icons.config.viewrow,'View all slices in 2 rows',@() segment('viewimage_Callback','montagerow'),1,2);
%     
%     iconcell{1,end+1}=myicon('undo',g.Handles.permanenticonholder,g.Icons.config.undo,'Undo last operation',@() tools('undosegmentation_Callback'),0);
%     iconcell{1,end+1}=myicon('refresh',g.Handles.permanenticonholder,g.Icons.config.refresh,'Refresh image view',@() segment('viewrefreshall_Callback'),0);
%     
%     iconcell{1,end+1}=myicon('play',g.Handles.permanenticonholder,g.Icons.config.play,'Play movie of all image stacks',@() segment('playall_Callback'),2);
%     iconcell{1,end+1}=myicon('next',g.Handles.permanenticonholder,g.Icons.config.next,'Next time frame for all image stacks',@() segment('nextallframe_Callback'),0);
%     iconcell{1,end+1}=myicon('prev',g.Handles.permanenticonholder,g.Icons.config.prev,'Previous time frame for all image stacks',@() segment('previousallframe_Callback'),0);
%     iconcell{1,end+1}=myicon('faster',g.Handles.permanenticonholder,g.Icons.config.faster,'Faster frame rate',@() segment('fasterframerate_Callback'),0);
%     iconcell{1,end+1}=myicon('slower',g.Handles.permanenticonholder,g.Icons.config.slower,'Slower frame rate',@() segment('slowerframerate_Callback'),0);
%     
%     iconcell{1,end+1}=myicon('hideall',g.Handles.permanenticonholder,g.Icons.config.hideall,'Hide all overlays (segmentation, point, text,...)',@() segment('viewhideall_Callback'),2);
%     iconcell{1,end+1}=myicon('clearall',g.Handles.permanenticonholder,g.Icons.config.clearall,'Clear all segmentation in current image stack',@() segment('segmentclearall_Callback'),0);
%     iconcell{1,end+1}=myicon('clearalledes',g.Handles.permanenticonholder,g.Icons.config.clearalledes,'Clear all segmentation except in ED and ES',@() segment('segmentclearalllvbutsystolediastole_Callback'),0); 
%     iconcell{1,end+1}=myicon('zoomin',g.Handles.permanenticonholder,g.Icons.config.zoomin,'Zoom in',@() segment('viewzoomin_Callback'),0);
%     iconcell{1,end+1}=myicon('zoomout',g.Handles.permanenticonholder,g.Icons.config.zoomout,'Zoom out',@() segment('viewzoomout_Callback'),0);
%     iconcell{1,end+1}=myicon('colorbar',g.Handles.permanenticonholder,g.Icons.config.colorbar,'Hide colorbar',@() segment('viewhidecolorbar_Callback'),2);
%     iconcell{1,end+1}=myicon('viewpixels',g.Handles.permanenticonholder,g.Icons.config.viewpixels,'Show image pixels',@() segment('viewinterp_Callback'),2);
%     iconcell{1,end+1}=myicon('reportsheet',g.Handles.permanenticonholder,g.Icons.config.reportsheet,'Open Report sheet generation',@() reporter.reportsheet,0);
%     iconcell{1,end+1}=myicon('savescreen',g.Handles.permanenticonholder,g.Icons.config.savescreen,'Save screen shot',@() export('screenshot_Callback'),0);
%     
%     iconcell{1,end+1}=myicon('settingsgeneral',g.Handles.permanenticonholder,g.Icons.config.settingsgeneral,'Set general preferences',@() segpref,0);
%     iconcell{1,end+1}=myicon('settingssystem',g.Handles.permanenticonholder,g.Icons.config.settingssystem,'Set patient database preferences',@() segpref('advancedsettings_Callback'),0);
%     iconcell{1,end+1}=myicon('settingspacs',g.Handles.permanenticonholder,g.Icons.config.settingspacs,'Set PACS connection preferences',@() pacspref,0);
%     g.Handles.permanenticonholder.add(iconcell);
%     
%     pos=plotboxpos(g.Handles.permanenticonholder.axeshandle);
%     currentpos=get(g.Handles.permanenticonholder.axeshandle,'position');
%     set(g.Handles.permanenticonholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);
%     set(g.Handles.iconuipanel,'visible','on')
%     
%     end
    
    %------------------------------------
    function startmodeplaceholders(varargin)
    %--------------------------------------
    % undents and disables all icons in toggleplaceholder. This is protocoll when no images
    %are loaded
    
    g=varargin{1};

    %empty config and permanent iconplaceholder
    if isfield(g.Handles,'configiconholder')
       delete(get(g.Handles.configaxes,'Children'))
       g.Handles.configiconholder.cdata=[];
    end
    
    if isfield(g.Handles,'hideiconholder')
       delete(get(g.Handles.hideaxes,'Children'))
       g.Handles.hideiconholder.cdata=[];
    end
    
    g.initpermanentplaceholder;
    iconcell=g.Handles.permanenticonholder.iconCell;
    
    icons2disable_ind=[2: numel(iconcell)-3];
    
    for i = icons2disable_ind
      iconcell{i}.disable;
    end
   
    iconcell=g.Handles.toggleiconholder.iconCell;
    
    for i = 1:numel(iconcell)
      iconcell{i}.disable;
    end
    g.Handles.toggleiconholder.disablepad=1;
    g.Handles.toggleiconholder.render;
    
           iconcell=g.Handles.toggleiconholder.iconCell;
    
    for i = 1:numel(iconcell)
      iconcell{i}.disable;
    end
    g.Handles.toggleiconholder.disablepad=1;
    g.Handles.toggleiconholder.render;
   % set(g.Handles.hideallpanelsradiobutton,'visible','off');

    
%     %Disable toggle buttons
%     set([g.Handles.viabilitypushbutton,g.Handles.analysispushbutton,...
%       g.Handles.leftventriclepushbutton,g.Handles.rightventriclepushbutton,g.Handles.roipushbutton,...
%       g.Handles.imagepushbutton],'Enable','off')
    
    end
    
    %---------------------------------------------
    function setviewbuttons(varargin)
    %---------------------------------------------
    g=varargin{1};
    
    if nargin~=1
      runbuttons=varargin{2};
    else
      runbuttons=1;
    end
    
    %First the panel group is found
    str=sprintf('%d%d',g.ViewMatrix(1),g.ViewMatrix(2));
    
    namelist={'panel1','panel2','panel2x1','panel3','panel3x1','panel4','panel6'};
    name=[];
    switch str
      case '11'
        name='panel1';
      case '12'
        name='panel2';
      case '21'
        name='panel2x1';
      case '13'
        name='panel3';
      case '31'
        name='panel3x1';
      case '22'
        name='panel4';
      case '23'
        name='panel6';
        %case '24'
        %name='panel8';
    end
    
    %if~isempty(name)
    ind=[];
    iconcell=g.Handles.permanenticonholder.iconCell;
    for i= 1:g.Handles.permanenticonholder.numberoficons
      if strcmp(iconcell{i}.name,name)
        ind=i;
        %shouldnt do anything if the icons are disabled i.e we are in start
        %mode
        if ~iconcell{i}.enabled
          return;
        end
      end
      if any(strcmp(iconcell{i}.name,namelist))
        iconcell{i}.undent
      end
    end
    
    if ~isempty(ind)
      iconcell{ind}.cdataDisplay=iconcell{ind}.cdataIndent;
      iconcell{ind}.isindented=1;
      
      if runbuttons
        feval(iconcell{ind}.execute);
      end
    end
    %end
    %if viewpanels isempty assume viewone
    %     if isempty(g.ViewPanelsType)
    %       name='viewone';
    %     else
    ind=[];
    %Then view mode
    switch g.ViewPanelsType{g.CurrentPanel}
      case 'montage'
        name='viewall';
      case 'one'
        name='viewone';
      case 'montagerow'
        name='viewrow';
    end
    %     end
    iconcell=g.Handles.permanenticonholder.iconCell;
    for i= 1:g.Handles.permanenticonholder.numberoficons
      if strcmp(iconcell{i}.name,name)
        ind=i;
      end
      if any(strcmp(iconcell{i}.name,{'viewall','viewone','viewrow'}))
        iconcell{i}.undent
      end
    end
    
    if ~isempty(ind)
      iconcell{ind}.cdataDisplay=iconcell{ind}.cdataIndent;
      iconcell{ind}.isindented=1;
      
      if runbuttons
        feval(iconcell{ind}.execute);
      end
    end
    g.Handles.permanenticonholder.render;
    end
    
    
    %------------------------------------
    function dataloadedplaceholders(varargin)
    %--------------------------------------
    % When data is loaded the placeholders are all enabled
    global SET
    
    g=varargin{1};
    
    updatetool('select');

    g.initpermanentplaceholder;
    g.inithideplaceholder;
    
    g.setviewbuttons(1);

    iconcell=g.Handles.toggleiconholder.iconCell;
    for i = 1:numel(iconcell)
      iconcell{i}.enable;
    end
    
    g.Handles.toggleiconholder.disablepad=0;
    
    %asserts that first button is toggle button.
    iconcell{1}.cdataDisplay=iconcell{1}.cdataIndent;
    iconcell{1}.isindented=1;
    feval(iconcell{1}.execute);
    g.Handles.toggleiconholder.clickedicon=iconcell{1};
    g.Handles.toggleiconholder.render;
    
    %Click the single frame mode button
    stateandicon = segment('iconson','singleframemode');
    icon = stateandicon{1,2};
    icon.cdataDisplay=icon.cdataIndent;
    icon.isindented=1;
    g.ThisFrameOnly = true;
    
    %synchronize click
    indent(g.Handles.hideiconholder,'synchronize',1)

    g.Handles.hideiconholder.render
    %g.thisframeonly(1)
    if ~isempty(SET(1).FileName)
      str = sprintf('%s\t%s\t%s\t%s',datestr(now,'yyyy-mm-dd HH:MM'),getenv('USERNAME'),sprintf('Loaded file %s',SET(1).FileName),SET(1).PatientInfo.ID);
    else
      str = sprintf('%s\t%s\t%s\t%s',datestr(now,'yyyy-mm-dd HH:MM'),getenv('USERNAME'),sprintf('Loaded file %s',SET(1).PathName),SET(1).PatientInfo.ID);
    end
    g.adduserevent(str);
    %g.adduserevent([' Loaded file: ', SET(1).FileName])
    %g.adduserevent(['Time:', datestr(now,'yyyymmddHHMMSS')])
    %g.adduserevent([' Patient name: ', SET(1).PatientInfo.Name])
    %g.adduserevent([' Patient ID: ', SET(1).PatientInfo.ID])
    
    end
    
    
    %---------------
    function init(g)
    %---------------
    %Create GUI
    g.initmaingui;
    
    loadpreferences;
    
    setuplicense; %This is done to set new license system in use
    
    if isequal(13,getmodule(12,'Ev',[],true))
      useragreement
      if ~ishandle(g.fig)
        %The software has been terminated
        return
      end
    end
    
    if g.Pref.UserLogging
      answer = yesno(sprintf('You are currently logged in as %s. Do you wish to proceed?',getenv('USERNAME')));
      
      if answer == 0
        filemenu('quit_Callback','Aborting initialization.');
        return
      else
        
        str = sprintf('%s\t%s\t%s\t%s',datestr(now,'yyyy-mm-dd HH:MM'),getenv('USERNAME'),'Software started.','-');
        fail = g.adduserevent(str);
        if fail
          
        end
        %str = sprintf(sprintf('Software started by %s',getenv('USERNAME'))
        %g.adduserevent(sprintf('Software started by %s',getenv('USERNAME'))) %starts user logging if preferences says so
        %g.adduserevent(['Time:', datestr(now,'yyyymmddHHMMSS')])
      end
    end
    
    setappdata(0,'UseNativeSystemDialogs',false); %this controls style of uigetdir

    % Generate a structure of handles and user data to pass to callbacks.
    %handles = guihandles(g.fig);

    % Check for new versions of software
    g.checkversion;
    
    %g.Handles=handles;
    g.Handles.imageaxes = [];
    g.Handles.boxaxes = [];
    g.imagefig = g.fig;
    
    segment('resetpreview');
    g.Run = false;
    g.Record = false;
    g.LastPointer = 'arrow';
    g.LastPointerShapeCData = [];
    g.VisibleThumbnails = [];    
    
    %g.ExcludePapilars = g.Pref.ExcludePapilars;
    %g.UseLight = g.Pref.UseLight;

    %Add some graphics
    %Install icons
    try
      load('newicons.mat','newicons')
    catch me
      myfailed('Critical error: Could not read icons.',g);
      mydispexception(me);
    end;
    
    %not used but idea is to speed up searching for active icons
    g.Icons.permanent = newicons;
    g.Icons.config = newicons;
    g.Icons.toggleicons = newicons;
    
    %added by eriks
    set(g.fig,'WindowScrollWheelFcn',@(h,e) segment('changewheel_Callback',h,e));

    %Create icon cell in which we place icons for desired appearance
    g.Handles.toggleiconholder = myiconplaceholder(g.Handles.ribbonaxes,1);
    g.Handles.permanenticonholder = myiconplaceholder(g.Handles.permanentaxes);
    g.Handles.configiconholder = myiconplaceholder(g.Handles.configaxes);
    g.Handles.hideiconholder = myiconplaceholder(g.Handles.hideaxes,0,2);%2 is that it drags to right
    
    %we must manage multiple axes so we need a hittest check for all
    %available axeses
    set(g.fig,'WindowButtonMotionFcn',@g.toggleplaceholdermotion);
   
    g.initribbonplaceholder;
    g.initconfigplaceholder;
    g.inithideplaceholder;
    g.startmodeplaceholders 
    g.inittoolbar;

    %--- Checks all parameterfiles, and fills listbox
    if isdeployed
      pathname = pwd;
    else
      temp = mfilename('fullpath');
      pathname = fileparts(temp);
    end;

    f = dir([pathname filesep '*.par']);
    names = cell(1,length(f));
    [names{:}] = deal(f(:).name);

    %Take away extension
    for loop=1:length(names)
      ind = find(names{loop}=='_');
      ind = ind(end);
      names{loop} = names{loop}(1:(ind-1));
    end;

    %Remove duplicates
    names = unique(names);

    %Load the files to check name
    stri = cell(1,length(names));

    for loop=1:length(names)
      try
        load([names{loop} '_endo.par'],'-mat');
      catch %#ok<CTCH>
        load([names{loop} '_epi.par'],'-mat');
      end;

      %Remove 'endo'
      tempstri = p.name;
      ind = findstr(tempstri,' endo');
      if not(isempty(ind))
        tempstri = tempstri(1:(ind-1));
      end;

      %Remove 'epi'
      ind = findstr(tempstri,' epi');
      if not(isempty(ind))
        tempstri = tempstri(1:(ind-1));
      end;

      %Store it
      stri{loop} = tempstri;
    end;

    g.ImagingTechniques = names;
    g.ImagingTechniquesFullNames = stri;

    g.NeedToSave = false;
    set([...
      g.Handles.filesaveicon ...
      g.Handles.pacsaddicon ...
      g.Handles.databaseaddicon],'enable','off');

    %Close restores defaults.
    g.filecloseall_Callback(true); %silent

    handles = [...
      g.Handles.filemenu ...
      g.Handles.utilitymenu ...
      g.Handles.pluginmenu ...
      g.Handles.helpmenu]; %...
    %g.Handles.fileopenicon ...
      %g.Handles.databaseicon ...
      %g.Handles.importfromcdicon ...
      %g.Handles.pacsicon];

    set(handles,'enable','off');

    g.versionhello;

    if ~g.Locked
      try
        set(handles,'enable','on');
      catch
      end;
    end;

    %This callback is set here to avoid error messages at startup
    set(g.GUI.Segment.fig,'ResizeFcn','segment(''mainresize_Callback'')');
    
    g.updatetitle;

    if ~g.Locked
      
      % Start sectramodule
      sectra('initsectra');
      %Start externalpacsmodule
      externalpacs('initexternalpacs');
    end;
    
    segpref('setbackgroundcolor',g.Pref.GUIBackgroundColor)
    %g.init_userlogging
    end
    
    %----------------------
    function initmaingui(g)
    %----------------------
    %Initiates main GUI. Overloaded in most other GUI's 
      
    g.GUI.Segment=mygui('segment.fig');%Load standard Segment GUI
    g.fig = g.GUI.Segment.fig;
    g.Handles = guihandles(g.fig);
		g.Handles = killhandles(g.fig,{'segment.fig'});
    
    colormap(gray(255)); %Moved this to here.    
    end

    %----------------------------
    function init_userlogging(g)
    %----------------------------
      val = mygetvalue(g.PrefHandlesAdvanced.userloggingcheckbox);
      
      if val == 1
        %Initiate log file.
        pathname =  g.Pref.UserLogPath;
        %g.UserFile = [pathname filesep sprintf('segmentuserlog.log')];
        fid = fopen(pathname,'a');
        if isequal(fid,-1)
          myfailed(dprintf('Could not write to user log path %s \n User logging is turned off.',pathname));
          g.Pref.UserLogPath=[];
          g.Pref.UserLogging = 0;
          %if isopen('segpref.fig')
            set(g.PrefHandlesAdvanced.userloggingcheckbox,'Value',0);
            set(g.PrefHandlesAdvanced.setuserlogpathpushbutton,'Enable','off');
          %end
          return;
        end;
        fclose(fid);
        adduserevent(g,sprintf('%s\t%s\t%s\t%s','DATE','USER','EVENT','PATIENT ID'));
        adduserevent(g,sprintf('%s\t%s\t%s\t%s',datestr(now,'yyyy-mm-dd HH:MM'),getenv('USERNAME'),'Started user log','-'));
      end
    end
    
    %---------------------------------
    function fail = adduserevent(g,str)
    %---------------------------------
    fail = 0;
    if g.Pref.UserLogging
      fid = fopen(g.Pref.UserLogPath,'a');
      if fid == -1
        fail=1;
        myfailed(dprintf('User logging turned off.\n Since unable to write to .txt file: \n %s.',g.Pref.UserLogPath));
        g.Pref.UserLogPath=[];
        g.Pref.UserLogging = 0;
        return
      end
      %mydisp(str)
      fprintf(fid,'%s\r\n',str);
      fclose(fid);
    end
    end
    
    %----------------------
    function initLogFile(g)
    %----------------------
    %Initiate log file. Overloaded in most other GUI's
    pathname = getpreferencespath;
    g.LogFile = [pathname filesep sprintf('segmentlog_%s.log',datestr(now,'yyyymmddHHMMSS'))];
    fid = fopen(g.LogFile,'w');
    if isequal(fid,-1)
      myfailed(dprintf('Could not create .log file %s.',g.LogFile));
    else
      fclose(fid); %Close the file, to make it empty
      g.startlog(g.LogFile); %Start diary.
    end;
    end
    
    %--------------------------
    function inittoolbar(g) 
    %--------------------------
    %Initiates toolbars. Overloaded in CVQgui, Segment CMR gui and RVQgui
    %g.initmaintoolbar;
    %segment('initviewtoolbar');
    segment('initmenu');
    end
    
    %--------------------------
    function initmaintoolbar(g)
    %--------------------------
    %Initiate main toolbar. Overloaded in any GUI with a toolbar
    end
     
    %----------------------
    function updatetitle(g)
    %----------------------
    %Updates titleline of main GUI. Overloaded in CVQgui and RVQgui
    global SET NO

    if g.Silent
      return;
    end;

    if (NO<1)||(NO>length(SET))
      no = 1;
    else
      no = NO;
    end;

    licensename = getmodule('getlicensename');
    if (isempty(licensename)|isnan(licensename)) %#ok<OR2>
      name = '(Academic Research only version)';
    else
      name = ['(' deblank(licensename') ')'];
    end;
    titlestring = strcat([g.ProgramName ...
      ' v' g.ProgramVersion ...
      ' ' name]);

    if g.Pref.AnonymMode
      set(g.imagefig,'name',...
        sprintf('%s',titlestring));
    else
      if isempty(SET)
        location= '';
      else
        location = sprintf('/ %s %s',SET(no).PatientInfo.Name,SET(no).PatientInfo.ID);
      end
      set(g.imagefig,'name',...
        sprintf('%s %s',titlestring,location));
    end;
    
    end
    
    
    %-------------------------------------
    function orthoview_Buttondown(g,panel)
    %-------------------------------------
    %Switch to slices corresponding to location of user click
    global NO
    
    segment_main('normal_Buttondown',panel);
    [y,x] = mygetcurrentpoint(g.Handles.imageaxes(panel));
    silent = true;
    switch g.ViewPanelsType{panel}
      case 'ortho'
        segment('switchtolongaxisslice',round(x),'HLA',silent);
        segment('switchtolongaxisslice',round(y),'VLA');
      case 'hla'
        segment('switchtolongaxisslice',round(y),'VLA',silent);
        segment('switchtoslice',round(x));
      case 'vla'
        segment('switchtolongaxisslice',round(y),'HLA',silent);
        segment('switchtoslice',round(x));
      case 'gla'
        [vslice,hslice] = calcfunctions('gla2sax',x,y,NO);
        segment('switchtolongaxisslice',round(vslice),'HLA',silent);
        segment('switchtolongaxisslice',round(hslice),'VLA',silent);
        segment('switchtoslice',round(x));
    end
    end
    
    %-------------------------------------
    function glarotatehandle_Buttondown(g)
    %-------------------------------------
    %Rotate GLA view using handle on intersection line
    global SET NO
    
    segment_main('normal_Buttondown',1);
    glaangle = SET(NO).GLA.angle;
    panel = g.CurrentPanel;
    
    xhan = get(g.Handles.glarotatehandle(panel),'XData');
    yhan = get(g.Handles.glarotatehandle(panel),'YData');
    xhan = xhan(2);
    yhan = yhan(2);
    set(g.Handles.glarotatehandle(panel), ...
      'XData',xhan+1000*cos(glaangle)*[-1 0 1], ...
      'YData',yhan+1000*sin(glaangle)*[-1 0 1], ...
      'LineStyle','-');
    set(g.fig,'WindowButtonMotionFcn','segment(''glarotatehandle_Motion'')');
    set(g.fig,'WindowButtonUpFcn','segment(''glarotatehandle_Buttonup'')');
    end
    
    %---------------------------------
    function glarotatehandle_Motion(g)
    %---------------------------------
    %Motion fcn for using GLA view rotation handle
    global SET NO
    panel = g.CurrentPanel;
    [x,y] = mygetcurrentpoint(g.Handles.imageaxes(panel));
    glaangle = mod(atan2(SET(NO).ResolutionX*(y-SET(NO).HLA.slice), ...
      SET(NO).ResolutionY*(x-SET(NO).VLA.slice))+pi/2,pi)-pi/2;
    set(g.Handles.glarotatehandle(panel), ...
      'XData',x+1000/SET(NO).ResolutionY*cos(glaangle)*[1 0 -1], ...
      'YData',y+1000/SET(NO).ResolutionX*sin(glaangle)*[1 0 -1]);
    SET(NO).GLA.angle = glaangle;
    end
    
    %-----------------------------------
    function glarotatehandle_Buttonup(g)
    %-----------------------------------
    %Update view according to new rotation
    global NO
    set(g.fig,'WindowButtonMotionFcn',[], ...
      'WindowButtonUpFcn',[]);
    set(g.Handles.glarotatehandle(g.CurrentPanel),'LineStyle','none');
    segment_main('updateoneim',NO);
    drawfunctions('updatenopanels',NO);
    drawfunctions('drawintersections');
    end
    
    %--------------------------
    function dispwelcometext(g)
    %--------------------------
    %Displays welcome text. Overloaded in CVQgui, SegmentCMR GUI
    mydisp('Welcome to Segment, Software for Cardiac Image Analysis');
    mydisp(' ');
    mydisp(['Software version ' g.ProgramVersion ' Copyright Einar Heiberg, Medviso AB (http://www.medviso.com)']);
    mydisp(' ');
    mydisp('In this window system and error messages will appear.');
    mydisp('Do not close this window.');
    mydisp(' ');
    end
    
    %---------------------------
    function loadguipositions(g)
    %---------------------------
    %Load prestored positions of GUI's.

    pathname = getpreferencespath;

    if ispc
      if not(exist([pathname filesep '.segment_guipositions.mat'],'file'))
        %mywarning('Moving system preferences file to new user location.');
        mydisp('Moving system preferences file to new user location.');
      end;
    end;

    %Try to find a more suitable place for the file
    %try
      load([pathname filesep '.segment_guipositions.mat'],'GUIPositions');
      %names={GUIPositions.FileName};  
      %ind=find(strcmp(names,'monitor'));
     
      mp = get(0,'MonitorPositions');
%       ss = get(0,'ScreenSize');
%       for mploop = 1:4, sameposition(mploop) = mp(1,mploop) == ss(1,mploop); end
%       if sum(sameposition) < 4
%           mp = ss;
%       end
      %return to default and introduce new guipositions field monitor
      if ~strcmp(GUIPositions(1).FileName,'monitor')
        DATA.GUIPositions=[];
        DATA.GUIPositions(1).FileName=monitor;
        DATA.GUIPositions(1).Position=mp;
        return;
      end
      
      prevmp=GUIPositions(1).Position;
      
      if all(size(mp)==size(prevmp)) && all(mp==prevmp)
        g.GUIPositions = GUIPositions;  %#ok<CPROP> 
      else
        DATA.GUIPositions=[];
        return
      end
%       
%       if monitor
%       %check if only one monitor and positions outside this monitor
%       %then change guipositions to be inside the monitor
%       monitorpositions = get(0,'MonitorPositions');
%       numbermonitors = size(monitorpositions,1);
%       if numbermonitors == 1
%         for loop=1:length(g.GUIPositions)
%           pos = g.GUIPositions(loop).Position;
%           if pos(1)>=1 || pos(1)<0 || pos(2)>=1 || pos(2)<0
%             pos(1)=pos(1)-floor(pos(1));
%             pos(2)=pos(2)-floor(pos(2));
%             if pos(3)>1
%               pos(3)=1;
%             end
%             if pos(4)>1
%               pos(4)=1;
%             end
%             g.GUIPositions(loop).Position=pos;
%           end
%         end
%       end
%     catch me
%       mydispexception(me);
%       mydisp('Could not read guipositions. Probably a new installation.');
%     end;
%     end
    end
    
    %------------------------------------
    function enableplaceholders(varargin)
    %--------------------------------------
    % Enables all icons in toggleplaceholder. This is protocoll when images
    %are loaded
    
    g=varargin{1};
    
    iconcell=g.Handles.toggleiconplaceholder.iconcell;
    
    for i = 1: numel(iconcell)
      iconcell{i}.enable;
    end
    end
    
    %------------------------------------
    function filecloseall_Callback(g,silent) %#ok<INUSD>
    %------------------------------------
    %Deletes all image stacks and resets Segment to original state.
    %Overloaded in CVQgui and RVQgui
    global SET NO 

    if (g.NeedToSave)&&(nargin==1)
      if ~yesno('Unsaved changes: are you sure you want to close all image stacks?',[],g.GUI.Segment);
        return;
      end;
    end;
    
    %Undent and disable all the toggle buttons

    g.startmodeplaceholders


    flushlog;

    g.Silent = false;
    g.LastSaved = now;
    g.ViewPanels = [];
    g.ViewPanelsType = {};
    g.ViewPanelsMatrix = {};
    g.ViewIM = {};
    g.Overlay = [];
    g.ViewMatrix = [];
    g.CurrentPanel = 1;
    g.Undo = [];
    g.UndoN = [];
    g.VisibleThumbnails = [];
    g.DATASETPREVIEW = [];
    
    %Close all associated GUIs
    g.closeallnoguis(1:numel(SET));

    %clear functions; %EH: Removed 2010-08-22

    %Delete window
    try
      delete(g.Handles.imageaxes);
      delete(g.Handles.boxaxes);
    catch %#ok<CTCH>
    end;

    %Disable non-valid options.
     set([... 
       g.Handles.measuremenu ...
       g.Handles.toolsmenu ...
       g.Handles.roimenu ...
       g.Handles.editmenu ...
       g.Handles.exportmenu ...
       g.Handles.viewmenu ....
       g.Handles.segmentmenu ...
       g.Handles.mrmenu ...
       g.Handles.viabilitymenu ...
       g.Handles.flowmenu ...
       g.Handles.ctmenu ...
       g.Handles.spectmenu ...
       g.Handles.lvmenu ...
       g.Handles.rvmenu ...
       g.Handles.analysismenu ...
       g.Handles.fileloadsegmentationmenu ...
       g.Handles.fileloadnextmenu ...
       g.Handles.filesavesegmentationmenu ...
       g.Handles.filesavecurrentmenu ...
       g.Handles.filesaveallmenu ...
       g.Handles.filesaveallasmenu ...
       g.Handles.filesavetodatabasemenu ...
       g.Handles.filesavetopacsmenu ...
       g.Handles.filesavesegdicom ...
       g.Handles.fileclosecurrentimagestack ...
       g.Handles.filecloseallmenu ...
       g.Handles.fileclosemultiplemenu ...
       g.Handles.filesavesubmenu ...
       g.Handles.thumbnailslider...
       ],'enable','off');

    set(g.imagefig,'keypressfcn','','name',...
      sprintf(strcat([g.ProgramName ' v' g.ProgramVersion])));

    %clear the measurementresultaxes before closing
    g.measurementreportclear;
    
    set([...
      g.Handles.datasetaxes ...
      g.Handles.reportpanel ...
      g.Handles.barpanel ...
      g.Handles.flowaxes ... 
      g.Handles.timebaraxes ... %g.Handles.excludepapilarscheckbox ... %g.Handles.uselightcheckbox ...g.Handles.thistimeframeonlycheckbox ... 
      ],...
      'visible','off');

     set([...
       g.Handles.thumbnailslider...
       ],'visible','off');
     
     %kill of children in timebar and volumeaxes.
      cla(g.Handles.timebaraxes);
      g.Handles.timebar = [];
      g.Handles.timebaraxeshelpers=[];
      g.Handles.edtimebartext=[];
      g.Handles.estimebartext=[];
      g.Handles.edtimebarline=[];
      g.Handles.estimebarline=[];
      
      cla(g.Handles.volumeaxes);
      g.Handles.timebarlv = [];
      g.Handles.volumecurve = [];
      g.Handles.masscurve = [];
      g.Handles.rvvcurve = [];
      g.Handles.volumeaxeshelpers = [];
      g.Handles.estext=[];
      g.Handles.edtext=[];
      g.Handles.esline=[];
      g.Handles.edline=[];
      g.Handles.lvmtext=[];

    if g.DataLoaded
      set(get(g.Handles.volumeaxes,'Children'),'visible','off');
      set(get(g.Handles.datasetaxes,'Children'),'visible','off');
      set(get(g.Handles.flowaxes,'Children'),'visible','off');
      set(get(g.Handles.timebaraxes,'Children'),'visible','off');
    end;

    g.ThisFrameOnly = true;
    g.Synchronize=0;
    g.GUISettings.AxesColor=[1,1,1];%g.GUISettings.AxesColorDefault;

    %Close potentially open figs. More might be added to this list
    if isa(g.GUI.GreyZoneHistAlt,'mygui')
      greyzonehist_alt('close_Callback');
    end
    
    %Setup preferences
    g.DataLoaded = 0;
    g.Run = 0;

    %Current view and slice
    g.StartFrame = 1;

    %Make sure data is removed, store memory
    g.BalloonLevel = 0;
    g.BALLOON = [];
		g.EndoBalloonForce=[];
		g.EpiBalloonForce=[];
    g.cleardatalevelset;
    g.DATASETPREVIEW = [];
    g.EndoEDGE0 = [];
    g.EndoEDGE1 = [];
    g.EndoEDGE2 = [];
    g.EndoEDGE3 = [];
    g.EpiEDGE0 = [];
    g.EpiEDGE1 = [];
    g.EpiEDGE2 = [];
    g.EpiEDGE3 = [];
    g.EndoEdgeDetected = false;
    g.EpiEdgeDetected = false;
    g.MovieFrame = [];
    g.NumPoints = g.Pref.NumPoints; %Constant...
    NO = 1;
    SET = [];
    g.NeedToSave = false;
    set(g.Handles.filesaveicon,'enable','off');

    set(g.fig,'pointer','arrow');

    drawnow;
    
    g.startmodeplaceholders
    % JU: Set mode back to select (from SpecialGUI)
    % this only matters if not opening a .mat file
    if ~isempty(g)
      if g.DataLoaded
        updatetool('select');
      end
    end
    flushlog;
    end
    
    %-----------------------------
    function closeallnoguis(g,nos)
    %-----------------------------
    %Close all GUIs associated with one of image stacks nos
    allguis = fieldnames(g.GUI);
    for i = 1:numel(allguis)
      gui = g.GUI.(allguis{i});
      if isa(gui,'mygui')
        for no = nos
          try
            didclose = closeifhasno(gui,no);
          catch %#ok<CTCH>
            %If failed assume it is already closed.
            g.GUI.(allguis{i}) = [];
          end;
          if didclose
            g.GUI.(allguis{i}) = [];
            break
          end
        end
      end
    end
    end
    
    %---------------------------------
    function setthisframeonly(g,value)
    %---------------------------------
    %Callback to set this frame only mode.
    
    if nargin<2
      value=g.ThisFrameOnly;
    end
		
		g.ThisFrameOnly=value;
    
%    if ~strcmp(DATA.ProgramName,'Segment CMR')
%       set(g.Handles.thistimeframeonlycheckbox,'value',value);
%     end
    
    end
    
%   -----------------------------
%   function framemode_Callback(varargin)
%   -----------------------------
%   global DATA
%   g=varargin{1};
%   singlebutton=g.handles.singleframemodepushbutton;
%   multibutton=g.handles.multiframemodepushbutton;
%   if g.ThisFrameOnly
%     thisframeonly_Callback(g,~val,silent)
%   else
%     
%   end
%   end
    
    %-----------------------------------------------
    function thisframeonly_Callback(g,thisframeonly,silent) 
    %-----------------------------------------------
    %Sets segment in 'this frame only' mode. Changes made to segmentation,
    %copying etc are performed only for the current timeframe.
            
    if nargin<3
      silent=false;
    end
    
    %We should change setting
		if not(isequal(thisframeonly,g.ThisFrameOnly))
			g.ThisFrameOnly=thisframeonly;

      if g.ThisFrameOnly
        g.GUISettings.AxesColor=[1 1 1];
        g.GUISettings.SliceLineSpec=g.GUISettings.SliceLineSpecOneSlice;
      else
        g.GUISettings.AxesColor=g.GUISettings.AxesColorDefault;
        g.GUISettings.SliceLineSpec=g.GUISettings.SliceLineSpecMultiSlice;
      end

			g.setthisframeonly(g.ThisFrameOnly);
			
      %Call to update drawimage to get frame around panel with correct
      %color.
			if not(silent)
        set(g.Handles.imageaxes(g.CurrentPanel), 'xcolor',g.GUISettings.AxesColor,...
          'ycolor',g.GUISettings.AxesColor, 'linewidth',2.5, 'visible','on');
			end
		end

    end
    
    %----------------------------
    function cleardatalevelset(g)
    %----------------------------
    %Clear the struct g.LevelSet properly called when switching image
    %stacks.
    g.LevelSet.SpeedIM = [];
    g.LevelSet.EDGE0 = [];
    g.LevelSet.EDGE1 = [];
    g.LevelSet.GradientPart = [];
    end

    %---------------------
    function togglescar(g)
    %---------------------
    %Called by segment('updateviewicons'). Overloaded in CVQgui.
    global SET NO
    if (SET(NO).TSize >1)
      set(g.Handles.scarpushbutton,'enable','off');
      newscar = 'off';
      viability('viabilitymenu');
    elseif (SET(NO).TSize == 1)
      set(g.Handles.scarpushbutton,'enable','on');
      newscar = 'on';
      viability('viabilitymenu');
    end
    if strcmp(g.CurrentTheme,'scar') && strcmp(newscar,'off')
%      g.updateicons('lv');
    end
    end
    
    %----------------------------
    function fileopen_Callback(g) 
    %----------------------------
    %Loads a set of files as a 4D volume, and also updates preferences.
    %Calls the fcn openfile which displays the fileloader GUI.
    %Overloaded in CVQgui, Segment CMR gui, Segment SPECT gui.
    global SET

    segment('stopmovie_Callback');
    %set(g.Handles.fileopenicon,'state','off');

    if ~isempty(SET)
      m = mymenu('Image stacks already exist',...
        {'Replace existing image stacks',...
        'Add image stacks to existing'},g.GUI.Segment);

      switch m
        case 0
          return;
        case 1
          %close existing image stacks before loading
          silent = true;
          g.filecloseall_Callback(); %not silent
        case 2
          % do nothing, image stacks should be added to the existing image
          % stacks
      end;
    end
    %Call loading helper program.
    oldpointer=get(g.imagefig,'pointer');
    set(g.imagefig,'pointer','arrow');
    openfile;
    set(g.imagefig,'pointer',oldpointer);
    
    end
    
    %-------------------------------------------------
    function fail = filesaveallas_Callback(g,varargin)
    %-------------------------------------------------
		%'Save as' callback. Overloaded in GUI's that save to database as
		%default
    fail = filemenu('saveallas_helper',varargin{:});%save image stacks and segmentations
		end
    
    %------------------------
    function really = quit(g)
    %------------------------
    %Quit segment, also ask user if he/she is sure.
    %Overloaded in CVQgui

    really = true;
    if (g.DataLoaded && g.GUISettings.AskToExitProgram)
      if g.NeedToSave && ...
          ~yesno('Unsaved changes: are you sure you want to quit the program?')
        really = false;
        return
      end
      if ~yesno('Are you sure that you want to close all image stacks?',[],g.GUI.Segment);
        really = false;
        return;
      end
    end;

    sectra('stopsectra');
    externalpacs('stopexternalpacs');

    flushlog;
    end
    
    %------------------------------
    function renderstacksfrommat(g)
    %------------------------------
    %This function displays stacks from a mat files in main gui.
    %Typically called upon loading;
    global SET

    %Do not mess with .Silent here since it will be taken care of by lower
    %routine images.

    %Check if viewing preferences are stored in file.
    if ~isempty(g.ViewPanels)
      %  addtopanels(no); % this is messy when mat files are loaded into mat
      %  files, etc. Just do thing.
      %  NO=g.ViewPanels(g.CurrentPanel);
      %  no=NO;
      drawfunctions('drawall');
      segment('switchtopanel',g.CurrentPanel);
    elseif isempty(SET(1).View)
      g.ViewPanels = [];
      g.ViewPanelsType = {};
      g.ViewPanelsMatrix = {};
      g.ViewIM = {};
      g.Overlay = [];
      g.ViewMatrix = [1 1];
      g.CurrentPanel = 1;
      g.CurrentTool = 'select';
      g.CurrentTheme = 'lv';
      segment('addtopanels',1);
      drawfunctions('drawall');
      segment('switchtopanel',g.CurrentPanel);
    else
      g.ViewPanels = SET(1).View.ViewPanels;
      g.ViewPanelsType = SET(1).View.ViewPanelsType;
      g.ViewPanelsMatrix = SET(1).View.ViewPanelsMatrix;
      g.ViewIM = cell(1,length(g.ViewPanels));
      g.Overlay = struct('alphadata',[],'cdata',[]);
      g.Overlay(length(g.ViewPanels)) = struct('alphadata',[],'cdata',[]);
      if isfield(SET(1).View,'ViewMatrix')
        g.ViewMatrix = SET(1).View.ViewMatrix;
      else
        g.ViewMatrix = [];
      end;
      if isfield(SET(1).View,'CurrentPanel')
        g.CurrentPanel = SET(1).View.CurrentPanel;
        %     g.CurrentTheme = SET(1).View.CurrentTheme;
        %     g.CurrentTool = SET(1).View.CurrentTool;
      else
        g.CurrentPanel = 1;
        %     g.CurrentTheme = 'lv';
        %     g.CurrentTool = 'select';
      end;
      drawfunctions('drawall',g.ViewMatrix);
      if isfield(SET(1).View,'ThisFrameOnly')
        %Now default is always on.
        SET(1).View.ThisFrameOnly=1;
        g.ThisFrameOnly = SET(1).View.ThisFrameOnly;
%          if g.ThisFrameOnly
        g.GUISettings.AxesColor=[1 1 1];
%          else
%            g.GUISettings.AxesColor=g.GUISettings.AxesColorDefault;
%          end
        
        g.setthisframeonly(g.ThisFrameOnly);
      end;
      segment('switchtopanel',g.CurrentPanel);
      %   temp = g.CurrentTool;
      %   updateicons(g.CurrentTheme);
      %   updatetool(temp);
      if isfield(SET(1).View,'CurrentTheme')
%        g.updateicons(SET(1).View.CurrentTheme);
        updatetool(SET(1).View.CurrentTool);
      else
        %g.updateicons('lv');
        updatetool('select');
      end;
    end;
    
    %start with no colorbar
    %segment('viewhidecolorbar_Callback')
    %set(g.Handles.colorbar,'visible','off')

    drawfunctions('drawthumbnails',isempty(g.DATASETPREVIEW)); %False means no calcpreview
    segment('switchtoimagestack',g.ViewPanels(g.CurrentPanel),true); %force
    %segment('switchtoimagestack',no,true); %force
    mydisp('Image stacks loaded.');
    %endoffcalculation;
    
    %Enable all icons
    g.dataloadedplaceholders
    flushlog;
    end
    
    %-------------------------
    function updateicons(g,mode)
    %-------------------------
    %Updates icons when new mode is selected.
    %Overloaded in Segment GUI, CVQgui, RVQgui and Segment CMR GUI
    % If making major changes, make sure to
    % update csegment('updateicons')!!! /JU

    global SET NO

    if nargin<2
      mode = g.CurrentTheme;
    end;

    % if ismember(mode,{'scar'}) && (SET(NO).TSize>1)
    %   if ~g.Silent
    %     myfailed('Scar tools not enabled for time-resolved data. Perhaps delete time frames?',g.GUI.Segment)
    %   end
    %   return;
    % end

    %Reset Tool structure
    g.Tools = [];

    %Always
    g.Tools.selectslices = g.Handles.selectslicespushbutton;
    g.Tools.contrastbrightness = g.Handles.contrastbrightnesspushbutton;
    g.Tools.move = g.Handles.movepushbutton;
    g.Tools.scale = g.Handles.scalepushbutton;
    g.Tools.undo = g.Handles.undopushbutton;
    g.Tools.click3d = [];


    %LV
    g.Tools.endopen = [];
    g.Tools.epipen = [];
    g.Tools.autolv = [];
    g.Tools.autoendo = [];
    g.Tools.autoepi = [];
    g.Tools.refineendo = [];
    g.Tools.refineepi = [];
    g.Tools.removepapilar = [];
    g.Tools.putendopin = [];
    g.Tools.putepipin = [];
    g.Tools.clearlv = [];
    g.Tools.interpendo = [];
    g.Tools.interpepi = [];
    g.Tools.expandendo = [];
    g.Tools.contractendo = [];
    g.Tools.expandepi = [];
    g.Tools.contractepi = [];

    %RV
    g.Tools.rvendopen = [];
    g.Tools.rvepipen = [];
    g.Tools.autorvendo = [];
    g.Tools.autorvepi = [];
    g.Tools.refinervendo = [];
    g.Tools.putrvendopin = [];
    g.Tools.clearrv = [];
    g.Tools.interprvendo = [];
    g.Tools.interprvepi = [];

    %ROI
    g.Tools.roipen = [];
    g.Tools.putroi = [];
    g.Tools.refineroi = [];
    g.Tools.refineroinext = [];
    g.Tools.trackroi = [];
    g.Tools.roicolor = [];
    g.Tools.roitext = [];
    g.Tools.roiclear = [];

    %Scar
    g.Tools.scarpen = [];
    g.Tools.mopen = [];
    g.Tools.autoscar = [];
    g.Tools.rubberpen = [];
    g.Tools.rubber = [];
    g.Tools.scarundo = [];
    g.Tools.hideoverlay = [];

    %Misc
    g.Tools.point = [];
    g.Tools.measure = [];
    g.Tools.crop = [];
    g.Tools.flipx = [];
    g.Tools.flipy = [];
    g.Tools.rotate90 = [];

    %View
    g.Tools.viewcinescar = [];
    g.Tools.viewcinescarperf = [];
    g.Tools.viewstress = [];
    g.Tools.viewflow = [];

    %MaR
    g.Tools.marauto = [];
    g.Tools.marhide = [];
    g.Tools.marpen = [];
    g.Tools.marrubberpen = [];
    g.Tools.marrubber = [];

    %Reset all
    if ~strcmp(mode,g.CurrentTheme)
      set([...
        g.Handles.icon01 ...
        g.Handles.icon02 ...
        g.Handles.icon03 ...
        g.Handles.icon04 ...
        g.Handles.icon05 ...
        g.Handles.icon06 ...
        g.Handles.icon07 ...
        g.Handles.icon08 ...
        g.Handles.icon09...
        g.Handles.icon10 ...
        g.Handles.icon11 ...
        g.Handles.icon12 ...
        g.Handles.icon13 ...
        g.Handles.icon14...
        g.Handles.icon15],...
        'cdata',[],'Callback','','Tooltipstring','', ...
        'Interruptible','on','backgroundcolor',g.GUISettings.ButtonColor);
    end

    set([...
      g.Handles.leftventriclepushbutton ...
      g.Handles.rightventriclepushbutton ...
      g.Handles.t2starpushbutton ...
      g.Handles.perfusionpushbutton ...
      g.Handles.roipushbutton ...
      g.Handles.scarpushbutton ...
      g.Handles.marpushbutton ...
      g.Handles.reservedpushbutton ...
      g.Handles.annotationspushbutton ...
      g.Handles.tipushbutton ...
      g.Handles.viewpushbutton ...
      ],...
      'backgroundcolor',g.GUISettings.ButtonColor);

    oldtheme = g.CurrentTheme;
    switch mode
      case 'lv'
        g.CurrentTheme = 'lv';
        set(g.Handles.leftventriclepushbutton,'backgroundcolor',g.GUISettings.ButtonSelectedColor);
        
        g.Tools.autolv     = g.iconcdatahelper(g.Handles.icon01,g.Icons.autolv,'Automatic segmentation of left ventricle (Ctrl-L)');
        g.Tools.refineendo = g.iconcdatahelper(g.Handles.icon02,g.Icons.refineendo,'Refine endocardium (Ctrl-R)');
        g.Tools.refineepi  = g.iconcdatahelper(g.Handles.icon03,g.Icons.refineepi,'Refine epicardium (Ctrl-Shift-R)');
        g.Tools.putendopin = g.iconcdatahelper(g.Handles.icon04,g.Icons.putendopin,'Put endocardial pins that attract the contour');
        g.Tools.putepipin  = g.iconcdatahelper(g.Handles.icon05,g.Icons.putepipin,'Put epicardial pins that attract the contour');
        %g.Tools.clearlv    = g.iconcdatahelper(g.Handles.icon05,g.Icons.trash,'Clear LV segmentation in selected slices');
        
        g.Tools.interpendo = g.iconcdatahelper(g.Handles.icon06,g.Icons.interpendo,'Draw endocardium by point interpolation (Shift-B)');
        g.Tools.endopen    = g.iconcdatahelper(g.Handles.icon07,g.Icons.endopen,'Manually draw left ventricle endocardium (Shift-N)');
        g.Tools.expandendo = g.iconcdatahelper(g.Handles.icon08,g.Icons.expandendo,'Expand endocardium (Ctrl-E)');
        g.Tools.contractendo=g.iconcdatahelper(g.Handles.icon09,g.Icons.contractendo,'Contract endocardium (Ctrl-K)');
        g.Tools.propagateendo = g.iconcdatahelper(g.Handles.icon10,g.Icons.propagateendo,'Propagate endocardium forward and refine (Ctrl-F)');
        
        g.Tools.interpepi  = g.iconcdatahelper(g.Handles.icon11,g.Icons.interpepi,'Draw epicardium by point interpolation (Shift-H)');
        g.Tools.epipen     = g.iconcdatahelper(g.Handles.icon12,g.Icons.epipen,'Manually draw left ventricle epicardium (Shift-G)');
        g.Tools.expandepi  = g.iconcdatahelper(g.Handles.icon13,g.Icons.expandepi,'Expand epicardium (Ctrl-Alt-E)');
        g.Tools.contractepi= g.iconcdatahelper(g.Handles.icon14,g.Icons.contractepi,'Contract epicardium (Ctrl-Alt-K)');
        g.Tools.propagateepi = g.iconcdatahelper(g.Handles.icon15,g.Icons.propagateepi,'Propagate epicardium forward and refine (Ctrl-Shift-F)');
        
        %callbacks
        g.iconcallbackhelper(g.Tools.endopen,'drawendo');
        g.iconcallbackhelper(g.Tools.epipen,'drawepi');
        g.iconcallbackhelper(g.Tools.putendopin,'endopin');
        g.iconcallbackhelper(g.Tools.putepipin,'epipin');
        g.iconcallbackhelper(g.Tools.interpendo,'interpendo');
        g.iconcallbackhelper(g.Tools.interpepi,'interpepi');

        g.iconcallbackhelper(g.Tools.autolv,'autolv');
        set(g.Tools.refineendo, ...
          'Callback','lvpeter(''segmentrefineendo_Callback'',false,false)', ...
          'Interruptible','off','BusyAction','cancel');
        set(g.Tools.refineepi, ...
          'Callback','lvpeter(''segmentrefineepi_Callback'',false)', ...
          'Interruptible','off','BusyAction','cancel');
        set(g.Tools.expandendo, ...
          'Callback','lv(''segmentexpandcontract_Callback'',1,''endo'')',...
          'Interruptible','off','BusyAction','queue');
        set(g.Tools.contractendo, ...
          'Callback','lv(''segmentexpandcontract_Callback'',-1,''endo'')',...
          'Interruptible','off','BusyAction','queue');
        set(g.Tools.expandepi, ...
          'Callback','lv(''segmentexpandcontract_Callback'',1,''epi'')',...
          'Interruptible','off','BusyAction','queue');
        set(g.Tools.contractepi, ...
          'Callback','lv(''segmentexpandcontract_Callback'',-1,''epi'')',...
          'Interruptible','off','BusyAction','queue');
        set(g.Tools.propagateendo, ...
          'Callback','lvpeter(''segmentpropagateendo_Callback'')', ...
          'Interruptible','off','BusyAction','queue');
        set(g.Tools.propagateepi, ...
          'Callback','lvpeter(''segmentpropagateepi_Callback'')', ...
          'Interruptible','off','BusyAction','queue');
        %set(g.Tools.clearlv,'Callback','segmentation(''clearslices_Callback'',true,true,false,false)');

        if not(strcmp(g.CurrentTheme, oldtheme))
          updatetool('drawendo');
        end
      case 'rv'
        g.CurrentTheme = 'rv';
        set(g.Handles.rightventriclepushbutton,'backgroundcolor',g.GUISettings.ButtonSelectedColor);
        g.Tools.interprvendo = g.iconcdatahelper(g.Handles.icon01,g.Icons.interprvendo,'Draw RV endocardium by point interpolation');
        g.Tools.rvendopen = g.iconcdatahelper(g.Handles.icon02,g.Icons.rvendopen,'Manually draw right ventricle endocardium');
        g.Tools.refinervendo = g.iconcdatahelper(g.Handles.icon03,g.Icons.refinervendo,'Refine RV endocardium (Ctrl-Alt-R)');
        g.Tools.putrvendopin = g.iconcdatahelper(g.Handles.icon04,g.Icons.putrvendopin,'Put endocardial pins that attracts the contour');
        g.Tools.interprvepi  = g.iconcdatahelper(g.Handles.icon06,g.Icons.interprvepi,'Draw RV epicardium by point interpolation');
        g.Tools.rvepipen     = g.iconcdatahelper(g.Handles.icon07,g.Icons.rvepipen,'Manually draw right ventricle epicardium');
        g.Tools.autorvendo   = g.iconcdatahelper(g.Handles.icon12,g.Icons.autorvendo,'Auto segment RV endocardium (Ctrl-Alt-M)');
        g.Tools.clearrv      = g.iconcdatahelper(g.Handles.icon15,g.Icons.trash,'Clear RV segmentation in selected slices');

        g.iconcallbackhelper(g.Tools.rvendopen,'drawrvendo');
        g.iconcallbackhelper(g.Tools.rvepipen,'drawrvepi');
        g.iconcallbackhelper(g.Tools.putrvendopin,'rvendopin');
        g.iconcallbackhelper(g.Tools.interprvendo,'interprvendo');
        g.iconcallbackhelper(g.Tools.interprvepi,'interprvepi');

        g.iconcallbackhelper(g.Tools.autorvendo,'autorvendo');
        %set(g.Tools.autorvepi,'Callback','lv(''segmentepi_Callback'',false,''rv'')');%finns
        %ï¿½nnnu ej implementerad
        set(g.Tools.refinervendo, ...
          'Callback','rv(''segmentrefinervendo_Callback''),false,false',...
          'Interruptible','off','BusyAction','cancel');
        set(g.Tools.clearrv,'Callback','segmentation(''clearslices_Callback'',false,false,true,true)');

        if not(strcmp(g.CurrentTheme, oldtheme))
          updatetool('drawrvendo');
        end
      case 'roi'
        g.CurrentTheme = 'roi';
        set(g.Handles.roipushbutton,'backgroundcolor',g.GUISettings.ButtonSelectedColor);
        g.Tools.roipen = g.iconcdatahelper(g.Handles.icon01,g.Icons.roipen,'Manually draw/adjust ROI');
        g.Tools.putroi = g.iconcdatahelper(g.Handles.icon02,g.Icons.putroi,'Place a circular ROI');
        g.Tools.refineroi = g.iconcdatahelper(g.Handles.icon03,g.Icons.refineroi,'Refine ROI (Alt-R)');
        g.Tools.refineroinext = g.iconcdatahelper(g.Handles.icon04,g.Icons.refineroinext,'Propagate and refine ROI forward (Alt-F)');
        g.Tools.trackroi = g.iconcdatahelper(g.Handles.icon05,g.Icons.trackroi,'Track a vessel in all timeframes (Alt-T)');
        g.Tools.roicolor = g.iconcdatahelper(g.Handles.icon06,g.Icons.palette,'Color a ROI');
        g.Tools.roitext = g.iconcdatahelper(g.Handles.icon07,g.Icons.text,'Set ROI label');
        g.Tools.roiclear = g.iconcdatahelper(g.Handles.icon08,g.Icons.trash,'Remove ROI');

        g.iconcallbackhelper(g.Tools.roipen,'drawroi');
        g.iconcallbackhelper(g.Tools.putroi,'putroi');

        set(g.Tools.refineroi, ...
          'Callback','flow(''flowrefine_Callback'')',...
          'Interruptible','off','BusyAction','cancel');
        set(g.Tools.refineroinext, ...
          'Callback','flow(''flowpropagate_Callback'')',...
          'Interruptible','off','BusyAction','queue');
        g.iconcallbackhelper(g.Tools.trackroi,'trackroi');%set(g.Tools.trackroi,'Callback','flow(''flowtrackroi_Callback'')');
        set(g.Tools.roicolor,'Callback','roi(''roisetcolor_Callback'')');
        set(g.Tools.roitext,'Callback','roi(''roisetlabel_Callback'')');
        set(g.Tools.roiclear,'Callback','roi(''roidelete_Callback'')');

        if not(strcmp(g.CurrentTheme, oldtheme))
          updatetool('drawroi');
        end
      case 'scar'
        g.CurrentTheme = 'scar';
        set(g.Handles.scarpushbutton,'backgroundcolor',g.GUISettings.ButtonSelectedColor);
        g.Tools.autoscar = g.iconcdatahelper(g.Handles.icon01,g.Icons.autoscar,'Automatically delineate scar regions');
        g.Tools.scarpen = g.iconcdatahelper(g.Handles.icon06,g.Icons.scarpen,'Manually add scar regions');
        g.Tools.mopen = g.iconcdatahelper(g.Handles.icon07,g.Icons.mopen,'Manually add region of microvascular obstruction');
        g.Tools.rubberpen = g.iconcdatahelper(g.Handles.icon08,g.Icons.rubberpen,'Manually remove scar regions');
        g.Tools.rubber = g.iconcdatahelper(g.Handles.icon09,g.Icons.rubber,'Draw to remove manual interaction');
        %    g.Tools.scarundo = g.iconcdatahelper(g.Handles.icon10,g.Icons.undo,'Undo manual interaction');

        g.iconcallbackhelper(g.Tools.autoscar,'autoscar');%set(g.Tools.autoscar,'Callback','viability(''viabilitycalc'');segment(''drawimageno'');');
        g.iconcallbackhelper(g.Tools.scarpen,'drawscar');
        g.iconcallbackhelper(g.Tools.mopen,'drawmo');
        g.iconcallbackhelper(g.Tools.rubberpen,'drawrubberpen');
        g.iconcallbackhelper(g.Tools.rubber,'drawrubber');

        if not(strcmp(g.CurrentTheme, oldtheme))
          if isempty(SET(NO).Scar)
            updatetool('select');
          else
            updatetool('drawscar');
          end;
        end;
      case 'misc'
        g.CurrentTheme = 'misc';
        set(g.Handles.annotationspushbutton,'backgroundcolor',g.GUISettings.ButtonSelectedColor);
        g.Tools.point = g.iconcdatahelper(g.Handles.icon01,g.Icons.point,'Place annotation point');
        g.Tools.measure = g.iconcdatahelper(g.Handles.icon02,g.Icons.measure,'Place measurement caliper');
        g.Tools.crop = g.iconcdatahelper(g.Handles.icon03,g.Icons.crop,'Crop image stack');
        g.Tools.autocropall = g.iconcdatahelper(g.Handles.icon04,g.Icons.autocropall,'Automatically crop all image stacks');        
        g.Tools.flipx = g.iconcdatahelper(g.Handles.icon06,g.Icons.flipx,'Flip image stack in x and z direction');
        g.Tools.flipy = g.iconcdatahelper(g.Handles.icon07,g.Icons.flipy,'Flip image stack in y and z direction');
        %g.Tools.rotate90 = g.iconcdatahelper(g.Handles.icon08,g.Icons.rotate90,'Rotate image stack 90 degrees');
        g.Tools.click3d = g.iconcdatahelper(g.Handles.icon11,g.Icons.click3d,'Click image to show point in all views');
        g.iconcallbackhelper(g.Tools.point,'point');
        g.iconcallbackhelper(g.Tools.measure,'measure');
        g.iconcallbackhelper(g.Tools.crop,'crop');
        g.iconcallbackhelper(g.Tools.autocropall,'autocropall');
        g.iconcallbackhelper(g.Tools.click3d,'click3d');
        set(g.Tools.flipx,'Callback','tools(''flipx_Callback'')');
        set(g.Tools.flipy,'Callback','tools(''flipy_Callback'')');
        %set(g.Tools.rotate90,'Callback','segment(''toolsrotate90right_Callback'')');

        if not(strcmp(g.CurrentTheme, oldtheme))
          updatetool('point');
        end
      case 'view'
        g.CurrentTheme = 'view';
        set(g.Handles.viewpushbutton,'backgroundcolor',g.GUISettings.ButtonSelectedColor);
        g.Tools.viewcine = g.iconcdatahelper(g.Handles.icon01,g.Icons.panel1,'View Cine images');
        g.Tools.viewcinescar = g.iconcdatahelper(g.Handles.icon02,g.Icons.row2,'View Cine and Scar');
        g.Tools.viewcinescarperf = g.iconcdatahelper(g.Handles.icon03,g.Icons.row3,'View Cine, Scar and Perfusion');
        g.Tools.viewstress = g.iconcdatahelper(g.Handles.icon04,g.Icons.row4,'View Stress');
        g.Tools.viewflow = g.iconcdatahelper(g.Handles.icon05,g.Icons.panel2,'View flow');
        set(g.Tools.viewcine,'Callback','segment(''viewspecial_Callback'',''cine'')');
        set(g.Tools.viewcinescar,'Callback','segment(''viewspecial_Callback'',''cinescar'')');
        set(g.Tools.viewcinescarperf,'Callback','segment(''viewspecial_Callback'',''cinescarperf'')');
        set(g.Tools.viewstress,'Callback','segment(''viewspecial_Callback'',''stress'')');
        set(g.Tools.viewflow,'Callback','segment(''viewspecial_Callback'',''flow'')');
        if not(strcmp(g.CurrentTheme, oldtheme))
          updatetool('select');
        end
      case 'mar'
        g.CurrentTheme = 'mar';
        set(g.Handles.marpushbutton,'backgroundcolor',g.GUISettings.ButtonSelectedColor);
        g.Tools.marauto = g.iconcdatahelper(g.Handles.icon01,g.Icons.automar,'Automatically delineate MaR regions');
        g.Tools.marpen = g.iconcdatahelper(g.Handles.icon03,g.Icons.marpen,'Manually add MaR regions');
        g.Tools.marrubberpen = g.iconcdatahelper(g.Handles.icon04,g.Icons.rubberpen,'Manually remove MaR region');
        g.Tools.marrubber = g.iconcdatahelper(g.Handles.icon05,g.Icons.rubber,'Draw to remove manual interaction');

        g.iconcallbackhelper(g.Tools.marauto,'automar');%set(g.Tools.marauto,'Callback','mar(''auto_Callback'');segment(''drawimageno'');');
        g.iconcallbackhelper(g.Tools.marpen,'drawmarpen');
        g.iconcallbackhelper(g.Tools.marrubberpen,'drawmarrubberpen');
        g.iconcallbackhelper(g.Tools.marrubber,'drawmarrubber');

        if not(strcmp(g.CurrentTheme, oldtheme))
          updatetool('drawmarpen');
        end
    end;
    if (isequal(oldtheme,'scar')&&not(isequal(g.CurrentTheme,'scar')))||...
        (isequal(oldtheme,'mar')&&not(isequal(g.CurrentTheme,'mar')))||...
        isequal(g.CurrentTheme,'scar')||isequal(g.CurrentTheme,'mar')
      drawfunctions('drawimageno');
    end
    end
    
    %---------------------------------
    function updatetimethings(g,no,mode) %#ok<INUSD>
    %---------------------------------
    %This fcn is called from switchimagestack and disables/enables
    %features that are not available depending of the image stack is
    %timeresolved or not.
    %Overloaded in CVQgui.
    global SET NO

%     if nargin==1
%       no = NO;
%     end;

%     %Check if timeresolved
%     if SET(no).TSize>1
%       set([g.Handles.mmodeviewicon],'enable','on');
%       set(g.Handles.volumeaxes,'visible','on');
%       set(g.Handles.flowaxes,'visible','on');
%       set(g.Handles.timebaraxes,'visible','on');
%       set([...
%         g.Handles.playmovieicon ...
%         g.Handles.playallicon ...
%         g.Handles.nextframeicon ...
%         g.Handles.nextallframeicon ...
%         g.Handles.previousframeicon ...
%         g.Handles.previousallframeicon ...
%         g.Handles.fasterframerateicon ...
%         g.Handles.slowerframerateicon ...
%         g.Handles.cinetoolicon ...
%         g.Handles.playmoviemenu ...
%         g.Handles.stopmoviemenu ...
%         g.Handles.reportslicemenu ...
%         g.Handles.reportradvelmenu ...
%         g.Handles.reportvolumecurvemenu],'enable','on');
%     else
%       %No time, disable things
%       set(g.Handles.volumeaxes,'Visible','off');
%       set(g.Handles.timebaraxes,'Visible','off');
%       set([...
%         g.Handles.mmodeviewicon ...
%         g.Handles.playmovieicon ...
%         g.Handles.playallicon ...
%         g.Handles.nextframeicon ...
%         g.Handles.nextallframeicon ...
%         g.Handles.previousframeicon ...
%         g.Handles.previousallframeicon ...
%         g.Handles.fasterframerateicon ...
%         g.Handles.slowerframerateicon ...
%         g.Handles.cinetoolicon ...
%         g.Handles.playmoviemenu ...
%         g.Handles.stopmoviemenu ...
%         g.Handles.reportslicemenu ...
%         g.Handles.reportradvelmenu ...
%         g.Handles.reportvolumecurvemenu],'enable','off');
%     end;
    
    end

    %---------------------------------------
    function switchtoimagestack(g, no,force)
    %---------------------------------------
    %This function makes no current image stack, and updates graphics.
    global NO

    if g.Silent %&& not(force) %PH:
      NO = no;
      return;
    end;

    %segment('updateviewicons',no);

    if isequal(NO,no)&&(nargin==2)
      return;
    end;

    if isequal(NO,no)&&(nargin==3)
      if ~force
        return;
      end;
    end;

    %--- Check if there are open report windows.

    %Check deleted windows
    ind = true(1,length(g.BlockingFigs));
    for loop=1:length(ind)
      try
        get(g.BlockingFigs);
      catch %#ok<CTCH>
        ind(loop) = false;
      end;
    end;

    %Remove deleted windows
    g.BlockingFigs = g.BlockingFigs(ind);
    g.switchtoimagestack_part2(no);
    
    %Play again if play button is indented
    stateandicon= segment('iconson','play');
    if stateandicon{1}
      segment('playall_Callback');
    end
    end
    
   
    %--------------------------------------
    function drawroiinpanel(g,panel,docalc)
    %--------------------------------------
    %Draw ROI's in one slice mode.
    %Not overloaded anywhere.
    global SET

    if nargin < 3
      docalc = true;
    end
    %First try to remove old ones
    try
      delete(g.Handles.roicontour{panel});
      delete(g.Handles.roitext{panel});
    catch %#ok<CTCH>
      %do nothing if failed.
    end

    %If flow then take ROI from flow.
    no = g.ViewPanels(panel);
    if ~isempty(SET(no).Parent)
      no = SET(no).Parent;
    end;

    % %Check on how many to plot
    % ind=(SET(no).Roi(:).Z==SET(no).CurrentSlice);
    % rois2plot = find(ind);
    % numrois2plot = length(rois2plot);
    % numrois=length(SET(no).RoiZ);


    %Draw with no forces phase images to use from mag images
    %if isequal(get(g.Handles.hideroiicon,'state'),'off')
      %Plot them
      g.Handles.roicontour{panel}= nan(1,SET(no).RoiN);
      g.Handles.roitext{panel}= nan(1,SET(no).RoiN);
      hold(g.Handles.imageaxes(panel),'on');
      for loop=1:SET(no).RoiN
        if docalc
          [~,area]=calcfunctions('calcroiarea',no,loop);
          SET(no).Roi(loop).Area = area;
          [m,sd]=calcfunctions('calcroiintensity',no,loop);
          SET(no).Roi(loop).Mean = m;
          SET(no).Roi(loop).StD = sd;
        end
        if SET(no).Roi(loop).Sign > 0
          roisign = '';
        else
          roisign = ' (-)';
        end
        %Input to force black background around corner text to improve visibility
        if g.Pref.BackgroundColor
          bgcolor = 'blue';
         else
          bgcolor = 'none';
        end
        g.Handles.roicontour{panel}(loop) = plot(g.Handles.imageaxes(panel),...
          SET(no).Roi(loop).Y(:,SET(no).CurrentTimeFrame),...
          SET(no).Roi(loop).X(:,SET(no).CurrentTimeFrame),...
          SET(no).Roi(loop).LineSpec);
        [ymin,ix] = min(SET(no).Roi(loop).Y(:,SET(no).CurrentTimeFrame));
        g.Handles.roitext{panel}(loop) = text(...
          ymin-1,SET(no).Roi(loop).X(ix,SET(no).CurrentTimeFrame),...
          {sprintf('%s%s',SET(no).Roi(loop).Name,roisign), ...
          sprintf('%3.1f [cm^2]', ...
          SET(no).Roi(loop).Area(SET(no).CurrentTimeFrame)), ...
          sprintf('%3.1f ± %3.1f', ...
          SET(no).Roi(loop).Mean(SET(no).CurrentTimeFrame), ...
          SET(no).Roi(loop).StD(SET(no).CurrentTimeFrame))},...
          'HorizontalAlignment','right','VerticalAlignment','middle', ...
          'parent',g.Handles.imageaxes(panel),'backgroundcolor',bgcolor);
        set(g.Handles.roitext{panel}(loop),'color',[1 1 1]);
        if ~(ismember(SET(no).CurrentSlice,SET(no).Roi(loop).Z) && ...
            ismember(SET(no).CurrentTimeFrame,SET(no).Roi(loop).T))
          set(g.Handles.roitext{panel}(loop),'Position',[nan nan]);
          set(g.Handles.roicontour{panel}(loop),'XData',nan,'YData',nan);
        end
      end;

      %Width of ROI's
      if g.Pref.LineWidth>0
        set(g.Handles.roicontour{panel},'linewidth',g.Pref.LineWidth);
        set(g.Handles.roicontour{panel}(SET(no).RoiCurrent), ...
          'linewidth',g.Pref.LineWidth+1);
      else
        set(g.Handles.roicontour{panel},'visible','off');
      end;

      if isequal(get(g.Handles.hidetexticon,'state'),'on')
        set(g.Handles.roitext{panel},'visible','off');
      end;

      hold(g.Handles.imageaxes(panel),'off');
    %else
%       %Set up handles do not plot
%       g.Handles.roicontour{panel} = zeros(1,SET(no).RoiN);
%       hold(g.Handles.imageaxes(panel),'on');
%       g.Handles.roicontour{panel}(:) = plot(g.Handles.imageaxes(panel),NaN,NaN);
%       hold(g.Handles.imageaxes(panel),'off');
%     end;
%     
    if ~isempty(g.FlowNO) && (no == g.FlowNO || ismember(g.FlowNO,SET(no).Linked)) && ~isempty(g.FlowROI)
      if SET(no).RoiN <1
        g.FlowROI = [];
      else
        calcfunctions('calcflow',no);
      end
    end
    g.updateaxestables('flow',no);

    end
    
    %-----------------------------------
    function plotrois(g,panel,no,docalc)
    %-----------------------------------
    %Plot roi's if existing. Overloaded in CVQgui.
    global DATA SET
    
    if nargin < 4
      docalc = true;
    end
    
    %First try to remove old ones
    try
      delete(g.Handles.roicontour{panel});
      delete(g.Handles.roitext{panel});
    catch %#ok<CTCH>
      %do nothing if failed.
    end
    
    g.Handles.roicontour{panel} = nan(1,SET(no).RoiN);
    g.Handles.roitext{panel} = [];
    hold(g.Handles.imageaxes(panel),'on');
    for loop=1:SET(no).RoiN
      xofs = SET(no).XSize*floor((round(SET(no).Roi(loop).Z)-1)/DATA.ViewPanelsMatrix{panel}(2));
      yofs = SET(no).YSize*mod(round(SET(no).Roi(loop).Z)-1,DATA.ViewPanelsMatrix{panel}(2));
      g.Handles.roicontour{panel}(loop) = plot(g.Handles.imageaxes(panel),...
        SET(no).Roi(loop).Y(:,SET(no).CurrentTimeFrame)+yofs,...
        SET(no).Roi(loop).X(:,SET(no).CurrentTimeFrame)+xofs,...
        SET(no).Roi(loop).LineSpec);
      if docalc
        [~,SET(no).Roi(loop).Area]=calcfunctions('calcroiarea',no,loop);
        [m,sd]=calcfunctions('calcroiintensity',no,loop);
        SET(no).Roi(loop).Mean = m;
        SET(no).Roi(loop).StD = sd;
      end
      if g.Pref.LineWidth>0
        if ismember(loop,SET(no).RoiCurrent)
          set(g.Handles.roicontour{panel}(loop),'linewidth',g.Pref.LineWidth+1);
        else
          set(g.Handles.roicontour{panel}(loop),'linewidth',g.Pref.LineWidth);
        end
      else
        set(g.Handles.roicontour{panel}(loop),'visible','off');
      end;
    end;
    hold(g.Handles.imageaxes(panel),'off');
    end
    
    %-------------------------------------------
    function segmentclearalllv_Callback(g,force) %#ok<INUSD,DEFNU>
    %-------------------------------------------
    %Clear all LV segmentation, both endo and epi
    %Overloaded in CVQgui

    if nargin==1
      if ~yesno('Do you really want to remove all LV segmentation ?',[],g.GUI.Segment);
          myfailed('Aborted by user.',g.GUI.Segment);
        return;
      end;
    end;
    
    segmentation('clearalllv_Callback');
    end
    
    %-------------------------------------------
    function segmentclearallrv_Callback(g,force) %#ok<INUSD,DEFNU>
    %-------------------------------------------
    %Clear all RV segmentation, both endo and epi
    %Overloaded in CVQgui

    if nargin==1
      if ~yesno('Do you really want to remove all RV segmentation ?',[],g.GUI.Segment);
          myfailed('Aborted by user.',g.GUI.Segment);
        return;
      end;
    end;
    
    segmentation('clearallrv_Callback');
    end
    
    %-----------------------------------------
    function segmentclearall_Callback(g,force) %#ok<INUSD>
    %-----------------------------------------
    %Clear all segmentation, both endo and epi, lv and rv mar and scar.
    %Overloaded in CVQgui and RVQgui.
    global SET NO
    
    if nargin==1
      msg='This removes all existing segmentation lv, rv, roi mar and scar. Are you sure?';
%       if isempty(SET(NO).Scar)
%         msg='Do you really want to remove all segmentation (both LV and RV) ?';
%       else
%         msg ='Do you really want to remove all segmentation (both LV, RV and Scar) ?';
%       end
      if ~yesno(msg,[],g.GUI.Segment);
        myfailed('Aborted by user.',g.GUI.Segment);
        return;
      end;
    end;
    
    viability('viabilityclear_Callback'); 
    roi('roiclearall_Callback')
    segmentation('clearall_Callback'); 
    mar('clearall_Callback');
    end
    
    %-----------------------------------------------------
    function segmentclearallbutsystolediastole_Callback(g)
    %-----------------------------------------------------
    %Clears all segmentation in all timeframes but systole and diastole.
    global SET NO

    if SET(NO).TSize<2
      myfailed('Not timeresolved data, aborting.',g.GUI.Segment);
      return;
    end;

    if isequal(SET(NO).EDT,SET(NO).EST)
      myfailed('Systole and diastole occurs at the same time frame, aborting.',g.GUI.Segment);
      return;
    end;

    segmentation('removeallpins_Callback',true); %side effect calls enableundo
    %Create index structure
    ind = true(1,SET(NO).TSize);
    ind(SET(NO).EDT) = false;
    ind(SET(NO).EST) = false;
    arg = struct('endo',true,'epi',true,'rvendo',true,'rvepi',true);
    indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind);
    segmentation('removeallinterp_Callback',true,[],arg,indarg);

    if ~isempty(SET(NO).EndoX)
      SET(NO).EndoX(:,ind,:) = NaN;
      SET(NO).EndoY(:,ind,:) = NaN;
    end;

    if ~isempty(SET(NO).EpiX)
      SET(NO).EpiX(:,ind,:) = NaN;
      SET(NO).EpiY(:,ind,:) = NaN;
    end;

    if ~isempty(SET(NO).RVEndoX)
      SET(NO).RVEndoX(:,ind,:) = NaN;
      SET(NO).RVEndoY(:,ind,:) = NaN;
    end;

    if ~isempty(SET(NO).RVEpiX)
      SET(NO).RVEpiX(:,ind,:) = NaN;
      SET(NO).RVEpiY(:,ind,:) = NaN;
    end;

    SET(NO).EndoDraged(ind,:) = false;
    SET(NO).EpiDraged(ind,:) = false;

    segment('updatevolume');
    segment('updatemodeldisplay');
    drawfunctions('drawsliceno');
    end

    %-------------------------------------------------------
    function segmentclearalllvbutsystolediastole_Callback(g)
    %-------------------------------------------------------
    %Clears all LV segmentation except in systole and diastole.
    %Overloaded in CVQgui
    global SET NO

    if SET(NO).TSize<2
      myfailed('Not timeresolved data, aborting.',g.GUI.Segment);
      return;
    end;

    if isequal(SET(NO).EDT,SET(NO).EST)
      myfailed('Systole and diastole occurs at the same time frame, aborting.',g.GUI.Segment);
      return;
    end;

    segmentation('removeallpins_Callback',true,1,1,0,0); %side effect calls enableundo
    %Create index structure
    ind = true(1,SET(NO).TSize);
    ind(SET(NO).EDT) = false;
    ind(SET(NO).EST) = false;
    arg = struct('endo',true,'epi',true,'rvendo',false,'rvepi',false);
    indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind);
    segmentation('removeallinterp_Callback',true,[],arg,indarg);

    %Create index structure
    ind = true(1,SET(NO).TSize);
    ind(SET(NO).EDT) = false;
    ind(SET(NO).EST) = false;

    if ~isempty(SET(NO).EndoX)
      SET(NO).EndoX(:,ind,:) = NaN;
      SET(NO).EndoY(:,ind,:) = NaN;
    end;

    if ~isempty(SET(NO).EpiX)
      SET(NO).EpiX(:,ind,:) = NaN;
      SET(NO).EpiY(:,ind,:) = NaN;
    end;

    SET(NO).EndoDraged(ind,:) = false;
    SET(NO).EpiDraged(ind,:) = false;

    lvsegchanged = true; segment('updatevolume',lvsegchanged);
    segment('updatemodeldisplay');
    drawfunctions('drawsliceno');

    end

    %-------------------------------------------------------
    function segmentclearallrvbutsystolediastole_Callback(g) 
    %-------------------------------------------------------
    %Clears all RV segmentation except in systole and diastole.
    %Overloaded in CVQgui
    global SET NO
    
    if SET(NO).TSize<2
      myfailed('Not timeresolved data, aborting.',g.GUI.Segment);
      return;
    end;

    if isequal(SET(NO).EDT,SET(NO).EST)
      myfailed('Systole and diastole occurs at the same time frame, aborting.',g.GUI.Segment);
      return;
    end;
    
    segmentation('removeallpins_Callback',true,0,0,1,1); %side effect calls enableundo
    %Create index structure
    ind = true(1,SET(NO).TSize);
    ind(SET(NO).EDT) = false;
    ind(SET(NO).EST) = false;
    arg = struct('endo',false,'epi',false,'rvendo',true,'rvepi',true);
    indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind);
    segmentation('removeallinterp_Callback',true,[],arg,indarg);

    if ~isempty(SET(NO).RVEndoX)
      SET(NO).RVEndoX(:,ind,:) = NaN;
      SET(NO).RVEndoY(:,ind,:) = NaN;
    end;

    if ~isempty(SET(NO).RVEpiX)
      SET(NO).RVEpiX(:,ind,:) = NaN;
      SET(NO).RVEpiY(:,ind,:) = NaN;
    end;

    SET(NO).EndoDraged(ind,:) = false;
    SET(NO).EpiDraged(ind,:) = false;

    segment('updatevolume');
    segment('updatemodeldisplay');
    drawfunctions('drawsliceno');

    end
    
    %---------------------------------
    function point_Buttondown(g,panel) 
    %---------------------------------
    %Button down function when point tool is active.
    %Overloaded in CVQgui
    global SET NO DATA

    if g.Interactionlock
      return;
    end;

    segment('switchtopanel',panel);
    no = NO;
    
    switch get(g.imagefig,'SelectionType')

      case 'alt'
%         g.contextmenu;
        [y,x,slice] = segment('getclickedcoords');
        ind = findfunctions('findclosestannotationpoint',x,y,slice);

        SET(NO).Point.X(ind)=x;
        SET(NO).Point.Y(ind)=y;
        SET(NO).Point.Z(ind)=slice;
        
        drawfunctions('drawimageno');
      case 'normal'

        %Use to point to mag data set
        if ~isempty(SET(NO).Parent)
          no = SET(NO).Parent;
        end;

        [x,y,slice] = segment('getclickedcoords');
        % If slice has changed, make sure montage/one are in sync
        if (slice>SET(no).ZSize)
          return;
        end
        segment('switchtoslice',slice);

        placetimeresolved = isequal(get(DATA.Handles.placetimeresolvedpoints,'checked'),'on');
        
        if isempty(SET(no).Point)
          annotationpoint('pointclearall'); % calls enableundo
        else
          tools('enableundo',no);
        end;

        set(g.Handles.hidepointsicon,'state','off');

        if (~placetimeresolved) || isempty(SET(no).Point) || (isempty(SET(no).Point.Label))
          %Ask for name if not placetimeresolved or first point
          menuitems = {...
            'Apex',...
            'RV insertion',...
            'AV plane',...
            'TV plane',...
            'P1',...
            'P2',...
            'General',...
            'Sector start',...
            'User defined ...'};

          c = mymenu('Name point',menuitems,g.GUI.Segment);
          if c<1
            return;
          end;
          if c <=( length(menuitems)-1)
            s = menuitems{c};
          else

            %last is User defined.
            if c == (length(menuitems))
              s = inputdlg({'Enter name'},'Name',1,{sprintf('Point_%d',length(SET(no).Point.X))});
              if isempty(s)
                myfailed('Invalid name.',g.GUI.Segment);
                return;
              end;
            end;
          end;
        else
          s = SET(no).Point.Label{end}; %Instead take name from last point
        end;
        
        if isfield(SET(no).StrainTagging,'LVupdated') && strcmp(s,'AV plane')
          SET(no).StrainTagging.LVupdated=1;
        end
        tvplanebool=strcmp(SET(no).Point.Label,'TV plane');
        if strcmp(s,'TV plane') && sum(tvplanebool)==2
          ind=find(tvplanebool,1);
          SET(no).Point.X(ind)=[];
          SET(no).Point.Y(ind)=[];
          SET(no).Point.T(ind)=[];
          SET(no).Point.Z(ind)=[];
          SET(no).Point.Label(ind)=[];
        end
        
        SET(no).Point.X = [SET(no).Point.X y];
        SET(no).Point.Y = [SET(no).Point.Y x];
        SET(no).Point.T = [SET(no).Point.T NaN]; %Default non time resolved.
        SET(no).Point.Z = [SET(no).Point.Z SET(no).CurrentSlice];
        SET(no).Point.Label = [SET(no).Point.Label s];
        
        drawfunctions('drawimageno');
        
        %If placetimeresolved option is enable then show this in this frame
        %only.
        if placetimeresolved
          g.pointshowthisframeonly_Callback;
        end;
        
    end
    
    %If straintagging initiated adjust LVupdated
    if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
      SET(no).StrainTagging.LVupdated = 1;
    end
    
    datacursormode off;    
    end; %end of point_Buttondown;
    
    %-----------------------------------------
    function measurerenamethis_Callback(g,arg) 
    %-----------------------------------------
    %Rename current measurement.
    %Overloaded in CVQgui.
    global SET NO
        
    ask = false;
    if nargin > 1 && strcmp(arg,'ask')
      ask = true;
    end

    %Use to point to mag data set
    no = NO;
    if ~isempty(SET(NO).Parent)
      no = SET(NO).Parent;
    end;

    if isempty(SET(no).Measure)
      myfailed('No measurement exist.');
      return;
    end
    
    if ask && length(SET(no).Measure) > 1
      m = mymenu('Select measurement',...
        {SET(no).Measure.LongName},g.GUI.Segment);
      if m > 0
        g.MeasureN = m;
      else
        return
      end
    end
    
    n = g.MeasureN;
    if (length(SET(no).Measure)<n)||(n<1)
      myfailed('Could not rename current measurement.',g.GUI.Segment);
      return;
    end;
    tools('enableundo',no);
    
    %Get new name
    [stri,lstr] = g.measureasklabel;
    if ~isempty(stri)
      SET(no).Measure(n).Name = stri;
      SET(no).Measure(n).LongName = lstr;
    else
      myfailed('Invalid name.',g.GUI.Segment);
      return;     
    end;

    drawfunctions('drawimageno');
    end
    
    %----------------------------------------
    function measureclearthis_Callback(g,arg) 
    %----------------------------------------
    %Clear current measurement.
    %Overloaded in RVQgui
    global SET NO
    
    ask = false;
    if nargin > 1 && strcmp(arg,'ask')
      ask = true;
    end
    
    %Use to point to mag data set
    no = NO;
    if ~isempty(SET(NO).Parent)
      no = SET(NO).Parent;
    end;
    
    if isempty(SET(no).Measure)
      myfailed('No measurement exist.');
      return;
    end
        
    if ask && length(SET(no).Measure) > 1
      m = mymenu('Select measurement',...
        {SET(no).Measure.LongName},g.GUI.Segment);
      if m > 0
        g.MeasureN = m;
      else
        return
      end
    end

    n = g.MeasureN;
    if (length(SET(no).Measure)<n)||(n<1)
      myfailed('Could not delete current measurement.',g.GUI.Segment);
      return;
    end;
    tools('enableundo',no);
    
    g.updateaxestables('measure',true,n);

    ind = true(1,length(SET(no).Measure));
    ind(n) = false;
    SET(no).Measure = SET(no).Measure(logical(ind));
    drawfunctions('drawimageno');
    end
    
    %-----------------------------------
    function measureclearall_Callback(g) 
    %-----------------------------------
    %Clear all measurements.
    global SET NO

    %Use to point to mag data set
    no = NO;
    if ~isempty(SET(NO).Parent)
      no = SET(NO).Parent;
    end;
    tools('enableundo',no);

    set(g.Handles.hidemeasuresicon,'state','off');
    SET(no).Measure = [];
    drawfunctions('drawimageno');

    g.updateaxestables('measureclearall');
    end
    
    %------------------------------
    function [stri,lstr] = measureasklabel(g)
    %------------------------------
    %Asks for a label of a measurement.
    %Overloaded in CVQgui and RVQgui
    
    m = mymenu('Select measurement type',...
      {'ASW (Anterior Septal Wall Thickness)',...
      'PLW (Posterior Lateral Wall Thickness)',...
      'ESD (End Systolic Dimension)',...
      'EDD (End Diastolic Dimension)',...
      'ESL (End Systolic Length)',...
      'EDL (End Diastolic Length)',...
      'AA (Ascending Aorta)',...
      'DA (Descending Aorta)',...
      'ARD (Aortic Root Diameter)',...
      'AL (Aortic Length)',...
      'RVmaj (RV Major Axis)',...
      'RVmin (RV Minor Axis)',...
      'User defined ...'},g.GUI.Segment);

    switch m
      case 0
        stri = '';
        lstr = '';
        return;
      case 1
        stri = 'ASW';
        lstr = 'Anterior Septal Wall Thickness';
      case 2
        stri = 'PLW';
        lstr = 'Posterior Lateral Wall Thickness';
      case 3
        stri = 'ESD';
        lstr = 'End Systolic Dimension';
      case 4
        stri = 'EDD';
        lstr = 'End Diastolic Dimension';
      case 5
        stri = 'ESL';
        lstr = 'End Systolic Length';
      case 6
        stri = 'EDL';
        lstr = 'End Diastolic Length';
      case 7
        stri = 'AA';
        lstr = 'Ascending Aorta';
      case 8
        stri = 'DA';
        lstr = 'Descending Aorta';
      case 9
        stri = 'ARD';
        lstr = 'Aortic Root Diameter';
      case 10
        stri = 'AL';
        lstr = 'Aortic Length';
      case 11
        stri = 'RVmaj';
        lstr = 'RV Major Axis';
      case 12
        stri = 'RVmin';
        lstr = 'RV Minor Axis';
      case 13
        s = inputdlg({'Enter name'},'Name',1,{sprintf('Measure_%d',g.MeasureN)});
        if isempty(s)
          stri = '';
          lstr = '';
        else
          stri = s{1};
          lstr = s{1};
        end;
    end;
    end
    
    %---------------------------
    function measure_Buttonup(g) 
    %---------------------------
    %Button up function for measurements.
    %Overloaded in CVQgui
    global SET NO

    %Use to point to mag data set
    no = NO;
    if ~isempty(SET(NO).Parent)
      no = SET(NO).Parent;
    end;
    tools('enableundo',no);

    %Restore so no motion is called
    set(g.imagefig,'WindowButtonMotionFcn','');

    %Restore main buttonup function
    set(g.imagefig,'WindowButtonUpFcn',...
      sprintf('%s(''buttonup_Callback'')','segment'));

    %Store template, ORDER is important!
    switch g.ViewPanelsType{g.CurrentPanel}
      case 'hla'
        M.X = g.MeasureZ; %SET(no).HLA.slice * [1; 1];
        M.Y = g.MeasureY;
        M.Z = g.MeasureX;
      case 'vla'
        M.X = g.MeasureY;
        M.Y = g.MeasureZ; %SET(no).VLA.slice * [1; 1];
        M.Z = g.MeasureX;
      case 'gla'
        [M.X,M.Y,M.Z] = calcfunctions('gla2sax',g.MeasureX,g.MeasureY,no);
      otherwise
        M.X = g.MeasureX;
        M.Y = g.MeasureY;
        M.Z = g.MeasureZ; %SET(no).CurrentSlice * [1; 1];
    end
    
    M.Length = 0;
    M.Name = '';
    M.LongName = '';
    % JU: T is last for compatibility with older files modified in loadfieldhelper
    M.T = g.MeasureT;

    %Calc length
    dist = sum(sqrt(...
      (SET(no).ResolutionX*diff(M.X)).^2+...
      (SET(no).ResolutionY*diff(M.Y)).^2+...
      ((SET(no).SliceThickness+SET(no).SliceGap)*diff(M.Z)).^2));
    M.Length = dist;

    %Check if need name = new measure
    if isempty(g.MeasureName)
      [stri,lstr] = g.measureasklabel;
      if isempty(stri)
        %    myfailed('Invalid name.',g.GUI.Segment);
        %ind = true(1,length(SET(no).Measure));
        %ind(g.MeasureN) = false;
        %SET(no).Measure = SET(no).Measure(ind);
        drawfunctions('drawimageno');
        return;
      else
        M.Name = stri;
        M.LongName = lstr;
      end;
    else
      M.Name = g.MeasureName;
      M.LongName = g.MeasureName;
    end;

    %Commit to set
    if isempty(SET(no).Measure)
      SET(no).Measure = M;
    else
      SET(no).Measure(g.MeasureN).X = M.X;
      SET(no).Measure(g.MeasureN).Y = M.Y;
      SET(no).Measure(g.MeasureN).Z = M.Z;
      SET(no).Measure(g.MeasureN).Length = M.Length;
      SET(no).Measure(g.MeasureN).Name = M.Name;
      SET(no).Measure(g.MeasureN).T = M.T;
      SET(no).Measure(g.MeasureN).LongName = M.LongName;
    end;

    drawfunctions('drawimageno');
    segment('updateselectedslices');
    g.updateaxestables('angio');
    end
    
    %----------------------
    function contextmenu(g)
    %----------------------
    %Context menu button down / callback
    %Overloaded in CVQgui, RVQgui and Segment SPECT GUI
    % Note: Can act on one slice (CurrentSlice) or range of slices
    % (StartSlice:EndSlice)

    global SET NO

    set(g.imagefig,'WindowButtonMotionFcn','');
    set(g.imagefig,'WindowButtonUpFcn','');

    [p(1),p(2)] = mygetcurrentpoint(g.imagefig);
    [px,py] = mygetcurrentpoint(g.Handles.imageaxes(g.CurrentPanel));
    
%     % Contrast (and Crop, which has no menu) has total override
%     if ismember(g.CurrentTool,{'contrast'})
%       set(g.Handles.contrastcontextmenu,...
%         'Position',p,'Visible','on');
%       return;
%     end

    h = gco; %handle to current object
    
    %Look for objects in proximity of click
    col = 1+floor((px-0.5)/SET(NO).YSize);
    row = 1+floor((py-0.5)/SET(NO).XSize);
    px = px-(col-1)*SET(NO).YSize;
    py = py-(row-1)*SET(NO).XSize;
    if isequal(get(h,'Type'),'image')
      hchildren = get(get(h,'Parent'),'Children');
      for hc = hchildren'
        if isequal(get(hc,'Type'),'line') && ~isequal(get(hc,'Marker'),'none')
          x = get(hc,'XData');
          y = get(hc,'YData');
          contdist = min(sqrt((px-x).^2+(py-y).^2));
          if contdist < g.Pref.ContourAdjustDistance
            h = hc;
            break
          end
        end
      end
    end
    
    %Handles that override tools/themes
    if ismember(h,[
        g.Handles.measureline{g.CurrentPanel}{:} ...
        g.Handles.measuretext{g.CurrentPanel} ])
      set(g.Handles.measurecontextmenu,...
        'Position',p,'Visible','on');
      return;
    elseif ismember(h,[
        g.Handles.pointo{g.CurrentPanel} ...
        g.Handles.pointp{g.CurrentPanel} ...
        g.Handles.pointtext{g.CurrentPanel} ])
      set(g.Handles.pointcontextmenu,...
        'Position',p,'Visible','on');
      %set(g.Handles.pointcontextmenu,...
       % 'Position',p,'Visible','on');
      return;
    elseif ismember(h,g.Handles.endopin(g.CurrentPanel))
      set(g.Handles.endopinmenu,...
        'Position',p,'Visible','on');
      return;
    elseif ismember(h,g.Handles.epipin(g.CurrentPanel))
      set(g.Handles.epipinmenu,...
        'Position',p,'Visible','on');
      return;
    elseif ismember(h,g.Handles.rvendopin(g.CurrentPanel))
      set(g.Handles.rvendopinmenu,...
        'Position',p,'Visible','on');
      return;

    elseif ismember(h,[...
        g.Handles.endointerp(g.CurrentPanel)...
        g.Handles.epiinterp(g.CurrentPanel)...
        g.Handles.rvendointerp(g.CurrentPanel)...
        g.Handles.rvepiinterp(g.CurrentPanel)...
        ])
      set(g.Handles.interppointmenu,...
        'Position',p,'Visible','on');
      return;
    elseif ~isempty(SET(NO).EndoInterpX)&&...
        ~isempty(SET(NO).EndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice})&&...
        ismember(h,g.Handles.endocontour(g.CurrentPanel))
      set(g.Handles.interppointmenu,...
        'Position',p,'Visible','on');
      return;
    elseif ismember(h,[g.Handles.roicontour{:}])
      set(g.Handles.roicontextmenu,...
        'Position',p,'Visible','on');
      return;
    end

    %Klas: only one contextmenu for not clicking on object now
%     switch g.CurrentTheme
%       case misc
%         set(g.Handles.annocontextmenu,...
%              'Position',p,'Visible','on');
%       otherwise
        set(g.Handles.selectcontextmenu,...
        'Position',p,'Visible','on');
%    end
%     
%     % Tools that override themes: select/scale/move, measure, point, pins
%     if ismember(g.CurrentTool,{...
%         'select' ...
%         'move' ...
%         'scale' })
%       % CVQ has THEMES overriding these. /JU
%       set(g.Handles.selectcontextmenu,...
%         'Position',p,'Visible','on');
%       return;
%     elseif ismember(g.CurrentTool,{'measure'})
%       set(g.Handles.measurecontextmenu,...
%         'Position',p,'Visible','on');
%       return;
%     elseif ismember(g.CurrentTool,{'point'})
%       set(g.Handles.pointcontextmenu,...
%         'Position',p,'Visible','on');
%       return;
%     end
%     %elseif ismember(g.CurrentTool,{'endopin'})
%     %  set(g.Handles.endopinmenu,...
%     %    'Position',p,'Visible','on');
%     %   return;
%     %elseif ismember(g.CurrentTool,{'epipin'})
%     %  set(g.Handles.epipinmenu,...
%     %    'Position',p,'Visible','on');
%     %   return;
%     %elseif ismember(g.CurrentTool,{'rvendopin'})
%     %  set(g.Handles.rvendopinmenu,...
%     %    'Position',p,'Visible','on');
%     %   return;
%     %elseif ismember(g.CurrentTool,{'rvepipin'})
%     %  set(g.Handles.rvepipinmenu,...
%     %    'Position',p,'Visible','on');
%     %   return;
%     %end
% 
%     %Themes
%     switch g.CurrentTheme
%       case 'roi'
%         set(g.Handles.roicontextmenu,...
%           'Position',p,'Visible','on');
%       case 'lv'
%         %Previously if not(isempty(SET(NO).Flow)) -> roicontext
%         set(g.Handles.lvcontextmenu,...
%           'Position',p,'Visible','on');
%       case 'rv'
%         set(g.Handles.rvcontextmenu,...
%           'Position',p,'Visible','on');
%       case 'scar'
%         set(g.Handles.scarcontextmenu,...
%           'Position',p,'Visible','on');
%       case 'misc'
%         % In Segment, all cases leading to misc are handled above
%         % Can get here in specialGUI.
%         set(g.Handles.annocontextmenu,...
%           'Position',p,'Visible','on');
% 
% 
%     end
    end
    
    %-------------------------
    function setprefhandles(g)
    %-------------------------
    %Called by segpref (main). Overloaded in CVQgui since uses another fig
    %Also overloaded in Segment CMR and Segment CT
    figname = 'segpref';
    fig = openfig(figname,'reuse');
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
  
    % Generate a structure of handles to pass to callbacks, and store it.
    g.PrefHandles = guihandles(fig);
    g.PrefHandles.fig = fig;
    
    end
    
    %----------------------------
    function updateprefhandles(g)
    %----------------------------
    %Updates PrefHandles. Called by segpref('update')
    %Overloaded in CVQgui and RVQgui.
    %Default folder locations
    set(g.PrefHandles.datapathtext,'String',g.Pref.datapath);
    set(g.PrefHandles.exportpathtext,'String',g.Pref.exportpath);
    set(g.PrefHandles.cdpathtext,'String',g.Pref.CDPath);

    %Drawing/Viewing
    set(g.PrefHandles.addpointscheckbox,'Value',g.Pref.AddPoints);
    set(g.PrefHandles.blackwhitecheckbox,'Value',g.Pref.BlackWhite);
    set(g.PrefHandles.linewidthedit,'String',sprintf('%0.5g',g.Pref.LineWidth));
    set(g.PrefHandles.markersizeedit,'String',sprintf('%0.5g',g.Pref.MarkerSize));
    set(g.PrefHandles.contouradjustdistanceedit,'String',sprintf('%0.5g',g.Pref.ContourAdjustDistance));
    set(g.PrefHandles.numpointsedit,'String',sprintf('%d',g.Pref.NumPoints));
    set(g.PrefHandles.numpointstext,'String',sprintf('%d',g.NumPoints));
    set(g.PrefHandles.interppointsedit,'String',sprintf('%d',g.Pref.NumInterpPoints));
    set(g.PrefHandles.numberthumbnailsedit,'String',sprintf('%0.5g',g.Pref.NumberVisibleThumbnails));
    set(g.PrefHandles.anonymcheckbox,'Value',g.Pref.AnonymMode);
    set(g.PrefHandles.viewinterpolatedcheckbox,'Value',g.Pref.ViewInterpolated);
    set(g.PrefHandles.bgcolorcheckbox,'Value',g.Pref.BackgroundColor);
    
    switch g.Pref.GUIBackgroundColor(1)
      case 0.94
        set(g.PrefHandles.backgroundcolorpopupmenu,'Value',1);
      case 0.34
        set(g.PrefHandles.backgroundcolorpopupmenu,'Value',2);
      case 0
        set(g.PrefHandles.backgroundcolorpopupmenu,'Value',3);
    end
    
    %Analysis
    set(g.PrefHandles.endocenterradiobutton,'Value',g.Pref.EndoCenter);
    set(g.PrefHandles.epicenterradiobutton,'Value',not(g.Pref.EndoCenter));
    set(g.PrefHandles.radialprofilesedit,'String',sprintf('%d',g.Pref.RadialProfiles));
    set(g.PrefHandles.includeallpixelsinroicheckbox,'value',g.Pref.IncludeAllPixelsInRoi);
    set(g.PrefHandles.uselightcheckbox,'Value',g.Pref.UseLight);

    %System
    set(g.PrefHandles.hidefilesunixcheckbox,'Value',g.Pref.HideFilesUnix); % /JT ticket #412
    set(g.PrefHandles.showseriesdescriptioncheckbox, 'Value', g.Pref.ShowSeriesDescription); % /JT #411
    set(g.PrefHandles.donotaskcheckbox,'Value',g.Pref.DoNotAsk);
    set(g.PrefHandles.fastpreviewloadcheckbox,'Value',g.Pref.FastPreviewLoad);
    set(g.PrefHandles.checkversioncheckbox,'Value',g.Pref.CheckVersion);
    set(g.PrefHandles.useproxyservercheckbox,'Value',g.Pref.UseProxyServer);
    if isequal(g.Pref.WebBrowser,'explorer')
      set(g.PrefHandles.webbrowserpopupmenu,'value',1);
    elseif isequal(g.Pref.WebBrowser,[getenv('ProgramFiles') '\Mozilla Firefox\firefox.exe']);
      set(g.PrefHandles.webbrowserpopupmenu,'value',2);
    else
      set(g.PrefHandles.webbrowserpopupmenu,'value',3);
    end
    end
    
    %------------------------------------
    function updateprefhandlesadvanced(g)
    %------------------------------------
    %Called by segpref('updateadvanced'). Overloaded in SegmentGUI and 
    %RVQGUI.
    
    %DICOM COmmunication Settings
    set(g.PrefHandlesAdvanced.imagebasepathtext,'String',g.Pref.Pacs.ImageBasePath);
    set(g.PrefHandlesAdvanced.tempstoragepathtext,'String',g.Pref.Pacs.TempStoragePath);
    set(g.PrefHandlesAdvanced.pafreportfoldertext,'String',g.Pref.Pacs.PAFPathname);
    set(g.PrefHandlesAdvanced.reportfoldertext,'String',g.Pref.Pacs.ReportsheetPath);
    set(g.PrefHandlesAdvanced.dicomportedit,'String',g.Pref.Server.DICOMPort);
    set(g.PrefHandlesAdvanced.aetitleedit,'String',g.Pref.Server.AETitle);
    set(g.PrefHandlesAdvanced.sendoptionsedit,'String',g.Pref.Pacs.SendOptions);
    set(g.PrefHandlesAdvanced.switchtagscheckbox,'Value',g.Pref.Pacs.SwitchTags);
    set(g.PrefHandlesAdvanced.receiveoptionsedit,'String',g.Pref.Server.ReceiveOptions);

    %DICOM Interpretation
    set(g.PrefHandlesAdvanced.fasterpreviewcheckbox,'Value',g.Pref.Dicom.FasterPreview);
    set(g.PrefHandlesAdvanced.force16bitcheckbox,'Value',g.Pref.Dicom.Force16Bit);
    if isempty(g.Pref.Dicom.NormalizePhase)
      set(g.PrefHandlesAdvanced.normalizephaseask,'Value',1);
      set(g.PrefHandlesAdvanced.normalizephaseno,'Value',0);
      set(g.PrefHandlesAdvanced.normalizephaseyes,'Value',0);
    else
      set(g.PrefHandlesAdvanced.normalizephaseask,'Value',0);
      set(g.PrefHandlesAdvanced.normalizephaseno,'Value',~g.Pref.Dicom.NormalizePhase);
      set(g.PrefHandlesAdvanced.normalizephaseyes,'Value',g.Pref.Dicom.NormalizePhase);
    end
  
    set(g.PrefHandlesAdvanced.userloggingcheckbox,'Value',g.Pref.UserLogging);
    val = mygetvalue(g.PrefHandlesAdvanced.userloggingcheckbox);
    
    if val==1
      set(g.PrefHandlesAdvanced.setuserlogpathpushbutton,'Enable','on');
    else
      set(g.PrefHandlesAdvanced.setuserlogpathpushbutton,'Enable','off');
    end
      set(g.PrefHandlesAdvanced.userlogpathtext,'String',g.Pref.UserLogPath);
    end
    
    %----------------------
    function defaultpref(g)
    %----------------------
    %Sets Pref to default values. Called by segpref('default_Callback')
    %Overloaded in CVQgui, RVQgui, Segment CMRgui.

    %fodle location
    g.Pref.datapath = pwd;
    g.Pref.exportpath = pwd;
    g.Pref.CDPath = 'F:';

    %drawing/viewing
    g.Pref.AddPoints = false;
    g.Pref.BlackWhite = false;
    g.Pref.LineWidth = 1;
    g.Pref.MarkerSize = 5;
    g.Pref.ContourAdjustDistance = 3;
    g.Pref.NumPoints = 80;
    g.Pref.NumInterpPoints = 15;
    g.Pref.NumberVisibleThumbnails = 7;
    g.Pref.AnonymMode = false;
    g.Pref.ViewInterpolated = false;
    g.Pref.UserLogging=0;
    g.Pref.UserLogPath=[];

    %Analysis
    g.Pref.EndoCenter = true;
    g.Pref.RadialProfiles = 80;
    g.Pref.IncludeAllPixelsInRoi = false;
    g.Pref.UseLight = false;

    %System
    g.Pref.HideFilesUnix = true;
    g.Pref.AllowDicomCache = 0;
    g.Pref.ShowSeriesDescription = false;
    g.Pref.DoNotAsk = false;
    g.Pref.WebBrowser = 'explorer';
    g.Pref.FastPreviewLoad = true;
    g.Pref.CheckVersion = true;
    
    %Communication
    g.Pref.Server.ReceiveOptions = ''; %'--prefer-little'; %This is used for storescp to configure.
    g.Pref.Pacs.SendOptions = ''; %This is used for storescu to configure.
    % g.Pref.LearnMode = false;
    % g.Pref.UndoHistory = 10;
    % g.Pref.AutoSave = false;
    end
    
    %-------------------------------
    function normalizephaseupdate(g)
    %-------------------------------
    %Called by segpref.m methods
    %Overloaded in CVQgui (different handle placement)
    handle = g.PrefHandlesAdvanced;

    if isempty(g.Pref.Dicom.NormalizePhase)
      set(handle.normalizephaseask,'Value',1);
      set(handle.normalizephaseno,'Value',0);
      set(handle.normalizephaseyes,'Value',0);
    else
      set(handle.normalizephaseask,'Value',0);
      set(handle.normalizephaseno,'Value',~g.Pref.Dicom.NormalizePhase);
      set(handle.normalizephaseyes,'Value',g.Pref.Dicom.NormalizePhase);
    end
    end
    
    %-----------------------------------
    function varargout = initopenfile(g)
    %-----------------------------------
    %Initializes the openfile GUI. Optional output is the fig
    %Called by openfile (main). Overloaded in CVQgui
    if isopengui('openfile.fig')
      openfile('exitgui')
    end
    g.GUI.OpenFile = mygui('openfile.fig'); 
    fig=g.GUI.OpenFile.fig;
    gui=g.GUI.OpenFile;
    
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

    % Generate a structure of handles to pass to callbacks, and store it.
    g.Preview.Handles = g.GUI.OpenFile.handles; %guihandles(fig);
    g.Preview.Handles.fig = fig;
    g.Preview.Silent = false;
    g.Preview.CanCrop = false;

    %Horisonal adjustment
    set(fig,'units','pixels');
    
    load('newicons.mat','newicons')
    gui.Icons=newicons;
    gui.iconholder = myiconplaceholder(gui.handles.iconaxes,0,3,fig);
    openfile('initiconholder');
    set(gui.fig,'WindowButtonMotionFcn','openfile(''motionfunc'')');
%     %Add some images
%     load('icons.mat')
%     set(g.Preview.Handles.browsepushbutton,'CData',icon.opendoc);
%     set(g.Preview.Handles.updirpushbutton,'CData',icon.updir);
%     set(g.Preview.Handles.refreshpushbutton,'CData',icon.refresh);

    %Fix with preview axis
    try
      axes(g.Preview.Handles.previewaxes);
    catch %#ok<CTCH>
      g.Preview.Handles.previewaxes = gca;
    end;
    g.Preview.Handles.image = image(repmat(0,[1 1 3]));
    g.Preview.PreviewFile = '';
    axis off;

    %Check directory and fill pathlistbox
    if ~exist(g.Preview.PathName,'dir')
      g.Preview.PathName = pwd;
      openfile('browsebutton_Callback');
    end;

    %set(g.Preview.Handles.fig,'keypressfcn',@keypressed);
    %set(g.Preview.Handles.pathlistbox,'keypressfcn',@keypressed);

    %Scans directory and update lister
    openfile('getpathinfo'); %Updates g.Preview.FileList
    if length(g.Preview.PathName)>130
      set(g.Preview.Handles.pathtext,'String',['... ' g.Preview.PathName((end-130):end)]);
    else
      set(g.Preview.Handles.pathtext,'String',g.Preview.PathName);
    end;
    
    [types,viewplanes,techniques] = g.imagedescription;
    g.ImageTypes = types;
    g.ImageViewPlanes = viewplanes;
    g.ImagingTechniques = techniques;
    if isfield(g.Preview.Handles,'imagingtechniquelistbox')
      set(g.Preview.Handles.imagingtechniquelistbox,'String',techniques);
      set(g.Preview.Handles.imagetypelistbox,'String',types);
      set(g.Preview.Handles.imageviewplanelistbox,'String',viewplanes);
    end

    g.Preview.ImageType = 'General';
    g.Preview.ImageViewPlane = 'Unspecified';
    g.Preview.ImagingTechnique = 'Unspecified';

    %Update region of interrest selection button.
    openfile('roisizelistbox_Callback', g.Preview.ROISize);

    %Update stable selection
    %set(g.Preview.Handles.stablecheckbox,'value',double(g.Preview.Stable));

    if nargout > 0
      varargout{1} = fig;
    end

    end
    
    %-------------------------
    function enableopenfile(g)
    %-------------------------
    %Called by openfile('enablesegmentgui'). Overloaded in CVQgui and RVQgui.

set([... %g.Handles.reportmenu ... g.Handles.t2starmenu ...g.Handles.t1analysismenu ...g.Handles.t2analysismenu ...
      g.Handles.measuremenu ...
      g.Handles.fileloadnextmenu ...
      g.Handles.toolsmenu ...
      g.Handles.roimenu ...
      g.Handles.editmenu ...
      g.Handles.exportmenu ...
      g.Handles.segmentmenu ...
      g.Handles.mrmenu ...
      g.Handles.viabilitymenu ...
      g.Handles.flowmenu ...
      g.Handles.ctmenu ...
      g.Handles.spectmenu ...
      g.Handles.viewmenu ...
      g.Handles.lvmenu ...
      g.Handles.rvmenu ...
      g.Handles.analysismenu ...
      g.Handles.fileloadsegmentationmenu ...
      g.Handles.fileloadnextmenu ...
      g.Handles.filesavesegmentationmenu ...
      g.Handles.filesavesubmenu ...
      g.Handles.filesavecurrentmenu ...
      g.Handles.filesaveallmenu ...
      g.Handles.filesaveallasmenu ...
      g.Handles.filesavetodatabasemenu ...
      g.Handles.filesavetopacsmenu ...
      g.Handles.filesavesegdicom ...
      g.Handles.fileclosecurrentimagestack ...
      g.Handles.filecloseallmenu ...
      g.Handles.fileclosemultiplemenu ...
      g.Handles.thumbnailslider...
      ],...
      'enable','on');
      
    set([...
      g.Handles.reportpanel ...
      g.Handles.barpanel ... %g.Handles.thistimeframeonlycheckbox ...%g.Handles.excludepapilarscheckbox ...%g.Handles.uselightcheckbox ...
      g.Handles.thumbnailslider...
      ],...
      'visible','on');
   
    end
    
    %-----------------------
    function checkversion(g) %#ok<MANU>
    %-----------------------
    %Check if new version is available. Define separately for each GUI.
    end
    
    %------------------------------------------
    function updateaxestables(g, arg, varargin) %#ok<INUSD>
    %------------------------------------------
    %Method to update AxesTables, in GUIs where present
    end
    
    %--------------------------------
    function roiputroi_helper(g,no,m)
    %--------------------------------
    %Called by roi('roiputroi_Buttondown')
    %Overloaded in CVQgui and RVQgui
    global SET DATA
    
    oldpref=DATA.ThisFrameOnly;
    
    if  ~isempty(SET(no).Flow)
     DATA.ThisFrameOnly = false;
    end
    
    [y,x,slice] = segment('getclickedcoords');
    
    o = linspace(0,2*pi,g.NumPoints)';
    SET(no).Roi(m).X = min(SET(no).XSize+0.5,max(0.5,...
      repmat(x+8*sin(o),1,SET(no).TSize) ));
    SET(no).Roi(m).Y = min(SET(no).YSize+0.5,max(0.5,...
      repmat(y+8*cos(o),1,SET(no).TSize) ));
    
    %Find name   
    if SET(no).RoiN==1 || isempty(SET(no).RoiCurrent)
      if isequal(SET(no).ImageViewPlane,'Aorta')
        SET(no).Roi(m).Name = 'Aortic ascending flow';
      else
        tempN = SET(no).RoiN; while ismember(sprintf('ROI-%d',tempN),{SET(no).Roi(1:end-1).Name}), tempN = tempN+1; end
        SET(no).Roi(m).Name = sprintf('ROI-%d',tempN);        
      end;
      SET(no).Roi(m).LineSpec = 'b-';
    else
      if length(SET(no).Roi(SET(no).RoiCurrent(end)).Name)>=4 && isequal(SET(no).Roi(SET(no).RoiCurrent(end)).Name(1:4),'ROI-')
        tempN = SET(no).RoiN; while ismember(sprintf('ROI-%d',tempN),{SET(no).Roi(1:end-1).Name}), tempN = tempN+1; end
        SET(no).Roi(m).Name = sprintf('ROI-%d',tempN);
      else
        SET(no).Roi(m).Name =SET(no).Roi(SET(no).RoiCurrent(end)).Name;
      end
      SET(no).Roi(m).LineSpec = SET(no).Roi(SET(no).RoiCurrent(end)).LineSpec;
    end;
    
    if g.ThisFrameOnly
      SET(no).Roi(m).T = SET(no).CurrentTimeFrame;
    else
      SET(no).Roi(m).T = 1:SET(no).TSize;
    end
    SET(no).Roi(m).Z = slice;
    SET(no).Roi(m).Sign = 1;
    %SET(no).Roi(m).LineSpec = g.GUISettings.DefaultROISpec;

    %define if ROI should updated flow result
    if ~isempty(SET(no).Flow)
      SET(no).Roi(m).Flow = true;
    else
      SET(no).Roi(m).Flow = false;
    end
    
    %Make as current
    SET(no).RoiCurrent = m;
    
    if oldpref~=DATA.ThisFrameOnly;
      DATA.ThisFrameOnly=oldpref;
    end
    end
    
    %--------------------------------------------------
    function name = roilabelmenu(g,roitoname,roinamein)
    %--------------------------------------------------
    %Prompt name of ROI from a menu selection
    %Overloaded in CVQ and RVQ
    
    global SET
    
    name = '';

    useroinname=true;
    if nargin < 2 || isempty(roitoname)
      useroinname=false;
    end
    
    temp = {...
      'Papillary',...
      'Remote ROI',...
      'Scar region ROI',...
      'Static tissue',...
      'Non-static tissue',...
      'Aortic ascending flow',...
      'Aortic descending flow',...
      'Abdominal Aorta',...
      'Pulmonary artery',...
      'Vena cava inf',...
      'Vena cava sup',...
      'ULPV', ...
      'LLPV', ...
      'URPV', ...
      'LRPV', ...
      'LPA', ...
      'RPA', ...
      'Sinus coronarius',...
      'Lung',...
      'Heart',...
      'Blood',...
      'Aortic valve',...
      'Pulmonary valve',...
      'Mitral valve',...
      'Tricuspid valve',...
      'Left atrium',...
      'Right atrium',...
      'Left ventricle',...
      'Right ventricle',...
      'Coronary sinus rest',...
      'Coronary sinus stress'};
    
    no=roi('roifindmag');
    
    if SET(no).RoiN>0
      roinames=cell(1,SET(no).RoiN);
      for rloop=1:SET(no).RoiN
        roinames{rloop}=SET(no).Roi(rloop).Name;
      end
      temp = cat(2,union(temp,roinames));
    end;

    if useroinname
      temp = cat(2,'User defined ...','ROI-n',temp);
    else
      temp = cat(2,'User defined ...',temp);
    end

    if nargin<3
      m = mymenu('Select a name',temp,g.GUI.Segment);
    else
      m = mymenu(strcat([translation.dictionary('Select a new name for "') roinamein '"']),temp,g.GUI.Segment);
    end

    if useroinname
      if m==1
        if iscell(roitoname)
          name = inputdlg({'Enter name for ROI'},...
            'New name',...
            1,...
            {'ROI'});
        else
          name = inputdlg({'Enter name for ROI'},...
            'New name',...
            1,...
            {sprintf('ROI-%d',roitoname)});
        end
        if isempty(name)
          name = sprintf('ROI-%d',roitoname);
        else
          name = name{1};
        end;
      end;

      if m==2
        if iscell(roitoname)
          for loop = 1:length(roitoname)
            name{loop}=sprintf('ROI-%d',roitoname{loop});
          end
        else
          name = sprintf('ROI-%d',roitoname);
        end
      end;

      if m>2
        name = temp{m};
      end;
    else
      if m==1
        name = inputdlg({'Enter name for ROI'},...
          'New name',...
          1,...
          {'ROI'});
        if isempty(name)
          name = 'ROI';
        else
          name = name{1};
        end;
      end;

      if m>1
        name = temp{m};
      end;
    end
    
    end
    
    %----------------------------
    function updatetimebaraxes(g)
    %----------------------------
    %update time bar axes in main interface
    global DATA SET NO
    
    if ~isempty(DATA.FlowNO) && ismember(NO,SET(DATA.FlowNO).Linked)
      no = DATA.FlowNO;
    else
      no = NO;
    end
    
    if isempty(no) || SET(no).TSize<2
      cla(g.Handles.timebaraxes);
      set(g.Handles.timebaraxes,'Visible','off')
      g.Handles.timebar = [];
      g.Handles.timebaraxeshelpers=[];
      g.Handles.edtimebartext=[];
      g.Handles.estimebartext=[];
      g.Handles.edtimebarline=[];
      g.Handles.estimebarline=[];
    else
      
      set(g.Handles.timebaraxes,'Visible','on')
      %this is speed problem ... try to fix by checking if 
      lc = DATA.GUISettings.BarColor;%[55 119 106]/255; %[1 0.47 0]; %[1 0.6 0.2]; % [0.5 0.5 0];    %JU
      t = SET(no).TimeVector*1000;
      if ~isfield(g.Handles,'timebar') || isempty(g.Handles.timebar)
      g.Handles.timebar = plot(g.Handles.timebaraxes,...
        [t(SET(no).CurrentTimeFrame) t(SET(no).CurrentTimeFrame)],...
        get(g.Handles.timebaraxes,'ylim'),'-','Color',lc);
      set(g.Handles.timebar,'linewidth',3);
      else
        set(g.Handles.timebar,'xdata',[t(SET(no).CurrentTimeFrame) t(SET(no).CurrentTimeFrame)],'color',lc);
      end
      hold(g.Handles.timebaraxes,'on');
      temp = get(g.Handles.timebaraxes,'ylim');
      ttemp = [t;t];
      temp = [repmat(temp(1),size(t));repmat(temp(1)+0.1*(temp(2)-temp(1)),size(t))];
      
       if  length(g.Handles.timebaraxeshelpers)~=length(t)
        delete(g.Handles.timebaraxeshelpers)
        g.Handles.timebaraxeshelpers=[];
      end
      
      if ~isfield(g.Handles,'timebaraxeshelpers') || isempty(g.Handles.timebaraxeshelpers)
        g.Handles.timebaraxeshelpers=plot(g.Handles.timebaraxes,ttemp,temp,'k-');
        set(g.Handles.timebaraxes,'YTick',[]);
      else
        instr=ones(1,length(t));
        xcell=mat2cell(ttemp,2,instr)';
        ycell=mat2cell(temp,2,instr)';
        set(g.Handles.timebaraxeshelpers,{'XData'},xcell,{'YData'},ycell)
      end
      
      set(g.Handles.timebaraxes,'ButtonDownFcn','segment(''timebaraxes_Buttondown'')');
      
      %ED and ES
      temp = get(g.Handles.timebaraxes,'ylim');
      tempmax = temp(end); %JU
      tempmin = temp(1);   %JU
      
      if ~isfield(g.Handles,'edtimebartext') || isempty(g.Handles.edtimebartext)
      g.Handles.edtimebartext = text(...
        'position',[t(SET(no).EDT)-t(end)*0.015 tempmax-0.2*(tempmax-tempmin)],...
        'string','ED',...
        'parent',g.Handles.timebaraxes,...
        'color',lc, ...
        'FontWeight','bold');
      else
        set(g.Handles.edtimebartext,'position',[t(SET(no).EDT)-t(end)*0.015 tempmax-0.2*(tempmax-tempmin)],'color',lc)  
      end
      
      if ~isfield(g.Handles,'edtimebarline') || isempty(g.Handles.edtimebarline)
      g.Handles.edtimebarline = plot(g.Handles.timebaraxes,...
        [t(SET(no).EDT) t(SET(no).EDT)],[tempmax tempmax*0.8],'color',lc);
      else
        set(g.Handles.edtimebarline,'xdata',[t(SET(no).EDT) t(SET(no).EDT)],'ydata',[tempmax tempmax*0.8],'color',lc); 
      end
      
      set([g.Handles.edtimebartext g.Handles.edtimebarline],'ButtonDownFcn','segment(''esedtimebar_Buttondown'',''ed'')','color',lc);
      
      if ~isfield(g.Handles,'estimebartext') || isempty(g.Handles.estimebartext)
      g.Handles.estimebartext = text(...
        'parent',g.Handles.timebaraxes,...
        'position',[t(SET(no).EST)-0.015*t(end) tempmin+0.2*(tempmax-tempmin)],...
        'string','ES',...
        'color',lc, ...
        'FontWeight','bold');
      else
        set(g.Handles.estimebartext,'position',[t(SET(no).EST)-0.015*t(end) tempmin+0.2*(tempmax-tempmin)],'color',lc) 
      end
      
      if ~isfield(g.Handles,'estimebarline') || isempty(g.Handles.estimebarline)
        g.Handles.estimebarline = plot(g.Handles.timebaraxes,[t(SET(no).EST) t(SET(no).EST)],...
          [tempmin tempmin+(0.2*(tempmax-tempmin))],'color',lc);
      else
        set(g.Handles.estimebarline,'xdata',[t(SET(no).EST) t(SET(no).EST)],'ydata',[tempmin tempmin+(0.2*(tempmax-tempmin))],'color',lc); 
      end
      
      
      hold(g.Handles.timebaraxes,'off');
      set([g.Handles.estimebartext g.Handles.estimebarline],'ButtonDownFcn','segment(''esedtimebar_Buttondown'',''es'')');
      set(g.Handles.timebar,'buttondownFcn','segment(''timebar_Buttondown'')');
      set(g.Handles.timebaraxes,'xlim',[t(1) t(end)]); 
      
      xlabel(g.Handles.timebaraxes,translation.dictionary('Time [ms]'),'color',g.GUISettings.TimebarAxesColor);
       set(g.Handles.timebaraxes,...
         'XColor',g.GUISettings.TimebarAxesColor,...
         'YColor',g.GUISettings.TimebarAxesColor);
      set(g.Handles.timebaraxes, ...
        'Color',g.GUISettings.VolumeColorGraph);
    
    if no == DATA.LVNO
      set(DATA.Handles.timebarlv,'xdata',...
        t(SET(no).CurrentTimeFrame)*[1 1])
    end
    if no == DATA.FlowNO
      set(DATA.Handles.timebarflow,'xdata',...
        t(SET(no).CurrentTimeFrame)*[1 1])
    end
    end
    end
    
    
    %---------------------------
    function updatevolumeaxes(g)
    %---------------------------
    %Called by segment_main('updatevolume').
    global DATA SET
    
    if ~isempty(DATA.LVNO)
      no = DATA.LVNO;
    elseif ~isempty(DATA.RVNO)
      no = DATA.RVNO;
    else
      no = [];
    end
    
    if isempty(no) || SET(no).TSize<2
      cla(g.Handles.volumeaxes);
      g.Handles.timebarlv = [];
      g.Handles.volumecurve = [];
      g.Handles.masscurve = [];
      g.Handles.rvvcurve = [];
      g.Handles.volumeaxeshelpers = [];
      g.Handles.estext=[];
      g.Handles.edtext=[];
      g.Handles.esline=[];
      g.Handles.edline=[]; 
      g.Handles.lvmtext=[];
      set(g.Handles.volumeaxes,'Visible','off');
    else
      lc = DATA.GUISettings.BarColor;%[55 119 106]/255;
      set(g.Handles.volumeaxes,'Visible','on');
      %LV and RV curve and lines
      t = SET(no).TimeVector*1000;
      oldvolpeak=0;
      oldrvpeak=0;
      ylim(g.Handles.volumeaxes,'auto')

      %This tells us if bar removal is needed. Needed being that if masscurve max is changed or the borders have changed
      if ~isempty(g.Handles.masscurve)
          %oldlvm=max(get(g.Handles.masscurve,'ydata'));
      end
      yborders=get(g.Handles.volumeaxes,'ylim');
      yroofprecurves=yborders(2);
      
      if isempty(g.Handles.volumecurve)||~ishandle(g.Handles.volumecurve)
          g.Handles.volumecurve = plot(g.Handles.volumeaxes,t,SET(no).LVV-SET(no).PV,'r.-');
          set(g.Handles.volumecurve,'markersize',5);
      else
          oldvolpeak=max(get(g.Handles.volumecurve,'ydata'));
          set(g.Handles.volumecurve,'xdata',t,'ydata',SET(no).LVV-SET(no).PV)
      end
  
      newvolpeak=max(get(g.Handles.volumecurve,'ydata'));
      hold(g.Handles.volumeaxes,'on');
      lvmplot = SET(no).EPV-SET(no).LVV+SET(no).PV;
      if sum(~isnan(lvmplot)) == length(lvmplot)
        linespec = 'b-';
        markersize = 3;
      else
        linespec = 'b.-';
        markersize = 5;
      end
      
      if isempty(g.Handles.masscurve) || ~ishandle(g.Handles.masscurve)
        g.Handles.masscurve = plot(g.Handles.volumeaxes,t,lvmplot,linespec);
        set(g.Handles.masscurve,'markersize',markersize);
      else
        set(g.Handles.masscurve,'xdata',t,'ydata',lvmplot)
      end
      
      if   isempty(g.Handles.rvvcurve) || ~ishandle(g.Handles.rvvcurve)
        %if ~isfield(g.Handles,'rvvcurve') || isempty(g.Handles.rvvcurve)
        %          if ishandle(g.Handles.rvvcurve)
        g.Handles.rvvcurve = plot(g.Handles.volumeaxes,t,SET(no).RVV,'m.-');
        set(g.Handles.rvvcurve,'markersize',5);
        
      else
        oldrvpeak=max(get(g.Handles.rvvcurve,'ydata'));
        set(g.Handles.rvvcurve,'xdata',t,'ydata',SET(no).RVV);
      end
      
      newrvpeak=max(get(g.Handles.rvvcurve,'ydata'));
      yborders=get(g.Handles.volumeaxes,'ylim');
      yroofpostcurves=yborders(2);
      
      if newvolpeak~=oldvolpeak || newrvpeak~=oldrvpeak || yroofprecurves~=yroofpostcurves || yborders(1)~=0
        %removeallbars
        delete(g.Handles.timebarlv);
        delete(g.Handles.volumeaxeshelpers);
        delete(g.Handles.estext);
        delete(g.Handles.edtext);
        delete(g.Handles.esline);
        delete(g.Handles.edline);
        delete(g.Handles.lvmtext)
        
        g.Handles.timebarlv=[];
        g.Handles.volumeaxeshelpers = [];
        g.Handles.estext=[];
        g.Handles.edtext=[];
        g.Handles.esline=[];
        g.Handles.edline=[];
        g.Handles.lvmtext=[];
      end
        
      %Texts and
      %bars-----------------------------------------------------------------------------
      
      % here we do not want to mess with upper limit this should be set by
      % the curves plotted so fixate the top limit
      temp = get(g.Handles.volumeaxes,'ylim');
      
      ylim(g.Handles.volumeaxes,'manual')
      ylim(g.Handles.volumeaxes,[0 temp(2)])
      ttemp = [t;t];
      temp = [zeros(size(t)); repmat(0.1*(temp(2)-temp(1)),size(t))];
      
      if  length(g.Handles.volumeaxeshelpers)~=length(t)
        delete(g.Handles.volumeaxeshelpers)
        g.Handles.volumeaxeshelpers=[];
      end
      
      if isempty(g.Handles.volumeaxeshelpers) || any(~ishandle(g.Handles.volumeaxeshelpers))
        g.Handles.volumeaxeshelpers=plot(g.Handles.volumeaxes,ttemp,temp,'k-');
        %set(g.Handles.timebaraxes,'YTick',[]);
      else
%         if length(t)~=length(g.Handles.volumeaxeshelpers)%find(~isempty(g.Handles.volumeaxeshelpers)))
%           delete(g.Handles.volumeaxeshelpers(length(t)+1:end))
%         end
        instr=ones(1,length(t));
        xcell=mat2cell(ttemp,2,instr)';
        ycell=mat2cell(temp,2,instr)';
        set(g.Handles.volumeaxeshelpers,{'XData'},xcell,{'YData'},ycell)
      end
      
      if isempty(g.Handles.timebarlv) || ~ishandle(g.Handles.timebarlv)
        g.Handles.timebarlv = plot(g.Handles.volumeaxes,...
          [t(SET(no).CurrentTimeFrame) t(SET(no).CurrentTimeFrame)],...
          get(g.Handles.volumeaxes,'ylim'),'-','Color',lc);
        set(g.Handles.timebarlv,'linewidth',2);
      else
        set(g.Handles.timebarlv,'xdata',[t(SET(no).CurrentTimeFrame) t(SET(no).CurrentTimeFrame)],'ydata',get(g.Handles.volumeaxes,'ylim'),'color',lc)
      end      
      %plot(g.Handles.volumeaxes,ttemp,temp,'k-');
      
      %ED and ES
      temp = get(g.Handles.volumeaxes,'ylim');
      tempmax = temp(end); %JU
      tempmin = temp(1);   %JU
      
      if isempty(g.Handles.edtext) || ~ishandle(g.Handles.edtext)
        g.Handles.edtext = text(...
          'position',[t(SET(no).EDT)-t(end)*0.05 tempmax-0.2*(tempmax-tempmin)],...
          'string','ED',...
          'parent',g.Handles.volumeaxes,...
          'color',lc,'FontWeight','bold');
      else
        set(g.Handles.edtext,'position',[t(SET(no).EDT)-t(end)*0.05 tempmax-0.2*(tempmax-tempmin)],'color',lc);
      end
      
      if isempty(g.Handles.edline) || ~ishandle(g.Handles.edline)
        g.Handles.edline = plot(g.Handles.volumeaxes,...
          [t(SET(no).EDT) t(SET(no).EDT)],[tempmax tempmax*0.8],'color',lc);
      else
        set(g.Handles.edline,'xdata',[t(SET(no).EDT) t(SET(no).EDT)],'ydata',[tempmax tempmax*0.8])
      end
      
      set([g.Handles.edtext g.Handles.edline],'ButtonDownFcn','segment(''esed_Buttondown'',''ed'')',...
        'color',lc,'linewidth',2);
      set(g.Handles.volumeaxes,'ButtonDownFcn','segment(''volumeaxes_Buttondown'')');
      
      
      if isempty(g.Handles.estext) || ~ishandle(g.Handles.estext)
        g.Handles.estext = text(...
          'parent',g.Handles.volumeaxes,...
          'position',[t(SET(no).EST)-0.05*t(end) tempmin+0.2*(tempmax-tempmin)],...
          'string','ES',...
          'color',lc,'FontWeight','bold');
      else
        set(g.Handles.estext,'position',[t(SET(no).EST)-0.05*t(end) tempmin+0.2*(tempmax-tempmin)],'color',lc);
      end
        
     if isempty(g.Handles.esline) || ~ishandle(g.Handles.esline) 
       g.Handles.esline = plot(g.Handles.volumeaxes,[t(SET(no).EST) t(SET(no).EST)],...
         [tempmin tempmin+(0.2*(tempmax-tempmin))],'color',lc,'linewidth',2);
     else
       set(g.Handles.esline,'xdata',[t(SET(no).EST) t(SET(no).EST)],'ydata',[tempmin tempmin+(0.2*(tempmax-tempmin))],'color',lc)
     end
     
      set([g.Handles.estext g.Handles.esline],'ButtonDownFcn','segment(''esed_Buttondown'',''es'')');
      maxlvm = max(SET(no).EPV-SET(no).LVV+SET(no).PV);
      if ~isnan(maxlvm)
        if ~isempty(g.Handles.lvmtext)
          set(g.Handles.lvmtext,'position',[0.85*t(end) maxlvm+0.06*(tempmax-tempmin)])
        else
          g.Handles.lvmtext = text(...
            'parent',g.Handles.volumeaxes,...
            'position',[0.85*t(end) maxlvm+0.06*(tempmax-tempmin)],...
            'string','LVM',...
            'color','blue','fontsize',8);
        end
      else
        delete(g.Handles.lvmtext)
        g.Handles.lvmtext=[];
      end
%------------------------------------------------------------------------------------------------------------
      
      hold(g.Handles.volumeaxes,'off');
      set(g.Handles.timebarlv,'buttondownFcn','segment(''timebarlv_Buttondown'')');
      set(g.Handles.volumeaxes,'xlim',[t(1) t(end)]); %,'ticklength',[0 0]);%1000*60/SET(no).HeartRate]);
      ylabel(g.Handles.volumeaxes,translation.dictionary('Volume [ml]'),'color',g.GUISettings.VolumeAxesColor);
      xlabel(g.Handles.volumeaxes,translation.dictionary('Time [ms]'),'color',g.GUISettings.VolumeAxesColor);
      set(g.Handles.volumeaxes,...
        'XColor',g.GUISettings.VolumeAxesColor,...
        'YColor',g.GUISettings.VolumeAxesColor);
      set(g.Handles.volumeaxes, ...
        'Color',g.GUISettings.VolumeColorGraph);
    end

    end
    
    
    %-------------------------
    function updateflowaxes(g)
    %-------------------------
    %Update the flow report axes in main GUI
    global DATA SET
    
    no = DATA.FlowNO;
    roinbr = DATA.FlowROI;
    
    if isempty(no) || SET(no).TSize<2 || isempty(SET(no).Flow) || isempty(SET(no).Flow.Result) || isempty(roinbr) || length(SET(no).Flow.Result)< roinbr || isempty(SET(no).Flow.Result(roinbr)) || isempty(SET(no).Flow.Result(roinbr).netflow) 
      cla(g.Handles.flowaxes);
      g.Handles.timebarflow = [];
      g.Handles.flowcurve = [];
      set(g.Handles.flowaxes,'Visible','off');
      %set(g.Handles.flowuipanel,'Visible','off');
      %set(g.Handles.flowresultaxes,'Visible','off');
    else
      lc = DATA.GUISettings.BarColor;%[55 119 106]/255;
      
      %set(g.Handles.flowuipanel,'Visible','on');
      set(g.Handles.flowaxes,'Visible','on');
      %set(g.Handles.flowresultaxes,'Visible','on');
      %Time resolved
      t = SET(no).TimeVector*1000;
      %flow curve
      netflow = SET(no).Flow.Result(roinbr).netflow;
      
      %plot in color of roi
      curvecolor=SET(no).Roi(roinbr).LineSpec(1);
      
      g.Handles.flowcurve = plot(g.Handles.flowaxes,t,netflow,[curvecolor,'.-']);%'b.-');
      hold(g.Handles.flowaxes,'on');
      set(g.Handles.flowcurve,'markersize',5);
      % set y-lim dependent on netflow values
      minflow = min(netflow);
      maxflow = max(netflow);
      set(g.Handles.flowaxes,'ylim',[floor((minflow-50)/100)*100 ceil((maxflow+50)/100)*100]);
      %plot line at zero
      plot(g.Handles.flowaxes,[t(1) t(end)],[0 0],'k:');
      %outer time bars
      plot(g.Handles.flowaxes,[t(1) t(1)],get(g.Handles.flowaxes,'ylim'),'color',lc)%'g-');
      plot(g.Handles.flowaxes,[t(end) t(end)],get(g.Handles.flowaxes,'ylim'),'color',lc)%'g-');
      %current time
      g.Handles.timebarflow = plot(g.Handles.flowaxes,...
        [t(SET(no).CurrentTimeFrame) t(SET(no).CurrentTimeFrame)],...
        get(g.Handles.flowaxes,'ylim'),'-','Color',lc);
      temp = get(g.Handles.flowaxes,'ylim');
      ttemp = [t;t];
      temp = [repmat(temp(1),size(t));repmat(temp(1)+0.1*(temp(2)-temp(1)),size(t))];
      plot(g.Handles.flowaxes,ttemp,temp,'k-');
      set(g.Handles.timebarflow,'linewidth',2);    
      set(g.Handles.timebarflow,'buttondownFcn','segment(''timebarflow_Buttondown'')');
      set(g.Handles.flowaxes,'xlim',[t(1) t(end)]); %,'ticklength',[0 0]);%1000*60/SET(no).HeartRate]);
      hold(g.Handles.flowaxes,'off');
      set(g.Handles.flowaxes,'ButtonDownFcn','segment(''flowaxes_Buttondown'')');      
      %labels and color
      ylabel(g.Handles.flowaxes,translation.dictionary('Flow [ml/s]'),'color',g.GUISettings.VolumeAxesColor);
      xlabel(g.Handles.flowaxes,translation.dictionary('Time [ms]'),'color',g.GUISettings.VolumeAxesColor);
      set(g.Handles.flowaxes,...
        'XColor',g.GUISettings.VolumeAxesColor,...
        'YColor',g.GUISettings.VolumeAxesColor);
      set(g.Handles.flowaxes, ...
        'Color',g.GUISettings.VolumeColorGraph);  
      if ~isempty(SET(SET(no).Flow.PhaseNo).Flow.PhaseCorr)
        if ~isfield(SET(SET(no).Flow.PhaseNo).Flow,'PhaseCorrTimeResolved')
          mywarning('Incompatible eddy current correction. Correction reset.',DATA.GUI.Segment);
        else
        g.Handles.flowtext = text(get(g.Handles.flowaxes,'XLim'),get(g.Handles.flowaxes,'YLim'),'[Eddy]', ...
          'Parent',g.Handles.flowaxes,'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','top');
        end
      else
%         g.Handles.flowtext = text(get(g.Handles.flowaxes,'XLim'),get(g.Handles.flowaxes,'YLim'),'hej', ...
%           'Parent',g.Handles.flowaxes,'FontSize',11,'HorizontalAlignment','right','VerticalAlignment','top');
      end
    end
    end
        
    %--------------------------------------
    function measurefontsize(g,panel,index)
    %--------------------------------------
    %Sets measure font size. Used in CVQgui.
    end
    
    %------------------------------------------------------------------
    function addmainicon_helper(g,callback,tooltip,cdata,tag,separator)
    %------------------------------------------------------------------
    %Helper function to add an icon

    if nargin<5
      myfailed('Too few input arguments.',g.GUI.Segment);
      return;
    end;

    if nargin<6
      separator = 'off';
    end;

    props = [];
    props.ClickedCallback = callback;
    props.ToolTip = translation.dictionary(tooltip);
    props.Tag = tag;
    props.CData = cdata;
    props.Separator = separator;
    g.Handles.(tag) = uitoggletool(g.Handles.maintoolbar,props);
    
    end
    
    %--------------------------------------
    function heartratezero(g, heartrateest)
    %--------------------------------------
    %React upon a heart rate read as zero from DICOM in openfile('initset')
    mywarning(dprintf(...
    'Specified Heart Rate is 0, assuming alive patient, guess on %0.5g (based on timeincrement)',...
    heartrateest),g.GUI.Segment);
    end
    
    %----------------------------------
    function drawsectorgrid(g,no,panel) 
    %----------------------------------
    %Overloaded in RVQ GUI
    end
    
    %----------------------------
    function buttonup_Callback(g) 
    %----------------------------
    %General buttonupfunction. Overloaded in RVQ GUI
    
    set(g.imagefig,'WindowButtonMotionFcn','');
    %drawfunctions('drawsliceno');
    end
      
    %-----------------------------------------
    function printthumbnailnumber(g,thumbsize)
    %-----------------------------------------
    %Overloaded in RVQ GUI
    global SET
    fontsize = thumbsize/23;
    
    g.Handles.datasetnumber = [];
   
    numletters=8;%thumbsize/fontsize;
    for loop=1:length(SET)
      dosplit=0;
      margin = fontsize;
      str = strtrim(SET(loop).ImageType);
      
      listmap={'Feature tracking','FT';...
        'Late enhancement','LGE';...
        'Perfusion rest', 'Perf R';...
        'Perfusion rest aligned', 'Perf R a';...
        'Perfusion stress', 'Perf S';...
        'Perfusion stress aligned', 'Perf S a.';...
        'Strain from tagging','Tagging';...
        'T1 map pre','T1 Pre';...
        'T1 map post', 'T1 Post';...
        'T2 map pre','T2 Pre';...
        'T2 map post', 'T2 Post';...
        'T2star', 'T2*';...
        'Flow (magnitude)','Flow';...
        'Flow (through plane)','Flow'};
      
      for i = 1:size(listmap,1)
        if strcmpi(listmap{i,1},str)
          str=listmap{i,2};
        end
      end
      
      imvpmap={'Short-axis', 'SAX';...
        'Short-axis basal', 'SAX basal';...
        'Short-axis', 'SAX mid';...
        'Short-axis', 'SAX apical';...
        'Sagittal','Sag';...
        'Coronal','Cor';...
        'Frontal','Front';...
        'Transversal', 'Trans'
        '2CH','2CH';...
        '3CH','3CH';...
        '4CH','4CH'};
      
      addimvp=0;
      if isempty(SET(loop).ImageViewPlane)
          SET(loop).ImageViewPlane = '';
      end;
      strimvp = strtrim(SET(loop).ImageViewPlane);
      for i = 1:size(imvpmap,1)
        if strcmpi(imvpmap{i,1},strimvp)
          strimvp=imvpmap{i,2};
          addimvp=1;
        end
      end
      
      if addimvp
        str=[str,', ',strimvp];
      end
      
      if length(str)>10
        str=str(1:8);
      end
      
%       if length(str)>numletters
%         %group into cell of strings
%         %check for spaces close to 10 line
%         formatstring='%d. %s';
%         strcell={};
%         str_tmp=str;
%         numbreaks=0;
%         
%         while true
%           if numbreaks>0
%             numletters=12;
%           end
%           
%           ind = find(isspace(strtrim(str_tmp)));
%           tmp=ind(find((ind>floor(1/2*numletters)).*(ind<numletters)));
%           if length(str_tmp)<=numletters
%             strcell = [strcell,str_tmp];
%             break;
%           elseif ~isempty(tmp)
%             formatstring=[formatstring,'\n%s'];
%             strcell = [strcell,str_tmp(1:tmp(end))];
%             str_tmp = str_tmp(tmp(end)+1:end);
%           else
%             formatstring=[formatstring,'\n%s'];
%             strcell = [strcell,str_tmp(1:numletters)];
%             str_tmp=str_tmp(numletters+1:end);
%           end
%           numbreaks = numbreaks + 1;
%           end
%         margin=(numbreaks+1)*fontsize;
%         dosplit=1;
%       end

      ypos = (loop-1)*thumbsize+1+2*margin;
      xpos = margin;
      
%       if dosplit
%       g.Handles.datasetnumber =  ...
%         [g.Handles.datasetnumber  ...
%         text(xpos,ypos,sprintf(formatstring,loop,strcell{:}),...
%         'parent',g.Handles.datasetaxes,...
%         'color',g.GUISettings.ThumbFlowLineColor,...
%         'fontsize',fontsize)...
%         ];
%       else
      g.Handles.datasetnumber =  ...
        [g.Handles.datasetnumber  ...
        text(xpos,ypos,sprintf('%d. %s',loop, str),...
        'parent',g.Handles.datasetaxes,...
        'color',g.GUISettings.ThumbFlowLineColor,...
        'fontsize',fontsize)...
        ];
%       end
    end;
    
    end
    
    %----------------------------------------
    function pathname = getpreferencespath(g)
    %----------------------------------------
    %Get path to preferences folder
    foldername = g.ProgramFolderName;
    
    pathname = '';
    customfile = fullfile(pwd,'preferencespath.txt');
    if exist(customfile,'file')
      try
        fid = fopen(customfile);
        pathname = fgetl(fid);
        fclose(fid);
      catch %#ok<CTCH>
      end
    end
    
    if ~isempty(pathname) && exist(pathname,'dir')
      return
    end
    
    if ispc
      %Is Windows, check for better
        temp = getenv('APPDATA');
        if not(isempty(temp))
          pathname = temp;
        else
          temp = getenv('USERPROFILE');
          if not(isempty(temp))
            pathname = temp;
          else
            temp = getenv('HOMEPATH');
            if not(isempty(temp))
              pathname = temp;
            end;
          end;
        end;
      
      %Check if subdirectory exists
      if not(exist([temp filesep foldername],'dir'))
        temppath = pwd;
        cd(pathname);
        disp(dprintf('Creating new folder %s%s%s',temp,filesep,foldername));
        suc = mkdir(foldername);
        cd(temppath);
        if not(suc)
          myfailed(dprintf('Could not create %s%s%s',temp,filesep,foldername));
        end;
      end;
      
      %Add Segment to the path
      pathname= [pathname filesep foldername];
    else
      pathname = pwd;
    end; %is pc
    end
    
    %-----------------------
    function standardpref(g)
    %-----------------------
    %Setup language. Overloaded in GUI's where translation is available.
    %Can also be used to set other preferences that are not optional.
    g.Pref.Language = 'English'; %Override language selection from other GUI
    end
    
    
    %----------------------------
    function updatepacslabels(g)
    %----------------------------
    %Change to name of PACS in labels if only one PACS connection exists.
    %Overloaded in Segment CMR GUI.
    end
    
    %------------------------------
    function updatesax3display(g,~)
    %------------------------------
    %Update data used to correctly display segmentation in SAX3 view.
    %Overloaded in Segment SPECT GUI
    end
    
    %------------------------
    function setglobalstacks(g)
    %------------------------
      %Find the best stacks for report LV, RV and Flow in main result
      
      global DATA
      
      cinestacks = findfunctions('findcineshortaxisno',true);
      if ~isempty(cinestacks)
        DATA.LVNO = cinestacks(1);
        DATA.RVNO = cinestacks(2);
      else
        DATA.LVNO = [];
        DATA.RVNO = [];
      end
      [DATA.FlowNO, DATA.FlowROI] = findfunctions('findflowaxisno');
    end
    
    
    %-----------------------
    function initbullseyeslices(g) %#ok<MANU>
    %-----------------------
    %Overloaded in Segment CMR GUI, Segment CT GUI and Segment Spect GUI
    end
    
    %-----------------------
    function initbullseyes(g) %#ok<MANU>
    %-----------------------
    %Overloaded in Segment CT GUI
    end
    
    %--------------------------
    function r = isorthoview(g)
    %-------------------------- 
    %Check if we have a proper setup for orthoview
    r = numel(g.ViewPanels) >= 3 && ...
      all(g.ViewPanels(1:3) == g.ViewPanels(1));
    end
    
%     %-------------------------------
%     function keyreleased(g)
%     %------------------------------- 
%     g.GUI.needrest = 1;
%     end
    %----------------------
    function synchronize(g)
    %-----------------------
    %toggle synchronize option currently makes 
    %slice toggling 
    %time toggling
    %go to ed
    %go to es
    %execute for all stacks.
    
    stateandicon = segment('iconson','synchronize');
    g.Synchronize=stateandicon{1};
    end
    %---------------------------------
    function keypressed(g,fignum,evnt)
    %---------------------------------
    %Called when key is pressed. Overloaded in Segment CT GUI.
    global SET NO DATA
    persistent buffer;
    persistent starttime;  
    tol=1e-6;
    diff=rem(now,1)-starttime;
    if ~isempty(diff) && diff<tol
      starttime=[];
      return
    else
      starttime = rem(now,1);
    end
    
    switch nargin
      case 1
        disp('No inputarguments.');
        return;
      case 2
        if g.Testing
          key = fignum;
        end;
      case 3
        key = getkey(evnt);
    end;
    
    if g.RecordMacro
      macro_helper('keypressed',['''' key '''']);
    end;
    
    isfkey = regexp(key,'[f]\d','once');
    if isfkey
      segmentview('keypressed',key);
      return
    end
        
    %Use same order as in manual
    switch key
      case 'scrolllock'
        if ~g.Synchronize
          indent(g.Handles.hideiconholder,'synchronize',1)
        else
          undent(g.Handles.hideiconholder,'synchronize',1)
        end
        %Tools for changing image frame or slice
      case 'd'
        if ~isempty(g.GUI.T2Star)
          %Call t2star centerpoint move funtion:
          fcn = 't2star.t2star(''point_right'')';
          eval(fcn);
        else
          %Diastole (d)
          if DATA.Synchronize
            tools('enddiastoleall_Callback');
          else
            tools('enddiastole_Callback');
          end
        end
      case 'n'
        if DATA.CurrentPanel-1~=0
          segment('switchtopanel',DATA.CurrentPanel-1)
        else
          segment('switchtopanel',numel(DATA.ViewPanels))
        end
      case 'm'
        if DATA.CurrentPanel+1<=numel(DATA.ViewPanels)
          segment('switchtopanel',DATA.CurrentPanel+1);
        else
          segment('switchtopanel',1);
        end
      case 's'
        if ~isempty(g.GUI.T2Star)
          %Call t2star centerpoint move funtion:
          fcn = 't2star.t2star(''point_down'')';
          eval(fcn);
        else
          %Systole (s)
          if DATA.Synchronize
           tools('endsystoleall_Callback');
          else
           tools('endsystole_Callback');
          end
        end
      case 'shift-d'
        tools('enddiastoleall_Callback');
      case 'shift-s'
        tools('endsystoleall_Callback');
      case 'leftarrow' %Left
        if DATA.Synchronize
          segment_main('previousallframe_Callback');
        else
          segment_main('previousframe_Callback');
        end
        
      case 'rightarrow' %Right
        if DATA.Synchronize
          segment_main('nextallframe_Callback');
        else
          segment_main('nextframe_Callback');
        end
       case 'uparrow' %Arrow up
       if DATA.Synchronize
         segment_main('movealltowardsbase_Callback');
       else
         segment_main('movetowardsbase_Callback');
       end
      case 'downarrow' %Arrow down
        if DATA.Synchronize
          segment_main('movealltowardsapex_Callback');
        else
          segment_main('movetowardsapex_Callback');
        end
        case 'shift-leftarrow'
        segment_main('previousallframe_Callback');
      case 'shift-rightarrow'
        segment_main('nextallframe_Callback');
      case 'shift-uparrow'
        segment_main('movealltowardsbase_Callback');
      case 'shift-downarrow'
        segment_main('movealltowardsapex_Callback');
      case 'alt-uparrow'
        segment_main('smoothendowall_Callback');
      case 'c'
        segment_main('cinetool_Callback');
      case 'p'
        %Play p
        stateandicon=segment('iconson','play');
        stateandicon{2}.isindented=~stateandicon{2}.isindented;
        if stateandicon{2}.isindented
          stateandicon{2}.cdataDisplay=stateandicon{2}.cdataIndent;
        else
          stateandicon{2}.cdataDisplay=stateandicon{2}.cdataIndent;
        end
        DATA.Handles.permanenticonholder.render;
        segment('playall_Callback');%,'keypressed');
        %segment_main('playmovie_Callback','keypressed');
      case 'shift-p'
        %Play all p
        segment_main('playmovie_Callback','keypressed');
        %Tools for viewing and selecting
      case 'l'
        iconcell=DATA.Handles.toggleiconholder.iconCell;
        for i=1:numel(iconcell)
          iconcell{i}.undent
          if strcmpi(iconcell{i}.name,'ribbonlv')
            modenbr = i;
          end
        end
         iconcell{modenbr}.cdataDisplay=iconcell{modenbr}.cdataIndent;
         iconcell{modenbr}.isindented=1;
        DATA.Handles.toggleiconholder.render;
        DATA.togglebuttonLV_Callback;
      case 'r'
        iconcell=DATA.Handles.toggleiconholder.iconCell;
        for i=1:numel(iconcell)
          iconcell{i}.undent
          if strcmpi(iconcell{i}.name,'ribbonrv')
            modenbr = i;
          end
        end
         iconcell{2}.cdataDisplay=iconcell{2}.cdataIndent;
         iconcell{2}.isindented=1;
        DATA.Handles.toggleiconholder.render;
        DATA.togglebuttonRV_Callback;
        %segment_main('viewrefresh_Callback');
      case 'f'
        iconcell=DATA.Handles.toggleiconholder.iconCell;
        for i=1:numel(iconcell)
          iconcell{i}.undent
           if strcmpi(iconcell{i}.name,'ribbonflow')
            modenbr = i;
          end
        end
         iconcell{modenbr}.cdataDisplay=iconcell{modenbr}.cdataIndent;
         iconcell{modenbr}.isindented=1;
        DATA.Handles.toggleiconholder.render;
        DATA.togglebuttonROIFLOW_Callback;
      case '3'
        iconcell=DATA.Handles.toggleiconholder.iconCell;
        for i=1:numel(iconcell)
          iconcell{i}.undent
          if strcmpi(iconcell{i}.name,'ribbon3dp')
            modenbr = i;
          end
        end
         iconcell{modenbr}.cdataDisplay=iconcell{modenbr}.cdataIndent;
         iconcell{modenbr}.isindented=1;
        DATA.Handles.toggleiconholder.render;
        DATA.togglebutton3dp_Callback;
      
      case 'u'
        
        if isfield(SET(NO).StrainTagging,'runningregistration') && SET(NO).StrainTagging.runningregistration
          mywarning('Image enhancement not available during registration.')
          return
        end
        
        if isempty(SET(NO).IntensityMapping.Compression)
        SET(NO).IntensityMapping.Compression=1;  
        end
        SET(NO).IntensityMapping.Compression=SET(NO).IntensityMapping.Compression+0.1;
   
        segment('makeviewim',DATA.CurrentPanel,NO)
        drawfunctions('drawimageno')

      case 'j'
        if isfield(SET(NO).StrainTagging,'runningregistration') && SET(NO).StrainTagging.runningregistration
          mywarning('Image enhancement not available during registration.')
          return
        end
        
        if isempty(SET(NO).IntensityMapping.Compression)
        SET(NO).IntensityMapping.Compression=1;  
        end
        SET(NO).IntensityMapping.Compression=SET(NO).IntensityMapping.Compression-0.1;
        segment('makeviewim',DATA.CurrentPanel,NO)

        drawfunctions('drawimageno')
      case 'w'
        iconcell=DATA.Handles.toggleiconholder.iconCell;
        for i=1:numel(iconcell)
          iconcell{i}.undent
          if strcmpi(iconcell{i}.name,'ribbonviability')
            modenbr = i;
          end
        end
         iconcell{modenbr}.cdataDisplay=iconcell{modenbr}.cdataIndent;
         iconcell{modenbr}.isindented=1;
        DATA.Handles.toggleiconholder.render;
        DATA.togglebuttonVia_Callback;
      case 'a'
          iconcell=DATA.Handles.toggleiconholder.iconCell;
        for i=1:numel(iconcell)
          iconcell{i}.undent
          if strcmpi(iconcell{i}.name,'ribbonanalysis')
            modenbr = i;
          end
        end
        iconcell{modenbr}.cdataDisplay=iconcell{modenbr}.cdataIndent;
        iconcell{modenbr}.isindented=1;
        DATA.Handles.toggleiconholder.render;
        DATA.togglebuttonAnalysis_Callback;
      case 'i'
        iconcell=DATA.Handles.toggleiconholder.iconCell;
        for i=1:numel(iconcell)
          iconcell{i}.undent
          if strcmpi(iconcell{i}.name,'ribbonimage')
            modenbr = i;
          end
        end
        
        iconcell{modenbr}.cdataDisplay=iconcell{modenbr}.cdataIndent;
        iconcell{modenbr}.isindented=1;
        DATA.Handles.toggleiconholder.render;
        DATA.togglebuttonImage_Callback;  
      case 't'
        iconcell=DATA.Handles.toggleiconholder.iconCell;
        
        for i=1:numel(iconcell)
          iconcell{i}.undent
          if strcmpi(iconcell{i}.name,'ribbontxmap')
            modenbr = i;
          end
        end
        
        iconcell{modenbr}.cdataDisplay=iconcell{modenbr}.cdataIndent;
        iconcell{modenbr}.isindented=1;
        DATA.Handles.toggleiconholder.render;
        DATA.togglebuttontxmap_Callback;
        
      case 'h'
        stateandicon=segment('iconson','hideall');
        if ~stateandicon{1}
          stateandicon{2}.cdataDisplay=stateandicon{2}.cdataIndent;
          stateandicon{2}.isindented=1;
        else
          stateandicon{2}.cdataDisplay=stateandicon{2}.cdata;
          stateandicon{2}.isindented=0;
        end
        segment_main('viewhideall_Callback')
        
        g.Handles.permanenticonholder.render;
%         h = [...
%           g.Handles.hidepinsicon ...
%           g.Handles.hideothercontouricon ...
%           g.Handles.hideinterpicon ...
%           g.Handles.hidelvicon ...
%           g.Handles.hidervicon ...
%           g.Handles.hidescaricon ...
%           g.Handles.hidemaricon ...
%           g.Handles.hideroiicon ...
%           g.Handles.hidemeasuresicon ...
%           g.Handles.hidepointsicon ...
%           g.Handles.hideplusicon ...
%           g.Handles.hideintersectionsicon ...
%           g.Handles.hidetexticon ...
%           g.Handles.hidepapicon ...
%           g.Handles.hideoverlayicon];
%         viewhidestatus = get(h,'state');
%         if all(strcmp(viewhidestatus,'on'))
%           set(h,'state','off');
%         else
%           set(h,'state','on');
%         end
%         drawfunctions('drawall');
      case 'v'
        %View mode (v)
        switch g.ViewPanelsType{g.CurrentPanel}
          case 'one'
            segment_main('viewimage_Callback','montage');
          case 'mmode'
            segment_main('viewimage_Callback','one');
          case 'montage'
            segment_main('viewimage_Callback','one');
          case 'montagerow'
            segment_main('viewimage_Callback','one');
          case 'montagefit'
            segment_main('viewimage_Callback','one');
          case 'sax3'
            segment_main('viewimage_Callbac','one');
        end;
      case 'y'
        %select select tool
        g.Handles.configiconholder.indent('select',1);
      case 'ctrl-a'
        segment_main('selectallslices_Callback');
      case 'shift-u' %unselect all image stacks
        segment_main('unselectallslices_Callback');
      case 'shift-a' %shift-a
        segment_main('viewallimagestacks_Callback');
      case 'shift-1'
        drawfunctions('drawall',1);
      case 'shift-2'
        drawfunctions('drawall',2);
      case 'alt-2'
        drawfunctions('drawall',2,1);
      case 'shift-3'
        drawfunctions('drawall',3);
      case 'alt-3'
        drawfunctions('drawall',3,1);
      case {'shift-4','alt-4'}
        drawfunctions('drawall',4);
      case 'shift-6'
        drawfunctions('drawall',6);
      case 'alt-6'
        drawfunctions('drawall',6,1);
      case 'shift-9'
        drawfunctions('drawall',9);
        % These are now set through GUIDE, linked by 'accelerator' for the menu
        % items (note that numbers can't normally be set, in the Menu editor you
        % have to select More Properties to set 1,2,3,4...
        % This so that they can be disabled by handles.
      case 'ctrl-1'
        segment_main('viewimage_Callback','one');
      case 'ctrl-2'
        segment_main('viewimage_Callback','mmode');
      case 'ctrl-3'
        segment_main('viewimage_Callback','montage');
      case 'ctrl-4'
        segment_main('viewimage_Callback','montagerow');
      case 'ctrl-5'
        segment_main('viewimage_Callback','montagefit');
        
        % Tools for segmentation
      case 'ctrl-l'
        lvsegmentation;%lv('segmentfullyautomatic_Callback');
      case 'ctrl-m'
        lvpeter('segmentendo_Callback');
      case 'shift-ctrl-m'
        lvpeter('segmentepi_Callback',false);
      case 'ctrl-alt-m'
        rv('segmentrvendo_Callback');
      case 'shift-alt-m'
        return; %reserved for segment RV epi
      case 'ctrl-r'
        lvpeter('segmentrefineendo_Callback',false,false);
      case 'shift-ctrl-r'
        lvpeter('segmentrefineepi_Callback',false);
      case 'ctrl-alt-r'
        rv('segmentrefinervendo_Callback',false,false);
      case 'alt-r'
        flow('flowrefine_Callback');
      case 'ctrl-f'
        lvpeter('segmentpropagateendo_Callback');
        %     lvpeter('segmentpropagateendoorepi_Callback');%Endo or epi
      case 'shift-ctrl-f'
        lvpeter('segmentpropagateepi_Callback');
        %      lv('segmentpropagateendoandepi_Callback');%Endo and epi
      case 'alt-f'
        flow('flowpropagate_Callback');%ROI
      case 'ctrl-alt-f'
        segment('copyrvendoforward');
      case 'shift-alt-f'
        return;%reserved for propagate RV epi
      case 'alt-t'
        flow('flowtrackroi_Callback');
      case 'ctrl-u'
        %tools('copyupward_Callback','endo',false,true); %LV endo %silent,dolv
        tools('copyupward_Callback','lv',false,true); %Copy upward both endo & Epi and refine
        %   case 'shift-ctrl-u'
        %     tools('copyupward_Callback','epi',false,true); %LV epi %silent,dolv
      case 'ctrl-alt-u'
        tools('copyupward_Callback','endo',false,false); %RV endo %silent,dolv
      case 'shift-alt-u'
        return; %reserved for copy upwards  for RV epi
%       case 'alt-u'
%         return;%reserved for copy upwards of ROI
%       case 'alt-d'
%         return;%reserved for copy downwards of ROI
      case 'ctrl-d'
        %tools('copydownward_Callback','endo',false,true); %LV endo %silent,dolv
        tools('copydownward_Callback','lv',false,true); %Copy downward both endo & epi
        %   case 'shift-ctrl-d'
        %     tools('copydownward_Callback','epi',false,true); %LV epi %silent,dolv
      case 'ctrl-alt-d'
        tools('copydownward_Callback','endo',false,false); %RV endo %silent,dolv
      case 'ctrl-alt-w'
        viability('weightedsliderautopushbutton_Callback'); %Auto position weighted slider, internal use only
      case 'shift-alt-d'
        return; %reserved for copy downwards  for RV epi
      case 'ctrl-v'
        lvpeter('segmentremovepapilars_Callback');
      case 'shift-alt-r'
        lv('segmentrefineendo_Callback',false,false);
      case 'ctrl-e' %expand endo
        lv('segmentexpandcontract_Callback',1,'endo');
      case 'ctrl-alt-e' %expand epi
        lv('segmentexpandcontract_Callback',1,'epi');
      case 'ctrl-k' %contract endo
        lv('segmentexpandcontract_Callback',-1,'endo');
      case 'ctrl-alt-k' %contract epi
        lv('segmentexpandcontract_Callback',-1,'epi');
      case 'o' %smooth latest segmentation
        tools('smoothsegmentation_Callback');
        
      %tools for define ED and ES
      case 'alt-d'
        tools('seted_Callback')
      case 'alt-s'
        tools('setes_Callback')
      case 'alt-h'
        val=get(DATA.Handles.hideallpanelscheckbox,'Value');
        set(DATA.Handles.hideallpanelscheckbox,'Value',~val);
        segment('hideallpanels');        
      %Tools for translating contour
      case 'shift-alt-a'
        tools('translatecontoursandimage',0,-1); %left
      case 'shift-alt-x'
        tools('translatecontoursandimage',0,1); %right
      case 'shift-alt-w'
        tools('translatecontoursandimage',-1,0); %up
      case 'shift-alt-z'
        tools('translatecontoursandimage',1,0); %down
      case 'alt-a'
        tools('translatecontours',0,-1); %left
      case 'alt-x'
        tools('translatecontours',0,1); %right
      case 'alt-w'
        tools('translatecontours',-1,0); %up
      case 'alt-z'
        tools('translatecontours',1,0); %down
        
      case 'space' %space
        %set(g.fig,'CurrentObject',g.Handles.volumeaxes);
        %set(0,'CallbackObject',g.Handles.volumeaxes);
        %get(get(g.fig,'CurrentObject'),'tag')
%           'drawendo', 'drawepi', ...
%           'drawepi', 'drawendo', ...
%           'drawrvendo', 'drawrvepi', ...
%           'drawrvepi', 'drawrvendo', ...

        toolstruct = struct(...
          'epipin', 'endopin', ...
          'endopin', 'epipin', ...
          'putroi', 'scale', ...
          'move', 'scale', ...
          'scale', 'move', ...
          'drawendo', 'drawepi', ...
          'drawepi', 'drawendo', ...
          'drawrvendo', 'drawrvepi', ...
          'drawrvepi', 'drawrvendo', ...
          'interpendo', 'interpepi', ...
          'interpepi', 'interpendo', ...
          'interprvendo', 'interprvepi', ...
          'interprvepi', 'interprvendo', ...
          'drawscar', 'drawrubberpen', ...
          'drawrubberpen', 'drawscar', ...
          'drawroi', 'select', ...
          'drawmarpen', 'drawmarrubberpen', ...
          'drawmarrubberpen', 'drawmarpen', ...
          'select', 'drawscar');
        
        %if isequal(g.CurrentTool, 'drawendo')
          drawnow;
        %end
        if isequal(g.CurrentTool, 'select') && isempty(SET(NO).Scar)
          %Loop over stacks to find next
          for loop=1:(length(g.ViewPanels)+1)
            
            %increase
            temp = g.CurrentPanel+1;
            if temp>length(g.ViewPanels)
              temp = 1;
            end;
            g.CurrentPanel = temp;
            
            if g.ViewPanels(temp)~=0
              %valid
              segment_main('switchtopanel',temp);
              g.CurrentPanel = temp;
              return;
            end;
          end;
          segment_main('switchtopanel',g.CurrentPanel)
        else
          %Toggle the gui
          %Given an update tool string toggle icons 2 the given tool
%           if strcmp(g.CurrentTool(1:4),'draw');
%             g.CurrentTool(1:4)=[];
%             g.CurrentTool=[g.CurrentTool,'pen'];
%           end
%
%      do buttonup before in case it didn't manage.
          %set(DATA.imagefig,'WindowButtonMotionFcn','')
%           buttonupfcn = get(DATA.imagefig,'WindowButtonUpFcn');
%           set(DATA.imagefig,'WindowButtonUpFcn','')
%           run(buttonupfcn);
%           
% %          this is safe since a buttonup will be toggled to!
%           set(DATA.imagefig,'WindowButtonDownFcn','')
%           
%           set(DATA.imagefig,'WindowButtonMotionFcn','')
%           set(DATA.imagefig,'WindowButtonUpFcn','')
          oldtool=g.CurrentTool;
          g.CurrentTool=toolstruct.(g.CurrentTool);
          g.guitoggletool(g.CurrentTool,oldtool);
%           donotchangetfo = true;
%           updatetool(toolstruct.(g.CurrentTool), [], donotchangetfo);
        end
        
%       case 'shift-l' 
%        % g.updateicons('lv'); %themes: L(LV),R(RV),V(Scar/Viability),F(ROI/Flow), M(MaR), I(Misc), E(perfusion), T(T2*)
%       case 'shift-r'
%         %g.updateicons('rv'); %themes: L,R,V,F
%       case 'shift-f'
%         %g.updateicons('roi'); %themes: L(LV),R(RV),V(Scar/Viability),F(ROI/Flow), M(MaR), I(Misc)
%       case 'shift-v'
%         %g.updateicons('scar'); %themes: L,R,V,F
%       case 'shift-m'
%         %g.updateicons('mar'); %themes: L(LV),R(RV),V(Scar/Viability),F(ROI/Flow), M(MaR), I(Misc), E(perfusion), T(T2*)
%       case 'shift-i'
%         %g.updateicons('misc'); %themes: L(LV),R(RV),V(Scar/Viability),F(ROI/Flow), M(MaR), I(Misc)
%       case 'shift-e'
%         return; %reserved for theme Perfusion
%       case 'shift-t'
%        return;% reserved for g.updateicons('t2'); %themes: L(LV),R(RV),V(Scar/Viability),F(ROI/Flow), M(MaR), I(Misc), E(perfusion), T(T2*)
%       case 'shift-n'  %Select LV Endo pen
%         %g.updateicons('lv');
%         updatetool('drawendo');
%       case 'shift-b'  %Select LV Epi pen
%         %g.updateicons('lv');
%         updatetool('drawepi');
%       case 'shift-g'  %select LV Endo interp
%         %g.updateicons('lv');
%         updatetool('interpendo');        
%       case 'shift-h'  %Select LV Epi interp
%         %g.updateicons('lv');
%         updatetool('interpepi');       
%         

      case 'e'
        %compress imagemapping
        
        %Miscellaneous commands
      case 'q'
        fcn = get(g.Handles.imagehandle,'ButtonDownFcn');
        set(g.imagefig,'SelectionType','normal')
        eval(fcn);
        fcn = get(g.imagefig,'WindowButtonUpFcn');
        pause(0.3);
        eval(fcn);
      case 'w'
        if ~isempty(g.GUI.T2Star)
          %Call t2star centerpoint move funtion:
          fcn = 't2star.t2star(''point_up'')';
          eval(fcn);
        else
          fcn = get(g.Handles.imagehandle,'ButtonDownFcn');
          set(g.imagefig,'SelectionType','alt')
          eval(fcn);
        end
      case 'ctrl-b'
        reportbullseye;
      case 'ctrl-n'
        filemenu('loadnext_Callback');
        %   case 'ctrl-o' open file loader set in GUI to be able to overload callback depending on gui
      case 'ctrl-p'
        patientdatabase;
      case 'shift-ctrl-w'
        g.filecloseall_Callback;
        %   case 'ctrl-q' quit program set in GUI to be able to overload callback depending on gui
      case 'ctrl-z'
        tools('undosegmentation_Callback');
        % Set in GUIDE, so they are disabled, but leave here for documentation!
      case 'ctrl-0'
        segment_main('viewzoomin_Callback');
      case 'ctrl-hyphen'
        segment_main('viewzoomout_Callback');        
      case 'shift-ctrl-b'
        segmenthelp('supportreq_Callback'); %Open supportreq GUI
      case 'ctrl-c'
        return;
      case 'shift-ctrl-e'
        segment_main('evalcommand_Callback');
      case 'ctrl-t' %plot flow
        reportflow();
        
      case 'k' % TEMPORARY HOTKEY, MAY BE CHANGED
        plugin_3dflow('calckineticenergy_Callback');
      case 'shift-k' % TEMPORARY HOTKEY, MAY BE CHANGED
        plugin_3dflow('createvelocityset_Callback');
      case 'shift-ctrl-k' % TEMPORARY HOTKEY, MAY BE CHANGED
        reportflow3d2d(NO);
      case 'y' %TEMPORARY
        v = get(g.Handles.imageaxes(g.CurrentPanel),'view');
        set(g.Handles.imageaxes(g.CurrentPanel),'view',[v(1)+1 90]);
      case 'shift-y' %TEMPORARY
        v = get(g.Handles.imageaxes(g.CurrentPanel),'view');
        set(g.Handles.imageaxes(g.CurrentPanel),'view',[v(1)-1 90]);
      case 'shift-ctrl-t' %TEMPORARY
        t2star.t2starnew('calc');
      case 'ctrl-alt-t'
        g.maketest;
      case 'ctrl-alt-a'
        %Fully Automatic Vessel Segmentation in 2D PC-MRI:
        vesselsnake_overview('segmentaorta')
      case 'ctrl-alt-v'
        %New method for semi-Automatic Vessel Segmentation in 2D PC-MRI:
        oldpref=DATA.ThisFrameOnly;
        vesselsnake_flowtrackroi('flowtrackroi');
        DATA.ThisFrameOnly=oldpref;
      case 'ctrl-alt-space'
        %Import T1map:
        cxt1map2segment()
      case 'ctrl-alt-z'
        spect.spectautobutton()
%       case 'tab'
%         ind=DATA.Handles.toggleiconholder.findindented;
%         indnext=ind+1;
%         iconcell = DATA.Handles.toggleiconholder.iconCell;
%         if indnext>DATA.Handles.toggleiconholder.numberoficons;
%           indnext=1;
%         end
%           iconcell{ind}.cdataDisplay=iconcell{ind}.cdata;
%           iconcell{ind}.isindented=0;
%           iconcell{indnext}.cdataDisplay=iconcell{indnext}.cdataIndent;
%           iconcell{indnext}.isindented=1;
% %           switch indnext
% %             case 1
% %               g.togglebuttonLV_Callback
% %             case 2
% %               g.togglebuttonRV_Callback
% %             case 3
% %               g.togglebuttonROIFLOW_Callback
% %             case 4
% %               g.togglebuttonVia_Callback
% %             case 5
% %               g.togglebuttonAnalysis_Callback
% %             case 6
% %               g.togglebuttonImage_Callback
% %           end
% %        run(iconcell{ind}.execute);
%           feval(iconcell{indnext}.execute);
%           g.Handles.toggleiconholder.render;
% %          g.Handles.configiconholder.render;

          %Set keypress functions
          %set(DATA.imagefig,'keypressfcn',@(fignum,evnt)segment('keypressed',fignum,evnt));
          %updatetool(g.CurrentTool)
%          gcf(DATA.fig)
    end;
    end
    
    %------------------------------------------------
    function guitoggletool(varargin)
    %------------------------------------------------
    g=varargin{1};
    tool2toggle2=varargin{2};
    oldtool=varargin{3};
    %Given an update tool string toggle icons 2 the given tool 
    if strcmp(tool2toggle2(1:4),'draw');
      tool2toggle2(1:4)=[];
      tool2toggle2=[tool2toggle2,'pen'];
    end
    
    if strcmp(oldtool(1:4),'draw');
      oldtool(1:4)=[];
      oldtool=[oldtool,'pen'];
    end
    
    g.Handles.configiconholder.indent(tool2toggle2,1)
    
%     icons={g.Icons.lviconcell;g.Icons.rviconcell;g.Icons.roiflowiconcell;g.Icons.viabilityiconcell;...
%       g.Icons.analysisiconcell; g.Icons.imageiconcell};
%     n=length(icons);
%     for k=1:n;
%       iconcell=icons{k};
%       N=length(iconcell);
%       %return icon state
%       for i=1:N
%         if strcmp(iconcell{i}.name,tool2toggle2)
%            holderind=k;
%            iconind=i;
%         end
%         %if strcmp(iconcell{i}.name,oldtool)
%          %  oldholderind=k;
%         %end
%       end
%     end
%     icon2indent=[];
%     icon2undent=[];
%     iconcell=g.Handles.configiconholder.iconCell;
%     N=length(iconcell);
%     for i=1:N
%       if strcmp(iconcell{i}.name,tool2toggle2)
%       icon2indent=iconcell{i};
%       end
%       if strcmp(iconcell{i}.name,oldtool)
%         icon2undent=iconcell{i};
%       end
%     end
    
   % if holderind==oldholderind
      %undent current tool this is safe since tools always are toggle
      %buttons.
      %stateandicon=segment('iconson',oldtool);
      %stateandicon{2}.isindented=0;
      %stateandicon{2}.cdataDisplay=stateandicon{2}.cdata;
    %end
    
    %g.Handles.toggleiconholder.clickedicon=icon2indent;%g.Handles.toggleiconholder.iconCell{holderind};
      %we want to go to the toggle menu with the tool2toggle2
      
%       
%       if ~isempty(icon2indent)
%         icon2indent.cdataDisplay=icon2indent.cdataIndent;%icons{holderind}{iconind}.cdataIndent;
%         icon2indent.isindented=1;
%       end
%       g.Handles.configiconholder.clickedicon=icon2indent;%{holderind}{iconind};
%       if ~isempty(icon2undent)
%         icon2undent.isindented=0;
%         icon2undent.cdataDisplay=icon2undent.cdata;
%       end
%     g.Handles.configiconholder.render;
    end
    
    %------------------------------------------------------
    function [measure,slice] = getmeasurecoords(g,no,panel)
    %------------------------------------------------------
    %Get coordinates of measurements with respect to the current view
    global SET
    
    switch g.ViewPanelsType{panel}
      case 'hla'
        measure = struct( ...
          'X',{SET(no).Measure.Z}, ...
          'Y',{SET(no).Measure.Y}, ...
          'Z',{SET(no).Measure.X});
        slice = SET(no).HLA.slice;
      case 'vla'
        measure = struct( ...
          'X',{SET(no).Measure.Z}, ...
          'Y',{SET(no).Measure.X}, ...
          'Z',{SET(no).Measure.Y});
        slice = SET(no).VLA.slice;
      case 'gla'
        ang = SET(no).GLA.angle;
        res = SET(no).ResolutionY*cos(ang)+SET(no).ResolutionX*abs(sin(ang));
        measure = struct( ...
          'X',{SET(no).Measure.Z}, ...
          'Y',cellfun(@(y,x)1/res*(SET(no).ResolutionY*(y-SET(no).GLA.y0)*cos(ang) + ...
          SET(no).ResolutionX*(x-SET(no).GLA.x0)*sin(ang)), ...
          {SET(no).Measure.Y},{SET(no).Measure.X},'UniformOutput',false), ...
          'Z',cellfun(@(y,x)1/res*(SET(no).ResolutionY*(y-SET(no).GLA.y0)*-sin(ang) + ...
          SET(no).ResolutionX*(x-SET(no).GLA.x0)*cos(ang)), ...
          {SET(no).Measure.Y},{SET(no).Measure.X},'UniformOutput',false));
        slice = 0;
      otherwise
        measure = struct( ...
          'X',{SET(no).Measure.X}, ...
          'Y',{SET(no).Measure.Y}, ...
          'Z',{SET(no).Measure.Z});
        slice = SET(no).CurrentSlice;
    end
    end
    
  end
  
  methods(Access='protected')
    
    %---------------------------------------
    function switchtoimagestack_part2(g, no)
    %---------------------------------------
    %switchtoimagestack continued. Overloaded in SegmentGUI and CVQgui
    global SET NO

    %Loop until there are no blocking windows left
    while ~isempty(g.BlockingFigs)

      %Check deleted windows
      ind = true(1,length(g.BlockingFigs));
      for loop=1:length(ind)
        try
          get(g.BlockingFigs);
        catch %#ok<CTCH>
          ind(loop) = false;
        end;
      end;

      %Remove deleted windows
      g.BlockingFigs = g.BlockingFigs(ind);

      if ~isempty(g.BlockingFigs)
        uiwait(msgbox('Close all reviewing windows, then press OK before changing current image stack.'));
      end;
    end;

    %if switchtostack is not in view, update current panel
     if ~ismember(no,g.ViewPanels)
       emptypanel = find(g.ViewPanels==0,1);
       if isempty(emptypanel)
         emptypanel = numel(g.ViewPanels);
       end
       g.ViewPanels(emptypanel) = no;
       g.ViewPanelsType{emptypanel} = 'one';
       g.ViewIM{emptypanel} = [];
       drawfunctions('drawimagepanel',emptypanel);
       segment('switchtopanel',emptypanel,false); %Do not update imagestack
     end
    
    %Store it
    NO = no;
    
    %Double check out of range
    if NO>length(SET)
      NO = length(SET);
    end;

    if NO<1
      NO = 1;
    end;

    %Variables that has to be reset in the change are:
    % endoedgedetect, epiedgedetect, endoedge, epiedge
    % BpInt, MInt, BALLOON, Montage, Wallthickness, Radius,
    % MaxRadialVel, BalloonLevel, run
    % interactionlock, LevelSet

    %Update viability menu if no scar is present then
    %options are hidden.
%    viability('viabilitymenu');
    if isempty(SET(NO).Scar) && strcmp(g.CurrentTheme,'scar')
      %g.updateicons('lv');
    end

    g.Run = 0;
    g.cleardatalevelset;

    g.EndoEdgeDetected = false;
    g.EpiEdgeDetected = false;
    g.EndoEDGE0 = [];
    g.EndoEDGE1 = [];
    g.EndoEDGE2 = [];
    g.EndoEDGE3 = [];
    g.EpiEDGE0 = [];
    g.EpiEDGE1 = [];
    g.EpiEDGE2 = [];
    g.EpiEDGE3 = [];
    g.BpInt = -1;
    g.MInt = -1;
    g.TInt = -1;
    g.BalloonLevel = -1;
    g.BALLOON = [];
		g.EndoBalloonForce = [];
		g.EpiBalloonForce = [];

    %Update that makes the function much slower, seems unnecessary
%     drawfunctions('drawimageno',NO);
% 
%     %Switch to a panel in that stack, if possible
%     if force
%       ind=find(g.ViewPanels==NO);
%       if not(isempty(ind))
%         g.CurrentPanel=ind(1);
%       end
%     end
    % end update
       
    drawfunctions('drawthumbnailframes');
    
    if ~isempty(g.Handles.imageaxes) && all(ishandle(g.Handles.imageaxes))
      %Make now one orange
      set(g.Handles.imageaxes,...
        'xcolor',[0 0 0],'ycolor',[0 0 0]);
      set(g.Handles.imageaxes(g.CurrentPanel),...
        'xcolor',g.GUISettings.AxesColor,'ycolor',g.GUISettings.AxesColor,...
        'linewidth',2.5);
    end

     %Enable/Disable play icon etc
    %g.updatetimethings;
    %g.updateaxestables;
%     %Update volumeaxes
    g.updateaxestables('allreports');
%     segment('updatevolume');
%     segment('updateflow');

    %Update title on window
    g.updatetitle;
    end
            
  end
  
  methods(Static)
    
    %--------------------
    function versionhello
    %--------------------
    %Versionhello, overloaded in CVQgui, RVQ, SegmentTransfer, Segment CMR
    %(soon)
    
    getmodule('versionhello');
    if verLessThan('matlab','8.3')
      mywarning(['Segment was developed using MATLAB 8.3 (R2014a). '...
        'You are now running on an older version, which may cause errors.']);
    end
    end
    
    %---------------------
    function orthoview(no)
    %---------------------
    % Create orthogonal view for stack no
    global SET NO 
    if nargin < 1
      no = NO;
    end
    
    SET(no).HLA.slice = round(SET(no).XSize/2);
    SET(no).HLA.maxslice = SET(no).XSize;
    SET(no).HLA.ZoomState = [];
    
    SET(no).VLA.slice = round(SET(no).YSize/2);
    SET(no).VLA.maxslice = SET(no).YSize;
    SET(no).VLA.ZoomState = [];
    
    ang = 2*pi/8;
    SET(no).GLA.angle = ang;
    SET(no).GLA.slice = 0;
    SET(no).GLA.ZoomState = [];
    
    updatetool('orthoview');
    %mywaitbarstart();
    drawfunctions('drawimageview',[no no no no],[2 2],{'ortho','hla','vla','gla'});
    %mywaitbarclose();
    %pause(0.1);
    %updatetool('orthoview');
    end
    
    %---------------------
    function startlog(log)
    %---------------------
    %Start Segment log.
    diary(log);
    end

    %---------------------
    function stoplog
    %---------------------
    %Stop Segment logging.
    diary off;
    end  
    
    %-----------------------------
    function setsectorrotation(no)
    %-----------------------------
    %Overloaded in RVQ GUI
    end
    
    %---------------------------------------
    function filename = generatesavefilename
    %---------------------------------------
    %Called by segment('filesaveallas_Callback'). Overloaded in CVQgui.
    global SET
    filename = [removeforbiddenchars([SET(1).PatientInfo.Name '-' datestr(now,'yyyy-mm-dd')]) '.mat'];
    end
    
    %-----------------------------------
    function ok = manualdraw_Buttonup_roi(no,xr,yr,slice)
    %-----------------------------------
    %Draw new ROI. Called by segment_main('manualdraw_Buttonup')
    %Overloaded in CVQgui and RVQgui.
    
    global SET DATA
       
    ok = true;
    if any(SET(no).RoiCurrent>SET(no).RoiN | SET(no).RoiCurrent<1)
      SET(no).RoiCurrent = SET(no).RoiN;
    end;
    if isempty(SET(no).RoiCurrent)
      SET(no).RoiCurrent = SET(no).RoiN;
    end;
    
    SET(no).RoiN = SET(no).RoiN+1;
    
    SET(no).Roi(SET(no).RoiN).X = repmat(yr',[1 SET(no).TSize]);
    SET(no).Roi(SET(no).RoiN).Y = repmat(xr',[1 SET(no).TSize]);
    SET(no).Roi(SET(no).RoiN).T = 1:SET(no).TSize;
    SET(no).Roi(SET(no).RoiN).Z = slice;
    SET(no).Roi(SET(no).RoiN).Sign = 1;
    if SET(no).RoiN==1
      SET(no).Roi(SET(no).RoiN).Name = sprintf('ROI-%d',SET(no).RoiN);
      SET(no).Roi(SET(no).RoiN).LineSpec = 'b-';
    else
      if length(SET(no).Roi(SET(no).RoiCurrent(end)).Name)>=4 && isequal(SET(no).Roi(SET(no).RoiCurrent(end)).Name(1:4),'ROI-')
        SET(no).Roi(SET(no).RoiN).Name = sprintf('ROI-%d',SET(no).RoiN);
      else
        SET(no).Roi(SET(no).RoiN).Name =SET(no).Roi(SET(no).RoiCurrent(end)).Name;
      end
      SET(no).Roi(SET(no).RoiN).LineSpec = SET(no).Roi(SET(no).RoiCurrent(end)).LineSpec;
    end;
    roi('roiforceapply');
    
    %Calculate area and intensity of ROI
    [~,SET(no).Roi(SET(no).RoiN).Area] = ...
      calcfunctions('calcroiarea',no,SET(no).RoiN);
    [m,sd]=calcfunctions('calcroiintensity',no,SET(no).RoiN);
    SET(no).Roi(SET(no).RoiN).Mean = m;
    SET(no).Roi(SET(no).RoiN).StD = sd;    
    
    SET(no).RoiCurrent = SET(no).RoiN;
    
    %update flow result panel
    if ~isempty(DATA.FlowNO) && ismember(DATA.FlowNO,SET(no).Linked) && ~isempty(DATA.FlowROI)
      if SET(no).RoiN <1
        DATA.FlowROI = [];
      else
        calcfunctions('calcflow',no);
      end
    end
    DATA.updateaxestables('area',no);
    
    end
    %-------------------------------------    
    function setribbonimages(language)
    %-------------------------------
    %This function is used for translation and needs to be updated when
    %adding ribbons.
    global DATA    
    g=DATA;
    
     
     if strcmp(language,'English')
       if strcmp(g.ProgramName,'Segment CMR')
         g.Handles.toggleiconholder.iconCell{3}.cdata=g.Icons.toggleicons.ribbonflowoff;
         g.Handles.toggleiconholder.iconCell{3}.cdataIndent=g.Icons.toggleicons.ribbonflowon;
         g.Handles.toggleiconholder.iconCell{4}.cdata=g.Icons.toggleicons.ribbonviabilityoff;
         g.Handles.toggleiconholder.iconCell{4}.cdataIndent=g.Icons.toggleicons.ribbonviabilityon;
         g.Handles.toggleiconholder.iconCell{6}.cdata=g.Icons.toggleicons.ribbonanalysisoff;
         g.Handles.toggleiconholder.iconCell{6}.cdataIndent=g.Icons.toggleicons.ribbonanalysison;
         g.Handles.toggleiconholder.iconCell{7}.cdata=g.Icons.toggleicons.ribbonimageoff;
         g.Handles.toggleiconholder.iconCell{7}.cdataIndent=g.Icons.toggleicons.ribbonimageon;
       else 
         g.Handles.toggleiconholder.iconCell{4}.cdata=g.Icons.toggleicons.ribbonanalysisoff;
         g.Handles.toggleiconholder.iconCell{4}.cdataIndent=g.Icons.toggleicons.ribbonanalysison;
         g.Handles.toggleiconholder.iconCell{5}.cdata=g.Icons.toggleicons.ribbonimageoff;
         g.Handles.toggleiconholder.iconCell{5}.cdataIndent=g.Icons.toggleicons.ribbonimageon;
       end
         
     
     else
       if strcmp(g.ProgramName,'Segment CMR')
         g.Handles.toggleiconholder.iconCell{3}.cdataIndent=g.Icons.toggleicons.(['ribbonflowon_',language]);
         g.Handles.toggleiconholder.iconCell{3}.cdata=g.Icons.toggleicons.(['ribbonflowoff_',language]);
         g.Handles.toggleiconholder.iconCell{4}.cdata=g.Icons.toggleicons.(['ribbonviabilityoff_',language]);
         g.Handles.toggleiconholder.iconCell{4}.cdataIndent=g.Icons.toggleicons.(['ribbonviabilityon_',language]);
         g.Handles.toggleiconholder.iconCell{6}.cdata=g.Icons.toggleicons.(['ribbonanalysisoff_',language]);
         g.Handles.toggleiconholder.iconCell{6}.cdataIndent=g.Icons.toggleicons.(['ribbonanalysison_',language]);
         g.Handles.toggleiconholder.iconCell{7}.cdata=g.Icons.toggleicons.(['ribbonimageoff_',language]);
         g.Handles.toggleiconholder.iconCell{7}.cdataIndent=g.Icons.toggleicons.(['ribbonimageon_',language]);
       else
         g.Handles.toggleiconholder.iconCell{4}.cdata=g.Icons.toggleicons.(['ribbonanalysisoff_',language]);
         g.Handles.toggleiconholder.iconCell{4}.cdataIndent=g.Icons.toggleicons.(['ribbonanalysison_',language]);
         g.Handles.toggleiconholder.iconCell{5}.cdata=g.Icons.toggleicons.(['ribbonimageoff_',language]);
         g.Handles.toggleiconholder.iconCell{5}.cdataIndent=g.Icons.toggleicons.(['ribbonimageon_',language]);
       end
     end
     
     for i = 1:length(g.Handles.toggleiconholder.iconCell)
       cdataIndent=g.Handles.toggleiconholder.iconCell{i}.cdataIndent;
       cdata=g.Handles.toggleiconholder.iconCell{i}.cdata;
       g.Handles.toggleiconholder.iconCell{i}.generateclickeddisabledandindent(cdata,cdataIndent);
       if ~g.Handles.toggleiconholder.iconCell{i}.enabled
         g.Handles.toggleiconholder.iconCell{i}.cdataDisplay=g.Handles.toggleiconholder.iconCell{i}.cdataDisabled;
       elseif g.Handles.toggleiconholder.iconCell{i}.isindented
         g.Handles.toggleiconholder.iconCell{i}.cdataDisplay=g.Handles.toggleiconholder.iconCell{i}.cdataIndent;
       else
         g.Handles.toggleiconholder.iconCell{i}.cdataDisplay=g.Handles.toggleiconholder.iconCell{i}.cdata;
       end
       
     end
     g.Handles.toggleiconholder.render;
    end
    %----------------------------------------------------------
    function handle = iconcdatahelper(handle,icondata,tiptext)
    %----------------------------------------------------------
    %Helper function to updateicons
    if nargin==2
      set(handle,'cdata',icondata);
    else
      set(handle,'cdata',icondata,'Tooltipstring',translation.dictionary(tiptext));
    end;
    end
    
    %-------------------------------------------
    function iconcallbackhelper(handle,callback)
    %-------------------------------------------
    %Helper function to updateicons.
    set(handle,'callback',sprintf('updatetool(''%s'')',callback));
    end
    
    %---------------------------------------
    function pointshowthisframeonly_Callback 
    %---------------------------------------
    %Callback to make the current point visible in this time frame only.
    %Overloaded in CVQgui.
    global SET NO

    %Use to point to mag data set
    no = NO;
    if ~isempty(SET(NO).Parent)
      no = SET(NO).Parent;
    end;

    % CVQ has this enabled by default for all points, therefore this
    % leads to duplicates in the undo history.
    tools('enableundo',no);

    pos = annotationpoint('pointfind');
    if isnan(pos)
      myfailed('No point found.');
      return;
    end
    SET(no).Point.T(pos) = SET(no).CurrentTimeFrame;
    drawfunctions('drawimageno');
    end
    
    %------------------
    function updateedes
    %------------------
    %Overloaded in CVQgui
    end
        
    %-------------------------------
    function initreportflow(handles) %#ok<INUSD>
    %-------------------------------
    %Overloaded in CVQgui.
    end
    
    %--------------------------
    function c = centercrossdef
    %--------------------------
    %Define color of image center cross
    c = [1 1 1];
    end
    
    %-----------------------------------------------------------
    function h = mywaitbarstart(iter,stri,ignorelimit,fighandle)
    %-----------------------------------------------------------
    %Makes it possible to overload the mywaitbarstart behaviour
    if nargin<2
      myfailed('Expected two or more input arguments.');
      return;
    end;
    
    stri = translation.dictionary(stri);

    myworkon;

    if length(iter)>1
      iter = length(iter);
    end;

    h = [];
    h.iter = iter;
    h.count = 0;
    h.oldfrac = 0;

    if nargin>2
      h.ignorelimit = ignorelimit;
    else
      h.ignorelimit = 1;
    end;

    %Check if empty vector supplied.
    if isempty(h.ignorelimit)
      h.ignorelimit = 1; %default
    end;

    if h.iter>h.ignorelimit
      h.h = waitbar(0,stri);
      %set(h.h,'visible','off');
    else
      h.h = [];
    end;

    %This code does not work with Matlab R2010, there waitbar is not a
    %handle, it is a class.
%    try
%      if nargin<4
%        myadjust(h.h); %Adjust horisontally!
%      else
%        myadjust(h.h,fighandle); %Adjust horisontally!
%      end;
%    catch %#ok<CTCH>
%    end

    %set(h.h,'visible','on');
    flushlog;
    end

    %-------------------------
    function mywaitbarclose(h)
    %------------------------
    %Makes it possible to overload mywaitbar behaviour
    
    if isempty(h.h)
    else
      close(h.h);
    end;

    flushlog;
    myworkoff;
    end;
    
    %--------------------------
    function h = mywaitbarupdate(h,stri)
    %--------------------------
    %Makes it possible to overload mywaitbar behaviour
    h.count = h.count+1;
    if isempty(h.h)
      return;
    end;
    newfrac = round(20*h.count/h.iter);
    if newfrac==h.oldfrac
    else
      h.oldfrac = newfrac;
      if nargin > 1 && ~isempty(stri)
        stri = translation.dictionary(stri);
        waitbar(h.count/h.iter,h.h,stri);
      else
        waitbar(h.count/h.iter,h.h);
      end
    end;

    flushlog;
    
    end;
    
    
    
    %-----------------------------------------------------------
    function h = mywaitbarmainstart(iter,stri,ignorelimit,waitbarfigure)
   %-----------------------------------------------------------
    %Overloaded functions to have the waitbar in the bottom of the main GUI
        
    global DATA
    
    if nargin<2
      myfailed('Expected two or more input arguments.');
      return;
    end;
    
    stri = translation.dictionary(stri);
    myworkon;
    if length(iter)>1
      iter = length(iter);
    end;

    %set parameters
    h = [];
    if nargin >= 4
      h.guihandle = guihandles(waitbarfigure);
    elseif ~isempty(gcbf)
      h.guihandle = guihandles(gcbf); 
    else
      h.guihandle = guihandles(DATA.GUI.Segment.fig);
    end
    h.iter = iter;
    h.count = 0;
    h.oldfrac = 0;
    if nargin>3
      h.ignorelimit = ignorelimit;
    else
      h.ignorelimit = 1;
    end;
    if isempty(h.ignorelimit) %Check if empty vector supplied.
      h.ignorelimit = 1; %default
    end;
    
    %initialize waitbar
    if isfield(h.guihandle,'waitbaraxes')
      axis(h.guihandle.waitbaraxes,'normal');
      h.guihandle.waitbarpatch = patch([0 0 0 0],[0 1 1 0],[1 1 1 1],...
        'parent',h.guihandle.waitbaraxes,'facecolor',[0 1 0]);
      hold(h.guihandle.waitbaraxes,'on');
      h.guihandle.waitbarline = plot(h.guihandle.waitbaraxes,[0 0 1 1 0],[0 1 1 0 0],'g-');
      set(h.guihandle.waitbarline,'linewidth',4);
      hold(h.guihandle.waitbaraxes,'off');
      xlim(h.guihandle.waitbaraxes,[0 1]);
      axis(h.guihandle.waitbaraxes,'off');
      
      %update waitbar
      if h.iter>=h.ignorelimit
        set(h.guihandle.waitbaraxes,'visible','on');
        set(h.guihandle.waitbarline,'visible','on');
        set(h.guihandle.waitbartext,'visible','on');
        set(h.guihandle.waitbarpatch,'visible','on');
        axis(h.guihandle.waitbaraxes,'off');
        set(h.guihandle.waitbartext,'backgroundcolor',DATA.GUISettings.BackgroundColor)
        set(h.guihandle.waitbartext,'foregroundcolor',DATA.GUISettings.TimebarAxesColor)
        set(h.guihandle.waitbartext,'string',stri);
        set(h.guihandle.waitbarpatch,'xdata',[0 0 0 0]);
        drawnow('expose');
      else
        h = [];
      end;
    else
      %if the waitbare axes do not exist in the current GUI, use matlabs
      %waitbar instead
      h.h = waitbar(0,stri);
    end

    flushlog;
    end
    
    %---------------------------------------------
    function h = mywaitbarmainupdate(h,stri)
    %---------------------------------------------
    %Overloaded functions to have the waitbar in the bottom of the main GUI
        
    h.count = h.count+1;
    if isempty(h)
      return;
    end;
    newfrac = round(20*h.count/h.iter);
    if newfrac==h.oldfrac
    else
      h.oldfrac = newfrac;
      f = h.count/h.iter;
      if isfield(h.guihandle,'waitbaraxes')
        set(h.guihandle.waitbarpatch,'xdata',[0 0 f f]);
        if nargin >= 3
          stri = translation.dictionary(stri);
          set(h.guihandle.waitbartext,'string',stri);
        end
        drawnow('expose');
      else
        if nargin >= 3 
          stri = translation.dictionary(stri);
          waitbar(h.count/h.iter,h.h,stri);
        else
          waitbar(h.count/h.iter,h.h);
        end
      end
    end;

    flushlog;
    
    end;
    
    %-----------------------------------
    function mywaitbarmainclose(h)
    %-----------------------------------
    %Overloaded functions to have the waitbar in the bottom of the main GUI
        
    if isempty(h)
    elseif isfield(h.guihandle,'waitbaraxes')
      set(h.guihandle.waitbartext,'visible','off');
      set(h.guihandle.waitbarline,'visible','off');
      set(h.guihandle.waitbarpatch,'visible','off');
      set(h.guihandle.waitbaraxes,'visible','off');
      drawnow('expose');
      clear('h');
    else
      close(h.h);
    end;

    flushlog;
    %myworkoff; mainwaitbar only used for registration which doesn't mess
    %with pointer
    end;
    
    
    %-------------------------
    function copyrvendoforward
    %-------------------------
    %Overloaded in SegmentGUI
    end
            
    %----------------
    function maketest
    %----------------
    %Make test. Overloaded in Segment CMR GUI and RVQ GUI.
    if (isequal(49,getmodule(48,'XT',[],true)))
      path = myuigetdir(pwd,'Select folder for writing test record');
      if ~isequal(path,0)
        internaltools.maketest('segment','allbutcallbacks',path);
      end
    end
    end
    
    %-----------------------------------------------------
    function [type,viewplane,technique] = imagedescription
    %-----------------------------------------------------
    %Define the image types, image view planes and imaging techniques
    %written by Helen Soneson 2009-05-25
    %Moved by Nils Lundahl to maingui.m 2012-07-20
    %Overloaded in RVQ GUI and Segment CT GUI
    
    type = {'General';
      'Cine';
      'Late enhancement';
      'Perfusion Rest';
      'Perfusion Stress';
      'Qflow';
      'Scout';
      'Strain FFE';
      'Strain TFE';
      'Strain from tagging';
      'T1BB';
      'T1 map'
      'T1 map Pre';
      'T1 map Post';
      'T2 map';
      'T2* map';
      'T2Stir';
      'User defined'};
 
    viewplane = {'Unspecified';
      '2CH';
      '3CH';
      '4CH';
      'Sagittal';
      'Coronal';
      'Frontal';
      'Transversal';
      'Short-axis';
      'Short-axis basal';
      'Short-axis mid-ventricular';
      'Short-axis apical';
      'RVOT';
      'Aorta';
      'Pulmonalis';
      'Vena cava inferior';
      'User defined'};
    
    technique = {'Unspecified';
      'MRSSFP';
      'MRDE';
      'MRBB';
      'MRPDW';
      'MRSTIR';
      'MRTOF';
      'MRGE';
      'NM';
      'PT';
      'CTheart';
      'US';
      'User defined'};
    end
    
    %------------------------------------------------
    function menuoptions = getmultiplepatientsoptions
    %------------------------------------------------
    %Return options for a user trying to load data from patient
    %with different name. Overloaded in Segment CMR GUI and RVQ GUI.
    menuoptions = {'Abort loading of new image stacks',...
      'Close previously loaded image stacks',...
      'Load all image stacks (not recommended)'};
    end
    
    %-------------------------------
    function investigationaluselabel
    %-------------------------------
    %In Segment do nothing. The overloaded function in Segment CMR GUI
    %displays warning message that the function is for investigational use
    %only   
    end
    
  end
end