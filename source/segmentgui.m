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
    
    %----------------------------------------
    function initpermanentplaceholder(varargin)
    %--------------------------------------
    g=varargin{1};
      
    iconcell=cell(1,1);
    iconcell{1,1}=myicon('database',g.Handles.permanenticonholder,g.Icons.config.database,'Open file loader',@() segment('fileopen_Callback'),0);
    iconcell{1,end+1}=myicon('databaseadd',g.Handles.permanenticonholder,g.Icons.config.databaseadd,'Save to disc',... 
      @() filemenu('saveall_Callback'),0);
%     iconcell{1,end+1}=myicon(g.Handles.permanenticonholder,g.Icons.config.connect,'Open PACS connection','pacs(''init_Callback'')',0);
%     iconcell{1,end+1}=myicon(g.Handles.permanenticonholder,g.Icons.config.connectadd,'Save to PACS','filemenu(''savetopacs_Callback'')',0);
    iconcell{1,end+1}=myicon('closeall',g.Handles.permanenticonholder,g.Icons.config.closeall,'Close all image stacks',@() segment('filecloseall_Callback'),0);
    
    iconcell{1,end+1}=myicon('panel1',g.Handles.permanenticonholder,g.Icons.config.panel1,'View one image panel',@() drawfunctions('drawimageview',[],[1,1]),1,1);
    iconcell{1,end+1}=myicon('panel2',g.Handles.permanenticonholder,g.Icons.config.panel2,'View two image panels',@() drawfunctions('drawimageview',[],[1,2]),1,1);
    iconcell{1,end+1}=myicon('panel2x1',g.Handles.permanenticonholder,g.Icons.config.panel2x1,'View two image panels',@() drawfunctions('drawimageview',[],[2,1]),1,1);
    iconcell{1,end+1}=myicon('panel3x1',g.Handles.permanenticonholder,g.Icons.config.panel3x1,'View three image panels',@() drawfunctions('drawimageview',[],[3,1]),1,1);
    iconcell{1,end+1}=myicon('panel3',g.Handles.permanenticonholder,g.Icons.config.panel3,'View three image panels',@() drawfunctions('drawimageview',[],[1,3]),1,1);
    iconcell{1,end+1}=myicon('panel4',g.Handles.permanenticonholder,g.Icons.config.panel4,'View four image panels',@() drawfunctions('drawimageview',[],[2,2]),1,1);
    iconcell{1,end+1}=myicon('panel6',g.Handles.permanenticonholder,g.Icons.config.panel6,'View six image panels',@()  drawfunctions('drawimageview',[],[2,3]),1,1);
    iconcell{1,end+1}=myicon('orthoview',g.Handles.permanenticonholder,g.Icons.config.orthoview,'Orthogonal view',@() segment('orthoview'),1,1);
    iconcell{1,end+1}=myicon('saveview',g.Handles.permanenticonholder,g.Icons.config.saveview,'Save view',@() segmentview,0);
    
    iconcell{1,end+1}=myicon('viewone',g.Handles.permanenticonholder,g.Icons.config.viewone,'View one slice',@() segment('viewimage_Callback','one'),1,2);
    iconcell{1,end+1}=myicon('viewall',g.Handles.permanenticonholder,g.Icons.config.viewall,'View all slices',@() segment('viewimage_Callback','montage'),1,2);
    iconcell{1,end+1}=myicon('viewrow',g.Handles.permanenticonholder,g.Icons.config.viewrow,'View all slices in 2 rows',@() segment('viewimage_Callback','montagerow'),1,2);
    
    iconcell{1,end+1}=myicon('undo',g.Handles.permanenticonholder,g.Icons.config.undo,'Undo last operation',@() tools('undosegmentation_Callback'),0);
    iconcell{1,end+1}=myicon('refresh',g.Handles.permanenticonholder,g.Icons.config.refresh,'Refresh image view',@() segment('viewrefreshall_Callback'),0);
    
    iconcell{1,end+1}=myicon('play',g.Handles.permanenticonholder,g.Icons.config.play,'Play movie of all image stacks',@() segment('playall_Callback'),2);
    iconcell{1,end+1}=myicon('next',g.Handles.permanenticonholder,g.Icons.config.next,'Next time frame for all image stacks',@() segment('nextallframe_Callback'),0);
    iconcell{1,end+1}=myicon('prev',g.Handles.permanenticonholder,g.Icons.config.prev,'Previous time frame for all image stacks',@() segment('previousallframe_Callback'),0);
    iconcell{1,end+1}=myicon('faster',g.Handles.permanenticonholder,g.Icons.config.faster,'Faster frame rate',@() segment('fasterframerate_Callback'),0);
    iconcell{1,end+1}=myicon('slower',g.Handles.permanenticonholder,g.Icons.config.slower,'Slower frame rate',@() segment('slowerframerate_Callback'),0);
    
    iconcell{1,end+1}=myicon('hideall',g.Handles.permanenticonholder,g.Icons.config.hideall,'Hide all overlays (segmentation, point, text,...)',@() segment('viewhideall_Callback'),2);
    iconcell{1,end+1}=myicon('clearall',g.Handles.permanenticonholder,g.Icons.config.clearall,'Clear all segmentation in current image stack',@() segment('segmentclearall_Callback'),0);
    iconcell{1,end+1}=myicon('clearalledes',g.Handles.permanenticonholder,g.Icons.config.clearalledes,'Clear all segmentation except in ED and ES',@() segment('segmentclearallbutsystolediastole_Callback'),0); 
    iconcell{1,end+1}=myicon('zoomin',g.Handles.permanenticonholder,g.Icons.config.zoomin,'Zoom in',@() segment('viewzoomin_Callback'),0);
    iconcell{1,end+1}=myicon('zoomout',g.Handles.permanenticonholder,g.Icons.config.zoomout,'Zoom out',@() segment('viewzoomout_Callback'),0);
    iconcell{1,end+1}=myicon('autozoom',g.Handles.permanenticonholder,g.Icons.config.autozoom,'Auto zoom',@() segment('autozoom'),0);
    iconcell{1,end+1}=myicon('colorbar',g.Handles.permanenticonholder,g.Icons.config.colorbar,'Hide colorbar',@() segment('viewhidecolorbar_Callback'),2);
    iconcell{1,end+1}=myicon('viewpixels',g.Handles.permanenticonholder,g.Icons.config.viewpixels,'Show image pixels',@() segment('viewinterp_Callback'),2);
%     iconcell{1,end+1}=myicon('reportsheet',g.Handles.permanenticonholder,g.Icons.config.reportsheet,'Open Report sheet generation',@() reporter.reportsheet,0);
    iconcell{1,end+1}=myicon('savescreen',g.Handles.permanenticonholder,g.Icons.config.savescreen,'Save screen shot',@() export('screenshot_Callback'),0);
    
    iconcell{1,end+1}=myicon('settingsgeneral',g.Handles.permanenticonholder,g.Icons.config.settingsgeneral,'Set general preferences',@() segpref,0);
    iconcell{1,end+1}=myicon('settingssystem',g.Handles.permanenticonholder,g.Icons.config.settingsdatabase,'Set patient database preferences',@() segpref('advancedsettings_Callback'),0);
    iconcell{1,end+1}=myicon('settingspacs',g.Handles.permanenticonholder,g.Icons.config.settingspacs,'Set PACS connection preferences',@() pacspref,0);
    g.Handles.permanenticonholder.add(iconcell);
    
    pos=plotboxpos(g.Handles.permanenticonholder.axeshandle);
    currentpos=get(g.Handles.permanenticonholder.axeshandle,'position');
    set(g.Handles.permanenticonholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);
    set(g.Handles.iconuipanel,'visible','on')
    
    end
    
    %---------------------------------------------
    function initconfigplaceholder(varargin)
      %--------------------------------------------
      g=varargin{1};
      
      %check if using new gui version
%       if all([isfield(g.Icons,'lviconcell'),isfield(g.Icons,'rviconcell'),...
%         isfield(g.Icons,'analysisviconcell'),isfield(g.Icons,'roiflowiconcell'),...
%         isfield(g.Icons,'viabilityiconcell'),isfield(g.Icons,'imageiconcell')]);
%        
%       g.lviconcell=cell(1,1);
%     g.rviconcell=cell(1,1);
%     g.analysisiconcell=cell(1,1);
%     g.roiflowiconcell=cell(1,1);
%     g.viabilityiconcell=cell(1,1);
%     g.imageiconcell=cell(1,1);
%       else
    %initcells
    lviconcell=cell(1,1);
    rviconcell=cell(1,1);
    analysisiconcell=cell(1,1);
    roiflowiconcell=cell(1,1);
    viabilityiconcell=cell(1,1);
    imageiconcell=cell(1,1);
    %LV
    lviconcell{1,1}=myicon('select',g.Handles.configiconholder,g.Icons.config.select,'Select image stack or object',@() updatetool('select'));
    lviconcell{1,end+1}=myicon('move',g.Handles.configiconholder,g.Icons.config.move,'Translate contour',@() updatetool('move'));
    lviconcell{1,end+1}=myicon('scale',g.Handles.configiconholder,g.Icons.config.scale,'Scale object',@() updatetool('scale'));
    lviconcell{1,end+1}=myicon('contrastbrightness',g.Handles.configiconholder,g.Icons.config.contrastbrightness,'Manually change contrast and brightness',@() updatetool('contrast'));
    lviconcell{1,end+1}=myicon('autocontrast',g.Handles.configiconholder,g.Icons.config.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    lviconcell{1,end+1}=myicon('lvstack',g.Handles.configiconholder,g.Icons.config.lvstack,'Go to LV stack',@() segment('viewspecial_Callback','lv'),0);
    lviconcell{1,end+1}=myicon('moveall',g.Handles.configiconholder,g.Icons.config.moveall,'Translate all contours',@() updatetool('moveall'));    
    lviconcell{1,end+1}=myicon('autolv',g.Handles.configiconholder,g.Icons.config.autolv,'Automatic LV segmentation',@() updatetool('autolv'),0);
    lviconcell{1,end+1}=myicon('endopen',g.Handles.configiconholder,g.Icons.config.endopen,'Endo pen',@() updatetool('drawendo'));
    lviconcell{1,end+1}=myicon('epipen',g.Handles.configiconholder,g.Icons.config.epipen,'Epi pen',@() updatetool('drawepi'));
    lviconcell{1,end+1}=myicon('interpendo',g.Handles.configiconholder,g.Icons.config.interpendo,'Set interpolation points for Endo',@() updatetool('interpendo'));
    lviconcell{1,end+1}=myicon('interpepi',g.Handles.configiconholder,g.Icons.config.interpepi,'Set interpolation points for Epi',@() updatetool('interpepi'));
    lviconcell{1,end+1}=myicon('smooth',g.Handles.configiconholder,g.Icons.config.smooth,'Smooth latest contour (O)',@() tools('smoothsegmentation_Callback'),0);
          
    lviconcell{1,end+1}=myicon('interpseginslice',g.Handles.configiconholder,g.Icons.config.interpseginslice,'Interpolate segmentation over slices',@() lv('interpolatedelineation_Callback'),0);
    lviconcell{1,end+1}=myicon('interpsegintime',g.Handles.configiconholder,g.Icons.config.interpsegintime,'Interpolate segmentation in time for selected slices',@() segmentation('interpolatedelineationovertime_Callback'),0);
    lviconcell{1,end+1}=myicon('refineendo',g.Handles.configiconholder,g.Icons.config.refineendo,'Refine Endo',@() lvpeter('segmentrefineendo_Callback'),0);
    lviconcell{1,end+1}=myicon('refineepi',g.Handles.configiconholder,g.Icons.config.refineepi,'Refine Epi',@() lvpeter('segmentrefineepi_Callback'),0);
    lviconcell{1,end+1}=myicon('propagateendo',g.Handles.configiconholder,g.Icons.config.propagateendo,'Propagate endo forward in time', @() lvpeter('segmentpropagateendo_Callback'),0);
    lviconcell{1,end+1}=myicon('propagateepi',g.Handles.configiconholder,g.Icons.config.propagateepi,'Propagate epi forward in time',@() lvpeter('segmentpropagateepi_Callback'),0);
    
    lviconcell{1,end+1}=myicon('contractendo',g.Handles.configiconholder,g.Icons.config.contractendo,'Contract Endo segmentation',@() lv('segmentexpandcontract_Callback',-1,'endo'),0);
    lviconcell{1,end+1}=myicon('expandendo',g.Handles.configiconholder,g.Icons.config.expandendo,'Expand Endo segmentation',@() lv('segmentexpandcontract_Callback',1,'endo'),0);
    lviconcell{1,end+1}=myicon('contractepi',g.Handles.configiconholder,g.Icons.config.contractepi,'Contract Epi segmentation',@() lv('segmentexpandcontract_Callback',-1,'epi'),0);
    lviconcell{1,end+1}=myicon('expandepi',g.Handles.configiconholder,g.Icons.config.expandepi,'Expand Epi segmentation',@() lv('segmentexpandcontract_Callback',1,'epi'),0);
    lviconcell{1,end+1}=myicon('evenoutwall',g.Handles.configiconholder,g.Icons.config.evenoutwall,'Even out wall',@() segment('smoothendowall_Callback'),0);
    lviconcell{1,end+1}=myicon('copylvup',g.Handles.configiconholder,g.Icons.config.copylvup,'Copy LV upwards and refine',@()tools('copyupward_Callback'),0);
    lviconcell{1,end+1}=myicon('copylvdown',g.Handles.configiconholder,g.Icons.config.copylvdown,'Copy LV downwards and refine',@()tools('copydownward_Callback'),0);
        
    lviconcell{1,end+1}=myicon('hidelv',g.Handles.configiconholder,g.Icons.config.hidelv,'Hide LV segmentation',@() segment('viewhidelv_Callback'),2);
    lviconcell{1,end+1}=myicon('hideinterp',g.Handles.configiconholder,g.Icons.config.hideinterp,'Hide interpolation points',@() segment('viewhideinterp_Callback'),2);
    lviconcell{1,end+1}=myicon('clearalllv',g.Handles.configiconholder,g.Icons.config.clearalllv,'Clear all LV segmentation',@() segment('segmentclearalllv_Callback'),0);
    lviconcell{1,end+1}=myicon('clearendo',g.Handles.configiconholder,g.Icons.config.clearendo,'Clear LV endo in selected slices according to mode',@() segmentation('clearslicesthis_Callback',1,0,0,0),0);
    lviconcell{1,end+1}=myicon('clearepi',g.Handles.configiconholder,g.Icons.config.clearepi,'Clear LV epi in selected slices according to mode',@() segmentation('clearslicesthis_Callback',0,1,0,0),0);
    lviconcell{1,end+1}=myicon('volumecurve',g.Handles.configiconholder,g.Icons.config.volumecurve,'Plot Volume Curve',@() lvpeter('plotvolumecurve'),0);
   
    g.Icons.lviconcell=lviconcell;
    
    %RV
    rviconcell{1,1}=myicon('select',g.Handles.configiconholder,g.Icons.config.select,'Select image stack or object',@() updatetool('select'));
    rviconcell{1,end+1}=myicon('move',g.Handles.configiconholder,g.Icons.config.move,'Translate contour', @() updatetool('move'));
    rviconcell{1,end+1}=myicon('scale',g.Handles.configiconholder,g.Icons.config.scale,'Scale object',@() updatetool('scale'));
    rviconcell{1,end+1}=myicon('contrastbrightness',g.Handles.configiconholder,g.Icons.config.contrastbrightness,'Manually change contrast and brightness',@() updatetool('contrast'));
    rviconcell{1,end+1}=myicon('autocontrast',g.Handles.configiconholder,g.Icons.config.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    rviconcell{1,end+1}=myicon('rvstack',g.Handles.configiconholder,g.Icons.config.rvstack,'Go to RV stack',@() segment('viewspecial_Callback','rv'),0);
    rviconcell{1,end+1} = lviconcell{1,7}; %moveall
    rviconcell{1,end+1}=myicon('autorvendo',g.Handles.configiconholder,g.Icons.config.autorvendo,'Automatic RV Endo segmentation',@() updatetool('autorvendo'),0);
    rviconcell{1,end+1}=myicon('rvendopen',g.Handles.configiconholder,g.Icons.config.rvendopen,'RV Endo pen',@() updatetool('drawrvendo'));
    rviconcell{1,end+1}=myicon('rvepipen',g.Handles.configiconholder,g.Icons.config.rvepipen,'RV Epi pen',@() updatetool('drawrvepi'));
    rviconcell{1,end+1}=myicon('interprvendo',g.Handles.configiconholder,g.Icons.config.interprvendo,'Set interpolation points for RV Endo',@() updatetool('interprvendo'));
    rviconcell{1,end+1}=myicon('interprvepi',g.Handles.configiconholder,g.Icons.config.interprvepi,'Set interpolation points for RV Epi',@() updatetool('interprvepi'));
    rviconcell{1,end+1}=myicon('refinervendo',g.Handles.configiconholder,g.Icons.config.refinervendo,'Refine RV Endo',@() rv('segmentrefinervendo_Callback'),0);
    %need icon
    rviconcell{1,end+1}=myicon('copyrvup',g.Handles.configiconholder,g.Icons.config.copyrvup,'Copy RV endo upwards',@()tools('copyupward_Callback','endo',false,false),0);
    rviconcell{1,end+1}=myicon('copyrvdown',g.Handles.configiconholder,g.Icons.config.copyrvdown,'Copy RV endo downwards',@()tools('copydownward_Callback','endo',false,false),0);
    
    rviconcell{1,end+1}=myicon('interpsegintime',g.Handles.configiconholder,g.Icons.config.interpsegintime,'Interpolate segmentation in time',@() segmentation('interpolatedelineationovertime_Callback'),0);
    rviconcell{1,end+1}=myicon('interpseginslice',g.Handles.configiconholder,g.Icons.config.interpseginslice,'Interpolate segmentation over slices',@() lv('interpolatedelineation_Callback'),0);
    
    rviconcell{1,end+1}=myicon('hiderv',g.Handles.configiconholder,g.Icons.config.hiderv,'Hide RV segmentation', @() segment('viewhiderv_Callback'),2);
    rviconcell{1,end+1}=lviconcell{end-4};%myicon('hideinterp',g.Handles.configiconholder,g.Icons.config.hideinterp,'Hide interpolation points',@() segment('viewhideinterp_Callback'),2);
    rviconcell{1,end+1}=myicon('clearallrv',g.Handles.configiconholder,g.Icons.config.clearallrv,'Clear all RV segmentation',@() segment('segmentclearallrv_Callback'),0);
    rviconcell{1,end+1}=myicon('clearrvendo',g.Handles.configiconholder,g.Icons.config.clearrv,'Clear RV endo in selected slices according to mode',@() segmentation('clearslicesthis_Callback',0,0,1,0),0);
    rviconcell{1,end+1}=myicon('clearrvepi',g.Handles.configiconholder,g.Icons.config.clearrvepi,'Clear RV epi in selected slices according to mode',@() segmentation('clearslicesthis_Callback',0,0,0,1),0);
    rviconcell{1,end+1}=myicon('volumecurve',g.Handles.configiconholder,g.Icons.config.volumecurve,'Plot Volume Curve',@() lvpeter('plotvolumecurve'),0);
    
    g.Icons.rviconcell=rviconcell;
    
    %ROIFLOW
    
    roiflowiconcell{1,1}=myicon('select',g.Handles.configiconholder,g.Icons.config.select,'Select image stack or object',@() updatetool('select'));
    roiflowiconcell{1,end+1}=myicon('move',g.Handles.configiconholder,g.Icons.config.move,'Translate contour',@() updatetool('move'));
    roiflowiconcell{1,end+1}=myicon('scale',g.Handles.configiconholder,g.Icons.config.scale,'Scale object',@() updatetool('scale'));
    %roiflowiconcell{1,end+1}=myicon('scaleROI',g.Handles.configiconholder,g.Icons.config.scaleROI,'Scale ROI',@() updatetool('scaleROI'));
    roiflowiconcell{1,end+1}=myicon('contrastbrightness',g.Handles.configiconholder,g.Icons.config.contrastbrightness,'Manually change contrast and brightness',@() updatetool('contrast'));
    roiflowiconcell{1,end+1}=myicon('autocontrast',g.Handles.configiconholder,g.Icons.config.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    roiflowiconcell{1,end+1}=myicon('flowstack',g.Handles.configiconholder,g.Icons.config.flowstack,'Go to flow stack',@() segment('viewspecial_Callback','flow'),0);
    roiflowiconcell{1,end+1}=myicon('putroi',g.Handles.configiconholder,g.Icons.config.putroi,'Place ROI',@() updatetool('putroi'));
    roiflowiconcell{1,end+1}=myicon('roipen',g.Handles.configiconholder,g.Icons.config.roipen,'ROI pen',@() updatetool('drawroi'));
    roiflowiconcell{1,end+1}=myicon('roipen',g.Handles.configiconholder,g.Icons.config.trackingvessel,'Track vessel in all time frames',@() vesselsnake_flowtrackroi('flowtrackroi'),0);
    roiflowiconcell{1,end+1}=myicon('refineroi',g.Handles.configiconholder,g.Icons.config.refineroi,'Refine ROI',@() flow('flowrefine_Callback'),0);
    roiflowiconcell{1,end+1}=myicon('refineroinext',g.Handles.configiconholder,g.Icons.config.refineroinext,'Propagate ROI to next timeframe',@() flow('flowpropagate_Callback'),0);
    roiflowiconcell{1,end+1}=myicon('contractroi',g.Handles.configiconholder,g.Icons.config.contractroi,'Contract ROI',@() roi('expandcontract_Callback',-1),0);
    roiflowiconcell{1,end+1}=myicon('expandroi',g.Handles.configiconholder,g.Icons.config.expandroi,'Expand ROI',@() roi('expandcontract_Callback',1),0);
    roiflowiconcell{1,end+1}=myicon('unwrap',g.Handles.configiconholder,g.Icons.config.unwrap,'Unwrap flow',@() flowunwrap,0);
    roiflowiconcell{1,end+1}=myicon('palette',g.Handles.configiconholder,g.Icons.config.palette,'Set ROI color',@() roi('roisetcolor_Callback'),0);
    roiflowiconcell{1,end+1}=myicon('text',g.Handles.configiconholder,g.Icons.config.text,'Set ROI label',@() roi('roisetlabel_Callback'),0);
    roiflowiconcell{1,end+1}=myicon('plotflow',g.Handles.configiconholder,g.Icons.config.plotflow,'Plot flow',@() reportflow,0);
    roiflowiconcell{1,end+1}=myicon('hideroi',g.Handles.configiconholder,g.Icons.config.hideroi,'Hide ROI',@() segment('viewhideroi_Callback'),2);
    roiflowiconcell{1,end+1}=myicon('clearroi',g.Handles.configiconholder,g.Icons.config.clearroi,'Clear selected ROIs',@() roi('roidelete_Callback'),0);  
    roiflowiconcell{1,end+1}=myicon('clearallroi',g.Handles.configiconholder,g.Icons.config.clearallroi,'Clear all ROIs',@() roi('roiclearall_Callback'),0); 
    g.Icons.roiflowiconcell=roiflowiconcell;
    
    %Viablility
    viabilityiconcell{1,1}=myicon('select',g.Handles.configiconholder,g.Icons.config.select,'Select image stack or object',@() updatetool('select'));
    viabilityiconcell{1,end+1}=myicon('move',g.Handles.configiconholder,g.Icons.config.move,'Translate contour',@() updatetool('move'));
    viabilityiconcell{1,end+1}=myicon('scale',g.Handles.configiconholder,g.Icons.config.scale,'Scale object',@() updatetool('scale'));
    viabilityiconcell{1,end+1}=myicon('contrastbrightness',g.Handles.configiconholder,g.Icons.config.contrastbrightness,'Manually change contrast and brightness',@() updatetool('contrast'));
    viabilityiconcell{1,end+1}=myicon('autocontrast',g.Handles.configiconholder,g.Icons.config.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    viabilityiconcell{1,end+1}=myicon('scarstack',g.Handles.configiconholder,g.Icons.config.scarstack,'Go to scar stack',@() segment('viewspecial_Callback','cinescar'),0);
    viabilityiconcell{1,end+1}=myicon('importfromother',g.Handles.configiconholder,g.Icons.config.importfromother,'Import LV segmentation from cine to scar image stack',@() segmentation('importfromcine2scar_Callback'),0);
    viabilityiconcell{1,end+1}=myicon('moveall',g.Handles.configiconholder,g.Icons.config.moveall,'Translate all contours',@() updatetool('moveall'));    
    viabilityiconcell{1,end+1}=myicon('endopen',g.Handles.configiconholder,g.Icons.config.endopen,'Endo pen',@() updatetool('drawendo'));
    viabilityiconcell{1,end+1}=myicon('epipen',g.Handles.configiconholder,g.Icons.config.epipen,'Epi pen',@() updatetool('drawepi'));
    viabilityiconcell{1,end+1}=myicon('interpendo',g.Handles.configiconholder,g.Icons.config.interpendo,'Set interpolation points for Endo',@() updatetool('interpendo'));
    viabilityiconcell{1,end+1}=myicon('interpepi',g.Handles.configiconholder,g.Icons.config.interpepi,'Set interpolation points for Epi',@() updatetool('interpepi'));
    viabilityiconcell{1,end+1}=myicon('autoscar',g.Handles.configiconholder,g.Icons.config.autoscar,'Auto scar',@() updatetool('autoscar'),0);
    viabilityiconcell{1,end+1}=myicon('scarpen',g.Handles.configiconholder,g.Icons.config.scarpen,'Draw scar',@() updatetool('drawscar'));
    viabilityiconcell{1,end+1}=myicon('mopen',g.Handles.configiconholder,g.Icons.config.mopen,'Draw MO',@() updatetool('drawmo'));
    viabilityiconcell{1,end+1}=myicon('rubberscar',g.Handles.configiconholder,g.Icons.config.rubberscar,'Manually remove scar segmentation',@() updatetool('drawrubberpen'));
    viabilityiconcell{1,end+1}=myicon('automar',g.Handles.configiconholder,g.Icons.config.automar,'Auto MaR',@() updatetool('automar'),0);
    viabilityiconcell{1,end+1}=myicon('marpen',g.Handles.configiconholder,g.Icons.config.marpen,'Draw MaR',@() updatetool('drawmarpen'));
    viabilityiconcell{1,end+1}=myicon('rubbermar',g.Handles.configiconholder,g.Icons.config.rubbermar,'Manually remove MaR segmentation',@() updatetool('drawmarrubberpen'));
    viabilityiconcell{1,end+1}=myicon('hidescar',g.Handles.configiconholder,g.Icons.config.hidescar,'Hide scar segmentation',@() segment('viewhidescar_Callback'),2);
    viabilityiconcell{1,end+1}=myicon('hidescarextent',g.Handles.configiconholder,g.Icons.config.hidescarextent,'Hide scar segmentation extent',@() segment('viewhidescarextent_Callback'),2);
    viabilityiconcell{1,end+1}=myicon('hidescarmanual',g.Handles.configiconholder,g.Icons.config.hidescarmanual,'Hide manual scar interaction',@() segment('viewhidemanualinteraction_Callback'),2);
    viabilityiconcell{1,end+1}=myicon('hidemar',g.Handles.configiconholder,g.Icons.config.hidemar,'Hide MaR segmentation',@() segment('viewhidemar_Callback'),2);
    viabilityiconcell{1,end+1}=myicon('clearscar',g.Handles.configiconholder,g.Icons.config.clearscar,'Clear scar segmentation',@() viability('viabilityclear_Callback'),0);  
    viabilityiconcell{1,end+1}=myicon('clearmar',g.Handles.configiconholder,g.Icons.config.clearmar,'Clear MaR segmentation',@() mar('clearall_Callback'),0);  
    g.Icons.viabilityiconcell=viabilityiconcell;
    
    %Analysis
    analysisiconcell{1,1}=myicon('select',g.Handles.configiconholder,g.Icons.config.select,'Select image stack or object',@() updatetool('select'));
    analysisiconcell{1,end+1}=myicon('move',g.Handles.configiconholder,g.Icons.config.move,'Translate contour',@() updatetool('move'));
    analysisiconcell{1,end+1}=myicon('scale',g.Handles.configiconholder,g.Icons.config.scale,'Scale object',@() updatetool('scale'));
    analysisiconcell{1,end+1}=myicon('contrastbrightness',g.Handles.configiconholder,g.Icons.config.contrastbrightness,'Manually change contrast and brightness',@() updatetool('contrast'));
    analysisiconcell{1,end+1}=myicon('autocontrast',g.Handles.configiconholder,g.Icons.config.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    analysisiconcell{1,end+1}=myicon('perfusionstack',g.Handles.configiconholder,g.Icons.config.perfusionstack,'Go to perfusion stacks',@() segment('viewspecial_Callback','perfusion'),0);
    analysisiconcell{1,end+1}=myicon('importfromother',g.Handles.configiconholder,g.Icons.config.importfromother,'Import LV segmentation from other image stack with snap',@() segmentation('importsegmentationwithsnap_Callback'),0);
    analysisiconcell{1,end+1}=myicon('measure',g.Handles.configiconholder,g.Icons.config.measure,'Place Measurement',@() updatetool('measure'));
    analysisiconcell{1,end+1}=myicon('point',g.Handles.configiconholder,g.Icons.config.point,'Place Annotation point',@() updatetool('point'));
    analysisiconcell{1,end+1}=myicon('roipen',g.Handles.configiconholder,g.Icons.config.roipen,'ROI pen',@() updatetool('drawroi'));
    analysisiconcell{1,end+1}=myicon('addroiinlv',g.Handles.configiconholder,g.Icons.config.addroiinlv,'Add ROIs to sector of LV wall in selected slices',@() roi('roiaddinsector_Callback'),0);
    analysisiconcell{1,end+1}=myicon('bullseye',g.Handles.configiconholder,g.Icons.config.bullseye,'Bullseye plot interface',@() reportbullseye,0);
    analysisiconcell{1,end+1}=myicon('AVPD',g.Handles.configiconholder,g.Icons.config.AVPD,'AV plane displacement',@() avplane,0);
    
    analysisiconcell{1,end+1}=myicon('perfusion',g.Handles.configiconholder,g.Icons.config.perfusion, 'Perfusion analysis',@() perfusion.perfusion,0);
    analysisiconcell{1,end+1}=myicon('perfusionscoring',g.Handles.configiconholder,g.Icons.config.perfusionscoring, 'Perfusion scoring',@() perfusion.perfusionscoring,0);
    analysisiconcell{1,end+1}=myicon('reportperslice',g.Handles.configiconholder,g.Icons.config.reportperslice, 'Report per slice',@() slicereport,0);
    analysisiconcell{1,end+1}=myicon('model3d',g.Handles.configiconholder,g.Icons.config.model3d, 'Show 3D model',@() report3dmodel,0);
    analysisiconcell{1,end+1}=myicon('hidemeasure',g.Handles.configiconholder,g.Icons.config.hidemeasure,'Hide measurements',@() segment('viewhidemeasures_Callback'),2);
    analysisiconcell{1,end+1}=myicon('hidepoint',g.Handles.configiconholder,g.Icons.config.hidepoint,'Hide annotation points',@() segment('viewhidepoints_Callback'),2);
    analysisiconcell{1,end+1}= roiflowiconcell{1,end-2};%myicon('hideroi',g.Handles.configiconholder,g.Icons.config.hideroi,'Hide ROI',@() segment('viewhideroi_Callback'),2);
    analysisiconcell{1,end+1}=myicon('clearmeasure',g.Handles.configiconholder,g.Icons.config.clearmeasure,'Clear all measurements',@() segment('measureclearall_Callback'),0);  
    analysisiconcell{1,end+1}=myicon('clearpoint',g.Handles.configiconholder,g.Icons.config.clearpoint,'Clear all annotation points',@() annotationpoint('pointclearall_Callback'),0); 
    analysisiconcell{1,end+1}=myicon('clearroi',g.Handles.configiconholder,g.Icons.config.clearroi,'Clear selected ROIs',@() roi('roidelete_Callback'),0);  
    analysisiconcell{1,end+1}=myicon('clearallroi',g.Handles.configiconholder,g.Icons.config.clearallroi,'Clear all ROIs',@() roi('roiclearall_Callback'),0); 
    g.Icons.analysisiconcell=analysisiconcell;
    
    %Image
    imageiconcell{1,1}=myicon('select',g.Handles.configiconholder,g.Icons.config.select,'Select image stack or object',@() updatetool('select'));
    imageiconcell{1,end+1}=myicon('move',g.Handles.configiconholder,g.Icons.config.move,'Translate contour',@() updatetool('move'));
    imageiconcell{1,end+1}=myicon('scale',g.Handles.configiconholder,g.Icons.config.scale,'Scale object',@() updatetool('scale'));
    imageiconcell{1,end+1}=myicon('contrastbrightness',g.Handles.configiconholder,g.Icons.config.contrastbrightness,'Manually change contrast and brightness',@() updatetool('contrast'));
    imageiconcell{1,end+1}=myicon('autocontrast',g.Handles.configiconholder,g.Icons.config.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    imageiconcell{1,end+1}=myicon('resetlight',g.Handles.configiconholder,g.Icons.config.resetlight,'Reset contrast and brightness',@() segment('resetlight_Callback'),0);
    imageiconcell{1,end+1}=myicon('autocontrastall',g.Handles.configiconholder,g.Icons.config.autocontrastall,'Set contrast and brightness to predefined values for all images',@() segment('autocontrastall_Callback'),0);
    imageiconcell{1,end+1}=myicon('resetlightall',g.Handles.configiconholder,g.Icons.config.resetlightall,'Reset contrast and brightness for all imagestacks',@() segment('resetlightall_Callback'),0);
    imageiconcell{1,end+1}=myicon('crop',g.Handles.configiconholder,g.Icons.config.crop,'Manual crop',@() updatetool('crop'));
    imageiconcell{1,end+1}=myicon('cropall',g.Handles.configiconholder,g.Icons.config.cropall,'Auto crop all',@() updatetool('autocropall'),0);
    imageiconcell{1,end+1}=myicon('croplv',g.Handles.configiconholder,g.Icons.config.croplv,'LV crop',@() lvsegmentation('croplvall',0),0);
    imageiconcell{1,end+1}=myicon('cineplay',g.Handles.configiconholder,g.Icons.config.cineplay,'Open cine tool',@() segment('cinetool_Callback'),2);
    imageiconcell{1,end+1}=myicon('movie',g.Handles.configiconholder,g.Icons.config.movie,'Open movie tool',@() export('exportmovierecorder_Callback'),0);
    imageiconcell{1,end+1}=myicon('click3d',g.Handles.configiconholder,g.Icons.config.click3d,'Set 3D point',@() updatetool('click3d'));
    imageiconcell{1,end+1}=myicon('rotate90',g.Handles.configiconholder,g.Icons.config.rotate90,'Rotate 90 degrees clockwise',@() tools('rotate90right_Callback'),0);
    imageiconcell{1,end+1}=myicon('mpr',g.Handles.configiconholder,g.Icons.config.mpr,'Reconstruct image stack',@() reformater,0);
    imageiconcell{1,end+1}=myicon('mergestacks',g.Handles.configiconholder,g.Icons.config.mergestacks,'Merge stacks',@() mergestacks,0);    
    imageiconcell{1,end+1}=myicon('imageinfo',g.Handles.configiconholder,g.Icons.config.imageinfo,'View and adjust image info',@() tools('imageinfo_Callback'),0);
    imageiconcell{1,end+1}=myicon('patientinfo',g.Handles.configiconholder,g.Icons.config.patientinfo,'View and adjust patient info',@() tools('viewpatientinfo_Callback'),0);    
    imageiconcell{1,end+1}=myicon('hidetext',g.Handles.configiconholder,g.Icons.config.hidetext,'Hide text',@() segment('viewhidetext_Callback'),2);
    imageiconcell{1,end+1}=myicon('hideplus',g.Handles.configiconholder,g.Icons.config.hideplus,'Hide center cross',@() segment('viewhideplus_Callback'),2);
    imageiconcell{1,end+1}=myicon('hideintersections',g.Handles.configiconholder,g.Icons.config.hideintersections,'Hide intersection lines',@() segment('viewhideinterp_Callback'),2);
    imageiconcell{1,end+1}=myicon('hideothercontour',g.Handles.configiconholder,g.Icons.config.hideothercontour,'Hide other contour points',@() segment('viewhideothercontour_Callback'),2);
    g.Icons.imageiconcell=imageiconcell;
    
    %TXmap
    txmapiconcell{1,1}=myicon('select',g.Handles.configiconholder,g.Icons.config.select,'Select image stack or object',@() updatetool('select'));
    txmapiconcell{1,end+1}=myicon('move',g.Handles.configiconholder,g.Icons.config.move,'Translate contour',@() updatetool('move'));
    txmapiconcell{1,end+1}=myicon('scale',g.Handles.configiconholder,g.Icons.config.scale,'Scale object',@() updatetool('scale'));
    txmapiconcell{1,end+1}=myicon('contrastbrightness',g.Handles.configiconholder,g.Icons.config.contrastbrightness,'Manually change contrast and brightness',@() updatetool('contrast'));
    txmapiconcell{1,end+1}=myicon('autocontrast',g.Handles.configiconholder,g.Icons.config.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    txmapiconcell{1,end+1}=myicon('resetlight',g.Handles.configiconholder,g.Icons.config.resetlight,'Reset contrast and brightness',@() segment('resetlight_Callback'),0);    
    txmapiconcell{1,end+1}=myicon('importfromother',g.Handles.configiconholder,g.Icons.config.importfromother,'Import LV segmentation from cine to Tx map image stack',@() segmentation('importfromcine2txmap_Callback'),0);
    txmapiconcell{1,end+1}=myicon('roipen',g.Handles.configiconholder,g.Icons.config.roipen,'ROI pen',@() updatetool('drawroi'));
    txmapiconcell{1,end+1}=myicon('addroiinlv',g.Handles.configiconholder,g.Icons.config.addroiinlv,'Add ROIs to sector of LV wall in selected slices',@() roi('roiaddinsector_Callback'),0);
    txmapiconcell{1,end+1}=myicon('T1',g.Handles.configiconholder,g.Icons.config.T1, 'T1 analysis',@() txmap('init',1),0);
    txmapiconcell{1,end+1}=myicon('T2',g.Handles.configiconholder,g.Icons.config.T2, 'T2 analysis',@() txmap('init',2),0);
    txmapiconcell{1,end+1}=myicon('T2star',g.Handles.configiconholder,g.Icons.config.T2star, 'T2* analysis',@() t2star.t2star,0);
    txmapiconcell{1,end+1}=myicon('ecv',g.Handles.configiconholder,g.Icons.config.ecv, 'ECV analysis',@() ecv('init_Callback'),0);
    txmapiconcell{1,end+1}=myicon('T1precolormap',g.Handles.configiconholder,g.Icons.config.T1precolormap, 'T1 pre colormap',@() tools('setcolormap_Callback','t1pre'),0);
    txmapiconcell{1,end+1}=myicon('T1postcolormap',g.Handles.configiconholder,g.Icons.config.T1postcolormap, 'T1 post colormap',@() tools('setcolormap_Callback','t1post'),0);
    txmapiconcell{1,end+1}=myicon('T2colormap',g.Handles.configiconholder,g.Icons.config.T2colormap, 'T2 colormap',@() tools('setcolormap_Callback','t2'),0);
    txmapiconcell{1,end+1}=myicon('T2starcolormap',g.Handles.configiconholder,g.Icons.config.T2starcolormap, 'T2* colormap',@() tools('setcolormap_Callback','t2star'),0);
    txmapiconcell{1,end+1}=myicon('ecvcolormap',g.Handles.configiconholder,g.Icons.config.ecvcolormap, 'ECV colormap',@() tools('setcolormap_Callback','ecv'),0);
    txmapiconcell{1,end+1}=myicon('greycolormap',g.Handles.configiconholder,g.Icons.config.greycolormap, 'Grey colormap',@() tools('setcolormap_Callback','gray'),0);
    txmapiconcell{1,end+1}= roiflowiconcell{1,end-2};
    txmapiconcell{1,end+1}=myicon('clearroi',g.Handles.configiconholder,g.Icons.config.clearroi,'Clear selected ROIs',@() roi('roidelete_Callback'),0);
    txmapiconcell{1,end+1}=myicon('clearallroi',g.Handles.configiconholder,g.Icons.config.clearallroi,'Clear all ROIs',@() roi('roiclearall_Callback'),0);

    g.Icons.txmapiconcell = txmapiconcell;
    end
    
     %----------------------------------------
    function initribbonplaceholder(varargin)
    %--------------------------------------

    g=varargin{1};
     %initiate iconholder for toggle axes 
    iconCell=cell(1,7);
    iconCell{1}=myicon('ribbonlv',g.Handles.toggleiconholder,g.Icons.toggleicons.ribbonlvoff,...
      'LV', @() segment('togglebuttonLV_Callback'),1,1,g.Icons.toggleicons.ribbonlvon);
    iconCell{2}=myicon('ribbonrv',g.Handles.toggleiconholder,g.Icons.toggleicons.ribbonrvoff,...
      'RV',@() segment('togglebuttonRV_Callback'),1,1,g.Icons.toggleicons.ribbonrvon);
    iconCell{3}=myicon('ribbonflow',g.Handles.toggleiconholder,g.Icons.toggleicons.ribbonflowoff,...
      'ROI/FLOW',@() segment('togglebuttonROIFLOW_Callback'),1,1,g.Icons.toggleicons.ribbonflowon);
        iconCell{4}=myicon('ribbonviability',g.Handles.toggleiconholder,g.Icons.toggleicons.ribbonviabilityoff,...
      'Viability',@() segment('togglebuttonVia_Callback'),1,1,g.Icons.toggleicons.ribbonviabilityon);
     iconCell{5}=myicon('ribbontxmap',g.Handles.toggleiconholder,g.Icons.toggleicons.ribbontxmapoff,...
      'TX maps and ECV',@() segment('togglebuttontxmap_Callback'),1,1,g.Icons.toggleicons.ribbontxmapon);
    iconCell{6}=myicon('ribbonanalysis',g.Handles.toggleiconholder,g.Icons.toggleicons.ribbonanalysisoff,...
      'Analysis',@() segment('togglebuttonAnalysis_Callback'),1,1,g.Icons.toggleicons.ribbonanalysison);
    iconCell{7}=myicon('ribbonimage',g.Handles.toggleiconholder,g.Icons.toggleicons.ribbonimageoff,...
      'Image',@() segment('togglebuttonImage_Callback'),1,1,g.Icons.toggleicons.ribbonimageon);
%   
     g.Handles.toggleiconholder.add(iconCell);   
    pos=plotboxpos(g.Handles.toggleiconholder.axeshandle);
    currentpos=get(g.Handles.toggleiconholder.axeshandle,'position');
    set(g.Handles.toggleiconholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);
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
    g.AxesTables.volume.addKey('LVCO','CO','l/min',bl);
    g.AxesTables.volume.addKey('LVHR','HR','bpm',bl);
    %g.AxesTables.volume.addKey('INFO','','','');
    
%     g.AxesTables.volume.addSpace();
    g.AxesTables.volume.addTable('RV',2,1,[0.6 0.4]);
    g.AxesTables.volume.addKey('RVM','ED/ES RVM','g',bl);
    g.AxesTables.volume.addKey('RVEDV','EDV','ml',bl);
    g.AxesTables.volume.addKey('RVESV','ESV','ml',bl);
    g.AxesTables.volume.addKey('RVSV','SV','ml',bl);
    g.AxesTables.volume.addKey('RVEF','EF','%',bl);
%     g.AxesTables.volume.addKey('RVCO','CO','l/min',bl);
%     g.AxesTables.volume.addKey('RVHR','HR','bpm',bl);
    
    %Flow report table
    g.AxesTables.flow.backgroundcolor = [0.94 0.94 0.94]; %[0 0 0];%
    g.AxesTables.flow.fontcolor = [0 0 0]; %[1 1 1];%
    g.AxesTables.flow.ystep = 20;%16;
    g.AxesTables.flow.fontsize = 10;
    g.AxesTables.flow.addTable('Flow',2,1,[0.65 0.35]);
    g.AxesTables.flow.addKey('ROI','ROI','','');
    g.AxesTables.flow.addKey('Netvol','Net vol','ml',bl);
   %g.AxesTables.flow.addKey('Forward','Forward','ml',bl);
    g.AxesTables.flow.addKey('Backward','Backward','ml',bl);
    g.AxesTables.flow.addKey('Regfrac','Regurg. frac.','%',bl);
    g.AxesTables.flow.addKey('FlowCO','FlowCO','l/min',bl);
    g.AxesTables.flow.addKey('FlowHR','HR','bpm',bl);
 
    %measurement report table
    g.AxesTables.measurement.backgroundcolor = [0.94 0.94 0.94]; %[0 0 0];%
    g.AxesTables.measurement.fontcolor = [0 0 0]; %[1 1 1];%
    g.AxesTables.measurement.ystep = 20;%16;
    g.AxesTables.measurement.fontsize = 10;
    g.AxesTables.measurement.addTable('Measurement',2,1,[0.6 0.4]);
    g.AxesTables.measurement.addKey('m1',bl,'',bl);
    g.AxesTables.measurement.addKey('m2',bl,'',bl);
    g.AxesTables.measurement.addKey('m3',bl,'',bl);
    g.AxesTables.measurement.addKey('m4',bl,'',bl);
            
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
    
    if isequal(36,getmodule(35,'KI',[],true));
      uimenu(g.Handles.analysismenu,'Label','RAMP', ...
        'Callback','plugin_RAMP(''init'')');
    end
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
    
    if ~isempty(SET) && ~isempty(SET(NO).Flow) && SET(NO).RoiN>0
      g.FlowNO = NO;
      g.FlowROI = SET(NO).RoiCurrent; %RoiN;
    elseif ~isempty(SET) && ~isempty(SET(NO).Flow) && ~isempty(SET(SET(NO).Flow.MagnitudeNo).Flow) && SET(SET(NO).Flow.MagnitudeNo).RoiN>0
      g.FlowNO = SET(NO).Flow.MagnitudeNo;
      g.FlowROI = SET(SET(NO).Flow.MagnitudeNo).RoiCurrent; %RoiN;
    else
      g.FlowNO = [];
      g.FlowROI = [];
      if strcmp(arg,'area')
        return;
      end
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
      case 'allreports'
        g.updatevolumeaxes;
        g.volumereportupdate(varargin{:});
        g.updateflowaxes;
        g.flowreportupdate;
        g.measurementreportupdate;        
    end
    
    end
    
    
    %-----------------------------
    function volumereportupdate(g)
    %-----------------------------
    %Update LV and RV report in axestables
    global SET DATA
    
%       m=g.();
      bla = {'---'};
      
      % -- LV --
      
      s = [];
      s.LVM = bla;
      s.LVEDV = bla;
      s.LVESV = bla;
      s.LVSV = bla;
      s.LVEF = bla;
      s.LVCO = bla;
      s.LVHR = bla;
      updatestruct = s;
      %       title{1} = 'LV';
      unitV='ml'; %default volume unit
      unitg='g';
      unitcorr=1; %correction for volume unit
      
      if ~isempty(DATA.LVNO) % isfield(DATA,'LVNO') && 
        no = DATA.LVNO(1);
        haslv = ~isempty(SET(no).EpiX) || ~isempty(SET(no).EndoX);
        
        if haslv          
          %Determine last character of key
          prepost = 'LV';
          smallanimal=false;
          if SET(no).EDV~=0 && ~isnan(SET(no).EDV)
            smallanimal = SET(no).EDV<=1;% && SET(no).EDV~=0;
          elseif SET(no).ESV~=0 && ~isnan(SET(no).ESV)
            smallanimal = SET(no).ESV<=1;
          end
          
          if smallanimal
            unitV=['\mu', 'l'];
            unitg='mg';
            unitcorr=1000;
          end
          if isnan(SET(no).LVM(SET(no).EDT)) && ~(SET(no).LVM(SET(no).EDT)>0)
            lvmedt = bla{1};
          else
            lvmedt = num2str(round(unitcorr*SET(no).LVM(SET(no).EDT)*1.05));
          end
          
          if isnan(SET(no).LVM(SET(no).EST))
            lvmest = bla{1};
          else
            lvmest = num2str(round(unitcorr*SET(no).LVM(SET(no).EST)*1.05));
          end
          updatestruct.([prepost 'M']) = sprintf('%s / %s',lvmedt,lvmest);
          updatestruct.([prepost 'EDV']) = round(unitcorr*SET(no).EDV);
          updatestruct.([prepost 'ESV']) = round(unitcorr*SET(no).ESV);
          updatestruct.([prepost 'SV']) = round(unitcorr*SET(no).SV);
          updatestruct.([prepost 'EF']) = round(100*SET(no).EF);
          updatestruct.([prepost 'CO']) = (SET(no).HeartRate*SET(no).SV)/1000;
          
          
          g.AxesTables.volume.updateUnit([prepost 'M'],unitg,true);
          g.AxesTables.volume.updateUnit([prepost 'EDV'],unitV,true);
          g.AxesTables.volume.updateUnit([prepost 'ESV'],unitV,true);
          g.AxesTables.volume.updateUnit([prepost 'SV'],unitV,true);

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
%       s.RVHR = bla;
      updatestruct = s;
      %       title{2} = 'RV';
      unitV='ml'; %default volume unit
      unitg='g';
      unitcorr=1; %correction for volume unit
        
      if ~isempty(DATA.RVNO) % isfield(DATA,'RVNO') && 
        no = DATA.RVNO;
        hasrv = ~isempty(SET(no).RVEpiX) || ~isempty(SET(no).RVEndoX);
        
        if hasrv          
          %Determine last character of key
          prepost = 'RV';   
          
%           smallanimal=SET(no).RVEDV<=1 && SET(no).RVEDV~=0;
          smallanimal=false;
          if SET(no).RVEDV~=0
            smallanimal = SET(no).RVEDV<=1;% && SET(no).EDV~=0;
          elseif SET(no).RVESV~=0
            smallanimal = SET(no).RVESV<=1;
          end
          
          if smallanimal
            unitV=['\mu' 'l'];
            unitg='mg';
            unitcorr=1000;
          end
          
          if isnan(SET(no).RVM(SET(no).EDT))
            rvmedt = bla{1};
          else
            rvmedt = num2str(round(unitcorr*SET(no).RVM(SET(no).EDT)*1.05));
          end
          
          if isnan(SET(no).RVM(SET(no).EST))
            rvmest = bla{1};
          else
            rvmest = num2str(round(unitcorr*SET(no).RVM(SET(no).EST)*1.05));
          end
          updatestruct.([prepost 'M']) = sprintf('%s / %s',rvmedt,rvmest);
          updatestruct.([prepost 'EDV']) = round(unitcorr*SET(no).RVEDV);
          updatestruct.([prepost 'ESV']) = round(unitcorr*SET(no).RVESV);
          updatestruct.([prepost 'SV']) = round(unitcorr*SET(no).RVSV);
          updatestruct.([prepost 'EF']) = round(100*SET(no).RVEF);
   
          g.AxesTables.volume.updateUnit([prepost 'M'],unitg,true);
          g.AxesTables.volume.updateUnit([prepost 'EDV'],unitV,true);
          g.AxesTables.volume.updateUnit([prepost 'ESV'],unitV,true);
          g.AxesTables.volume.updateUnit([prepost 'SV'],unitV,true);

          
%           title{2} = sprintf('RV  [Image Stack %d]',no);
        end
%         updatestruct.RVHR = round(SET(no).HeartRate);
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
    
    global SET DATA NO
    
    updatestruct = struct(...
      'ROI', '',...
      'Netvol', '---',... % 'Forward', '---',...
      'Backward', '---',...
      'Regfrac', '---',...
      'FlowCO', '---' ,...
      'FlowHR', '---');
    name.ROI = 'ROI';
    %set(DATA.Handles.flowstackpushbutton,'String','Set stack');
%     title{1} = 'Flow';
    
    if ~isempty(DATA.FlowNO) && ~isempty(DATA.FlowROI)
      no = DATA.FlowNO;
      roinbr = DATA.FlowROI;
      %Condition for non-time resolved Phase Contrast Flow measurement
      if length(SET(no).Roi(roinbr).T) == 1 && ~isempty(SET) && ~isempty(no) && SET(no).RoiN > 0 && ~isempty(SET(no).Flow.Result) && not(isempty(roinbr)) && length(SET(no).Flow.Result)>=roinbr && isfield(SET(no).Flow.Result(roinbr),'nettotvol') %For no time resolved Flow data
        updatestruct.Netvol = SET(no).Flow.Result(roinbr).nettotvol;
        updatestruct.Backward = SET(no).Flow.Result(roinbr).negflow;
        updatestruct.Regfrac ='---';
        updatestruct.FlowCO='---';
        updatestruct.FlowHR='---';
         name.ROI = SET(no).Roi(roinbr).Name;
         
         %Update names in tables structures and units
          g.AxesTables.flow.updateName('Netvol','NetFlow',true);
          g.AxesTables.flow.updateUnit('Netvol','ml/s',true);
          g.AxesTables.flow.updateName('Backward','NegFlow',true);
        g.AxesTables.flow.updateUnit('Backward','ml/s',true);
        
      elseif ~isempty(SET) && ~isempty(no) && SET(no).RoiN > 0 && ~isempty(SET(no).Flow.Result) && not(isempty(roinbr)) && length(SET(no).Flow.Result)>=roinbr && isfield(SET(no).Flow.Result(roinbr),'nettotvol')
        
        updatestruct.ROI = '';
        updatestruct.Netvol = SET(no).Flow.Result(roinbr).nettotvol;
%         updatestruct.Forward = SET(no).Flow.Result(roinbr).netforwardvol;
        updatestruct.Backward = SET(no).Flow.Result(roinbr).netbackwardvol;
        updatestruct.Regfrac = SET(no).Flow.Result(roinbr).regfrac;
        %define heart rate
        if not(isfield(SET(no).Flow,'HeartRate')) || isempty(SET(no).Flow.HeartRate)
          temphr = (60/(SET(no).TSize*SET(no).TIncr));
          if temphr < 35
            reportflow('defineheartrate',no);
          else
            SET(no).Flow.HeartRate = temphr;
          end
        end
        updatestruct.FlowCO = SET(no).Flow.Result(roinbr).nettotvol*SET(no).Flow.HeartRate/1000; 
        updatestruct.FlowHR = SET(no).Flow.HeartRate; %round((60/(SET(no).TSize*SET(no).TIncr))); %SET(no).HeartRate;
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
    for floop = 2:numel(fnames)
      fname = fnames{floop};
      g.AxesTables.flow.updateKey(fname,updatestruct.(fname),true);
    end
%     g.AxesTables.flow.updateTitle(title);
   g.AxesTables.flow.draw();
   %set(g.Handles.flowuipanel,'Visible','on');
            
    end
    
    
    %---------------------------------
    function measurementreportupdate(g)
    %---------------------------------
    %Update scar report in axestables
    
    global SET NO
                
      no = NO;
      bla = {'---'};
      hasscar = ~isempty(SET(no).Scar) && ~isempty(SET(no).Scar.Percentage);
      hasmar = ~isempty(SET(no).MaR) && ~isempty(SET(no).MaR.Result);
      hasflow = ~isempty(SET) && SET(no).RoiN > 0 && isfield(SET(no).Flow,'nettotvol');
      hasmmode = (strcmp(g.ViewPanelsType{g.CurrentPanel},'mmodespatial') ...
        || strcmp(g.ViewPanelsType{g.CurrentPanel},'mmodetemporal'));
      hasatrialscar = isfield(SET(no),'AtrialScar') && ~isempty(SET(no).AtrialScar) && ~isempty(SET(no).AtrialScar.Percentage);
      updatestruct = [];
      
      if hasscar
        scarpro2mlcoeff=sum(SET(no).Scar.MyocardMask(:))/100/1000*SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceGap+SET(no).SliceThickness);
        scarml=round(10*SET(no).Scar.Percentage*scarpro2mlcoeff)/10;     
      end
        
      if hasmar && hasscar
        
        updatestruct.m1 = scarml;%round(SET(no).Scar.Percentage);%([SET(no).EDT, SET(no).EST]));
        updatestruct.m2 = round(10*SET(no).Scar.Percentage)/10;%round(SET(no).Scar.Percentage);%([SET(no).EDT, SET(no).EST]));
        updatestruct.m3 = round(10*SET(no).Scar.MOPercentage)/10;
        updatestruct.m4 = [num2str(round(SET(no).MaR.Percentage(SET(no).EDT))),'/',num2str(round(SET(no).MaR.Percentage(SET(no).EST)))];
        
        name.m1 = 'Scar';
        unit.m1 = 'ml';
        name.m2 = 'Scar';
        unit.m2 = '%';
        name.m3 = 'MO';
        unit.m3 = '%';
        name.m4 = 'ED/ES MaR';
        unit.m4 = '%';
        title{1} = sprintf('Scar & MaR');
           
      elseif hasscar
        %updatestruct.m1 = round(SET(no).LVM(SET(no).EDT)*1.05);
        updatestruct.m1 = scarml;%round(SET(no).Scar.Percentage);%([SET(no).EDT, SET(no).EST]));
        updatestruct.m2 = round(10*SET(no).Scar.Percentage)/10;
        updatestruct.m3 = round(10*SET(no).Scar.MOPercentage)/10;
        
        %name.m1 = 'LVM';
        name.m1 = 'Scar';
        %unit.m1 = 'g';
        unit.m1 = 'ml';
        
        name.m2 = 'Scar';
        unit.m2 = '%';
        
        name.m3 = 'MO';
        unit.m3 = '%';
        
        updatestruct.m4 = bla;
        name.m4 = bla;
        unit.m4 = '';
        
        title{1} = sprintf('Scar');
        
      elseif hasmar
        %updatestruct.m1 = round(SET(no).LVM(SET(no).EDT)*1.05);
        updatestruct.m1 = [num2str(round(SET(no).MaR.Percentage(SET(no).EDT))),'/',num2str(round(SET(no).MaR.Percentage(SET(no).EST)))];%round(SET(no).MaR.Percentage([SET(no).EDT,SET(no).EST]));
        %name.m1 = 'LVM';
        name.m1 = 'ED/ES MaR';
        %unit.m1 = 'g';
        unit.m1 = '%';
        
        updatestruct.m2 = bla;
        name.m2 = bla;
        unit.m2 = '';
        
        updatestruct.m3 = bla;
        name.m3 = bla;
        unit.m3 = '';
         
        updatestruct.m4 = bla;
        name.m4 = bla;
        unit.m4 = '';
        title{1} = sprintf('MaR');  
        
      elseif hasatrialscar
        updatestruct.m1 = round(10*SET(no).AtrialScar.Percentage)/10;
        
        name.m1 = 'Atrial Scar';
        unit.m1 = '%';
                
        updatestruct.m2 = bla;
        name.m2 = bla;
        unit.m2 = '';
        
        updatestruct.m3 = bla;
        name.m3 = bla;
        unit.m3 = '';
        
        updatestruct.m4 = bla;
        name.m4 = bla;
        unit.m4 = '';
        
        title{1} = sprintf('Atrial Scar');
        
%       elseif hasflow
%         rloop = 1;
%         updatestruct.m1 = SET(no).Flow.nettotvol(rloop);
%         updatestruct.m2 = SET(no).Flow.regfrac(rloop);
%         name.m1 = 'Net vol';
%         name.m2 = 'Regurg. frac.';
%         unit.m1 = 'ml';
%         unit.m2 = '%';
%         title{1} = sprintf('Flow');
      elseif hasmmode && ~hasflow
        [dist,timedist] = calcfunctions('calcmmodedists',no);
        updatestruct.m1 = dist;
        updatestruct.m2 = timedist;
        name.m1 = 'Distance';
        name.m2 = 'Time';
        unit.m1 = 'mm';
        unit.m2 = 'ms';
        updatestruct.m3 = bla;
        name.m3 = bla;
        unit.m3 = '';
        
        updatestruct.m4 = bla;
        name.m4 = bla;
        unit.m4 = '';
        title{1} = sprintf('Mmode: Distance between lines');
      else
        %bla = {'---'};
        updatestruct.m1 = bla;
        updatestruct.m2 = bla;
        updatestruct.m3 = bla;
        updatestruct.m4 = bla;
        name.m1 = bla;
        name.m2 = bla;
        name.m3 = bla;
        name.m4 = bla;
        unit.m1 = '';
        unit.m2 = '';
        unit.m3 = '';
        unit.m4 = '';
        title{1} = '';
      end
      

      if hasscar || hasmar || hasatrialscar || (hasmmode && ~hasflow) 
        set(g.Handles.measurementuipanel,'Visible','on');
        %set(g.Handles.measurementslideruipanel,'Visible','on');
        set(g.Handles.flowuipanel,'Visible','off');
        %set(g.Handles.flowaxes,'Visible','off');
        %set(g.Handles.flowresultaxes,'Visible','off');
        %set(g.Handles.measurementresultaxes,'Visible','on');
        g.AxesTables.measurement.show;
        %viability('sliderupdate')
%         if hasscar
%           set([g.Handles.slider1,g.Handles.slider2, g.Handles.slider1edit,g.Handles.slider2edit],'Visible', 'on')
%         else
%           set([g.Handles.slider1,g.Handles.slider2, g.Handles.slider1edit,g.Handles.slider2edit],'Visible', 'off')
%         end
      %else
        %   viability('sliderupdate')
      end
      viability('sliderupdate')
      
      fnames = fieldnames(updatestruct);
      for floop = 1:numel(fnames)
        fname = fnames{floop};
        g.AxesTables.measurement.updateKey(fname,updatestruct.(fname),true);
        g.AxesTables.measurement.updateName(fname,name.(fname),true);
        g.AxesTables.measurement.updateUnit(fname,unit.(fname),true);
      end
      g.AxesTables.measurement.updateTitle(title);
      g.AxesTables.measurement.draw();
      
      if not(hasscar) && not(hasmar) && not(hasatrialscar) && not(hasflow) && not(hasmmode)
        g.AxesTables.measurement.hide;
        set(g.Handles.measurementuipanel,'Visible','off');
%        val=get(g.Handles.hideallpanelscheckbox,'value');
%        if val==0
%          set(g.Handles.flowuipanel,'Visible','on');
 %       end
%       else        
%         xval = get(g.Handles.measurementresultaxes,'xlim');
%         yval = get(g.Handles.measurementresultaxes,'ylim');
%         hold(g.Handles.measurementresultaxes,'on');
%         plot(g.Handles.measurementresultaxes,[xval(1) xval(2) xval(2) xval(1) xval(1)],[yval(1) yval(1) yval(2) yval(2) yval(1)],'k-');
%         hold(g.Handles.measurementresultaxes,'off');
      end
      
    end
    
    %---------------------------------
    function measurementreportclear(g)
    %---------------------------------
    %Clear measurement report in axestables
    
    updatestruct = [];
    bla = {'---'};
    updatestruct.m1 = bla;
    updatestruct.m2 = bla;
    name.m1 = bla;
    name.m2 = bla;
    unit.m1 = '';
    unit.m2 = '';
    title{1} = '';
    
    fnames = fieldnames(updatestruct);
    for floop = 1:numel(fnames)
      fname = fnames{floop};
      g.AxesTables.measurement.updateKey(fname,updatestruct.(fname),true);
      g.AxesTables.measurement.updateName(fname,name.(fname),true);
      g.AxesTables.measurement.updateUnit(fname,unit.(fname),true);
    end
    g.AxesTables.measurement.updateTitle(title);
    g.AxesTables.measurement.draw();
    
    g.AxesTables.measurement.hide;
    set(g.Handles.measurementuipanel,'Visible','off');
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
    
    %-----------------------
    function initbullseye(g)
    %-----------------------
    %init the listbox in the bullseye reporter
    end
    
    %-----------------------
    function initbullseyeslices(g) %#ok<MANU>
    %-----------------------
    %Autoselects slices with myocardium in them before GUI is opened
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