classdef segmentgui < maingui
%Contains methods and properties that are specific to (original) Segment.
 
properties
    AxesTables = [];
    LVNO = [];
    RVNO = [];
    FlowNO = [];
    FlowROI = [];
    fontsizeincm = [];
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
    g = varargin{1};
    
    gHnd = g.Handles;
    gIcn = g.Icons.config;
    iconcell = cell(1,1);
    iconcell{1,1} = myicon('database',gHnd.permanenticonholder,gIcn.database,'Open from disc [Ctrl-O]',@() segment('fileopen_Callback'),0);
    iconcell{1,end+1} = myicon('databaseadd',gHnd.permanenticonholder,gIcn.databaseadd,'Save to disc [Ctrl-S]',... 
      @() filemenu('saveall_Callback'),0);
    iconcell{1,end+1} = myicon('closeall',gHnd.permanenticonholder,gIcn.closeall,'Close all image stacks [Ctrl-Shift-W]',@() segment('filecloseall_Callback'),0);
    
    iconcell{1,end+1} = myicon('panel1',gHnd.permanenticonholder,gIcn.panel1,'View one image panel [Shift-1]',@() viewfunctions('setview',1,1),1,1);
    iconcell{1,end+1} = myicon('panel2',gHnd.permanenticonholder,gIcn.panel2,'View two image panels [Shift-2]',@() viewfunctions('setview',1,2),1,1);
    iconcell{1,end+1} = myicon('panel2x1',gHnd.permanenticonholder,gIcn.panel2x1,'View two image panels [Alt-2]',@() viewfunctions('setview',2,1),1,1);
    iconcell{1,end+1} = myicon('panel3x1',gHnd.permanenticonholder,gIcn.panel3x1,'View three image panels [Alt-3]',@() viewfunctions('setview',3,1),1,1);
    iconcell{1,end+1} = myicon('panel3',gHnd.permanenticonholder,gIcn.panel3,'View three image panels [Shift-3]',@() viewfunctions('setview',1,3),1,1);
    iconcell{1,end+1} = myicon('panel4',gHnd.permanenticonholder,gIcn.panel4,'View four image panels [Shift-4]',@() viewfunctions('setview',2,2),1,1);
    iconcell{1,end+1} = myicon('panel6',gHnd.permanenticonholder,gIcn.panel6,'View six image panels [Shift-6]',@()  viewfunctions('setview',2,3),1,1);
    iconcell{1,end+1} = myicon('orthoview',gHnd.permanenticonholder,gIcn.orthoview,'Orthogonal view',@() callbackfunctions('orthoview_Callback'),1,1);
    iconcell{1,end+1} = myicon('saveview',gHnd.permanenticonholder,gIcn.saveview,'Save view',@() segmentview,0);
    
    iconcell{1,end+1} = myicon('viewone',gHnd.permanenticonholder,gIcn.viewone,'View one slice [Ctrl-1]',@() viewfunctions('setviewtype','one'),1,2);
    iconcell{1,end+1} = myicon('viewall',gHnd.permanenticonholder,gIcn.viewall,'View all slices [Ctrl-3]',@() viewfunctions('setviewtype','montage'),1,2);
    iconcell{1,end+1} = myicon('viewrow',gHnd.permanenticonholder,gIcn.viewrow,'View all slices in 2 rows [Ctrl-4]',@() viewfunctions('setviewtype','montagerow'),1,2);
    
    iconcell{1,end+1} = myicon('undo',gHnd.permanenticonholder,gIcn.undo,'Undo last operation [Ctrl-Z]',@() tools('undosegmentation_Callback'),0);
    iconcell{1,end+1} = myicon('refresh',gHnd.permanenticonholder,gIcn.refresh,'Refresh image view',@() viewfunctions('setview'),0);
    
    iconcell{1,end+1} = myicon('play',gHnd.permanenticonholder,gIcn.play,'Play movie [P]',@() callbackfunctions('play_Callback'),2);
    iconcell{1,end+1} = myicon('next',gHnd.permanenticonholder,gIcn.next,'Next frame [Rigth arrow]',@() viewfunctions('switchtimeframe',1,true),0);
    iconcell{1,end+1} = myicon('prev',gHnd.permanenticonholder,gIcn.prev,'Previous frame [Left arrow]',@() viewfunctions('switchtimeframe',-1,true),0);
    iconcell{1,end+1} = myicon('faster',gHnd.permanenticonholder,gIcn.faster,'Faster frame rate',@() segment('fasterframerate_Callback'),0);
    iconcell{1,end+1} = myicon('slower',gHnd.permanenticonholder,gIcn.slower,'Slower frame rate',@() segment('slowerframerate_Callback'),0);
    
    iconcell{1,end+1} = myicon('hideall',gHnd.permanenticonholder,gIcn.hideall,'Hide/show all overlays [H]',@() viewfunctions('viewhideall_Callback'),2);
    iconcell{1,end+1} = myicon('clearall',gHnd.permanenticonholder,gIcn.clearall,'Clear all segmentation in current image stack',@() callbackfunctions('segmentclearall_Callback'),0);
    iconcell{1,end+1} = myicon('clearalledes',gHnd.permanenticonholder,gIcn.clearalledes,'Clear all segmentation except in ED and ES',@() callbackfunctions('segmentclearallbutsystolediastole_Callback'),0); 
    iconcell{1,end+1} = myicon('zoomin',gHnd.permanenticonholder,gIcn.zoomin,'Zoom in [Ctrl-plus]',@() viewfunctions('zoom',1),0);%segment('viewzoomin_Callback'),0);
    iconcell{1,end+1} = myicon('zoomout',gHnd.permanenticonholder,gIcn.zoomout,'Zoom out [Ctrl-minus]',@() viewfunctions('zoom',-1),0);%segment('viewzoomout_Callback'),0);
    iconcell{1,end+1} = myicon('autozoom',gHnd.permanenticonholder,gIcn.autozoom,'Auto zoom',@() segment('autozoom'),0);
    iconcell{1,end+1} = myicon('colorbar',gHnd.permanenticonholder,gIcn.colorbar,'Hide colorbar',@() viewfunctions('viewhidecolorbar_Callback'),2);
    iconcell{1,end+1} = myicon('viewpixels',gHnd.permanenticonholder,gIcn.viewpixels,'Show image pixels',@() viewfunctions('viewinterp_Callback'),2);
%     iconcell{1,end+1} = myicon('reportsheet',gHnd.permanenticonholder,gIcn.reportsheet,'Open Report sheet generation',@() reporter.reportsheet,0);
    iconcell{1,end+1} = myicon('savescreen',gHnd.permanenticonholder,gIcn.savescreen,'Save screen shot',@() export('screenshot_Callback'),0);
    
    iconcell{1,end+1} = myicon('settingsgeneral',gHnd.permanenticonholder,gIcn.settingsgeneral,'Set general preferences',@() segpref,0);
    iconcell{1,end+1} = myicon('settingssystem',gHnd.permanenticonholder,gIcn.settingsdatabase,'Set patient database preferences',@() segpref('advancedsettings_Callback'),0);
    iconcell{1,end+1} = myicon('settingspacs',gHnd.permanenticonholder,gIcn.settingspacs,'Set PACS connection preferences',@() pacspref,0);
    gHnd.permanenticonholder.add(iconcell);
    
    pos = plotboxpos(gHnd.permanenticonholder.axeshandle);
    currentpos=get(gHnd.permanenticonholder.axeshandle,'position');
    set(gHnd.permanenticonholder.axeshandle,'position',currentpos-[pos(1),0,0,0]);
    set(gHnd.iconuipanel,'visible','on')
    
    end
    
    %---------------------------------------------
    function initconfigplaceholder(varargin)
    %--------------------------------------------
    g = varargin{1};
    gHnd = g.Handles;
    gIcn = g.Icons.config;
      
    %initcells
    lviconcell = cell(1,1);
    rviconcell = cell(1,1);
    strainiconcell = cell(1,1);
    analysisiconcell = cell(1,1);
    roiflowiconcell = cell(1,1);
    viabilityiconcell = cell(1,1);
    imageiconcell = cell(1,1);
    %LV
    lviconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'));
    lviconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'));
    lviconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'));
    lviconcell{1,end+1} = myicon('contrastbrightness',gHnd.configiconholder,gIcn.contrastbrightness,'Manually change contrast and brightness',@() buttondownfunctions('updatebuttondowns','Contrast'));
    lviconcell{1,end+1} = myicon('autocontrast',gHnd.configiconholder,gIcn.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    lviconcell{1,end+1} = myicon('lvstack',gHnd.configiconholder,gIcn.lvstack,'Go to LV stack',@() segment('viewspecial_Callback','lv'),0);
    lviconcell{1,end+1} = myicon('selectoneall',gHnd.configiconholder,gIcn.selectoneall,'Select one or all frames mode [1]',@()segment_main('singleframemode_Callback'),2,1,gIcn.selectallone);
    
    lviconcell{1,end+1} = myicon('moveall',gHnd.configiconholder,gIcn.moveall,'Translate all contours',@() buttondownfunctions('updatebuttondowns','moveall'));    
    lviconcell{1,end+1} = myicon('autolv',gHnd.configiconholder,gIcn.autolv,'AI-based semi-automatic LV segmentation',@() lvsegmentationml,0);
    % check if GPU available   
    if gpuDeviceCount > 0
      lviconcell{1,end+1} = myicon('autolvwand',gHnd.configiconholder,gIcn.autolvwand,'AI-based automatic LV segmentation [Ctrl-L]',@() lvsegmentationml('doaltlv2_Callback'),0);
    else
      lviconcell{1,end+1} = myicon('autolvedes',gHnd.configiconholder,gIcn.autolvedes,'AI-based automatic LV segmentation for ED/ES [Ctrl-L]',@() lvsegmentationml('doaltlv2_Callback',true),0);
    end
    lviconcell{1,end+1} = myicon('endopen',gHnd.configiconholder,gIcn.endopen,'Endo pen',@() buttondownfunctions('updatebuttondowns','Endo'));
    lviconcell{1,end+1} = myicon('epipen',gHnd.configiconholder,gIcn.epipen,'Epi pen',@() buttondownfunctions('updatebuttondowns','Epi'));
    lviconcell{1,end+1} = myicon('interpendo',gHnd.configiconholder,gIcn.interpendo,'Set interpolation points for Endo',@() buttondownfunctions('updatebuttondowns','EndoInterp'));
    lviconcell{1,end+1} = myicon('interpepi',gHnd.configiconholder,gIcn.interpepi,'Set interpolation points for Epi',@() buttondownfunctions('updatebuttondowns','EpiInterp'));
    lviconcell{1,end+1} = myicon('balloonendo',gHnd.configiconholder,gIcn.balloonendo,'Semi automatic Endo tool',@() buttondownfunctions('updatebuttondowns','EndoBalloon'));
%     lviconcell{1,end+1} = myicon('balloonepi',gHnd.configiconholder,gIcn.balloonepi,'Semi automatic Epi tool',@() lvsegmentation('smartepi_Callback'),0);
    lviconcell{1,end+1} = myicon('smooth',gHnd.configiconholder,gIcn.smooth,'Smooth current contour [O]',@() tools('smoothsegmentation_Callback'),0);
          
    lviconcell{1,end+1} = myicon('interpseginslice',gHnd.configiconholder,gIcn.interpseginslice,'Interpolate segmentation over slices',@() lv('interpolatedelineation_Callback'),0);
    lviconcell{1,end+1} = myicon('interpsegintime',gHnd.configiconholder,gIcn.interpsegintime,'Interpolate segmentation in time for selected slices',@() segmentation('interpolatedelineationovertime_Callback'),0);
    lviconcell{1,end+1} = myicon('refineendo',gHnd.configiconholder,gIcn.refineendo,'Refine LV Endo [Ctrl-R]',@() lvpeter('segmentrefineendo_Callback'),0);
    lviconcell{1,end+1} = myicon('refineepi',gHnd.configiconholder,gIcn.refineepi,'Refine LV Epi [Ctrl-Shift-R]',@() lvpeter('segmentrefineepi_Callback'),0);
    lviconcell{1,end+1} = myicon('propagateendo',gHnd.configiconholder,gIcn.propagateendo,'Propagate LV Endo forward and refine [Ctrl-F]', @() lvpeter('segmentpropagateendo_Callback'),0);
    lviconcell{1,end+1} = myicon('propagateepi',gHnd.configiconholder,gIcn.propagateepi,'Propagate LV Epi forward and refine [Ctrl-Shift-F]',@() lvpeter('segmentpropagateepi_Callback'),0);
    
    lviconcell{1,end+1} = myicon('contractendo',gHnd.configiconholder,gIcn.contractendo,'Contract LV Endo [Ctrl-K]',@() lv('segmentexpandcontract_Callback',-1,'endo'),0);
    lviconcell{1,end+1} = myicon('expandendo',gHnd.configiconholder,gIcn.expandendo,'Expand LV Endo [Ctrl-E]',@() lv('segmentexpandcontract_Callback',1,'endo'),0);
    lviconcell{1,end+1} = myicon('contractepi',gHnd.configiconholder,gIcn.contractepi,'Contract LV Epi [Ctrl-Alt-K]',@() lv('segmentexpandcontract_Callback',-1,'epi'),0);
    lviconcell{1,end+1} = myicon('expandepi',gHnd.configiconholder,gIcn.expandepi,'Expand LV Epi [Ctrl-Alt-E]',@() lv('segmentexpandcontract_Callback',1,'epi'),0);
    lviconcell{1,end+1} = myicon('evenoutwall',gHnd.configiconholder,gIcn.evenoutwall,'Even out myocardium wall [Alt-Up]',@() segment('smoothendowall_Callback'),0);
    lviconcell{1,end+1} = myicon('copylvup',gHnd.configiconholder,gIcn.copylvup,'Copy LV upwards and refine [Ctrl-U]',@()tools('copyupward_Callback'),0);
    lviconcell{1,end+1} = myicon('copylvdown',gHnd.configiconholder,gIcn.copylvdown,'Copy LV downwards and refine [Ctrl-D]',@()tools('copydownward_Callback'),0);
        
    lviconcell{1,end+1} = myicon('hidelv',gHnd.configiconholder,gIcn.hidelv,'Hide LV segmentation',@() viewfunctions('hide_Callback'),2);
    lviconcell{1,end+1} = myicon('hideinterp',gHnd.configiconholder,gIcn.hideinterp,'Hide interpolation points',@() viewfunctions('hide_Callback'),2);
    lviconcell{1,end+1} = myicon('clearalllv',gHnd.configiconholder,gIcn.clearalllv,'Clear all LV segmentation',@() callbackfunctions('segmentclearalllv_Callback'),0);
    lviconcell{1,end+1} = myicon('clearendo',gHnd.configiconholder,gIcn.clearendo,'Clear LV Endo in selected slices according to mode',@() segmentation('clearslicesthis_Callback',1,0,0,0),0);
    lviconcell{1,end+1} = myicon('clearepi',gHnd.configiconholder,gIcn.clearepi,'Clear LV Epi in selected slices according to mode',@() segmentation('clearslicesthis_Callback',0,1,0,0),0);
    lviconcell{1,end+1} = myicon('volumecurve',gHnd.configiconholder,gIcn.volumecurve,'Plot volume curve',@() lvpeter('plotvolumecurve'),0);
   
    g.Icons.lviconcell = lviconcell;
    
    %RV
    rviconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'));
    rviconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour', @() buttondownfunctions('updatebuttondowns','move'));
    rviconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'));
    rviconcell{1,end+1} = myicon('contrastbrightness',gHnd.configiconholder,gIcn.contrastbrightness,'Manually change contrast and brightness',@() buttondownfunctions('updatebuttondowns','Contrast'));
    rviconcell{1,end+1} = myicon('autocontrast',gHnd.configiconholder,gIcn.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    rviconcell{1,end+1} = myicon('rvstack',gHnd.configiconholder,gIcn.rvstack,'Go to RV stack',@() segment('viewspecial_Callback','rv'),0);
    rviconcell{1,end+1} = lviconcell{1,7};

    rviconcell{1,end+1} = lviconcell{1,8}; %moveall
    rviconcell{1,end+1} = myicon('autorvendo',gHnd.configiconholder,gIcn.autorvendo,'Automatic RV Endo segmentation [Ctrl-Alt-M]',@() rvsegmentation,0);
    rviconcell{1,end+1} = myicon('rvendopen',gHnd.configiconholder,gIcn.rvendopen,'RV Endo pen',@() buttondownfunctions('updatebuttondowns','RVEndo'));
    rviconcell{1,end+1} = myicon('rvepipen',gHnd.configiconholder,gIcn.rvepipen,'RV Epi pen',@() buttondownfunctions('updatebuttondowns','RVEpi'));
    rviconcell{1,end+1} = myicon('interprvendo',gHnd.configiconholder,gIcn.interprvendo,'Set interpolation points for RV Endo',@() buttondownfunctions('updatebuttondowns','RVEndoInterp'));
    rviconcell{1,end+1} = myicon('interprvepi',gHnd.configiconholder,gIcn.interprvepi,'Set interpolation points for RV Epi',@() buttondownfunctions('updatebuttondowns','RVEpiInterp'));
    rviconcell{1,end+1} = myicon('balloonrvendo',gHnd.configiconholder,gIcn.balloonrvendo,'Semi automatic RV Endo tool',@() buttondownfunctions('updatebuttondowns','RVEndoBalloon'));
    rviconcell{1,end+1} = myicon('refinervendo',gHnd.configiconholder,gIcn.refinervendo,'Refine RV Endo [Ctrl-Alt-R]',@() rv('segmentrefinervendo_Callback'),0);
    %need icon
    rviconcell{1,end+1} = myicon('copyrvup',gHnd.configiconholder,gIcn.copyrvup,'Copy RV Endo upwards and refine [Ctrl-Alt-U]',@()tools('copyupward_Callback','rvendo',false,false),0);
    rviconcell{1,end+1} = myicon('copyrvdown',gHnd.configiconholder,gIcn.copyrvdown,'Copy RV Endo downwards and refine [Ctrl-Alt-D]',@()tools('copydownward_Callback','rvendo',false,false),0);
    
    rviconcell{1,end+1} = myicon('interpsegintime',gHnd.configiconholder,gIcn.interpsegintime,'Interpolate segmentation in time',@() segmentation('interpolatedelineationovertime_Callback'),0);
    rviconcell{1,end+1} = myicon('interpseginslice',gHnd.configiconholder,gIcn.interpseginslice,'Interpolate segmentation over slices',@() lv('interpolatedelineation_Callback'),0);
    
    rviconcell{1,end+1} = myicon('hiderv',gHnd.configiconholder,gIcn.hiderv,'Hide RV segmentation', @() viewfunctions('hide_Callback'),2);
    rviconcell{1,end+1} = lviconcell{end-4};%myicon('hideinterp',gHnd.configiconholder,gIcn.hideinterp,'Hide interpolation points',@() segment('viewhideinterp_Callback'),2);
    rviconcell{1,end+1} = myicon('clearallrv',gHnd.configiconholder,gIcn.clearallrv,'Clear all RV segmentation',@() callbackfunctions('segmentclearallrv_Callback'),0);
    rviconcell{1,end+1} = myicon('clearrvendo',gHnd.configiconholder,gIcn.clearrv,'Clear RV endo in selected slices according to mode',@() segmentation('clearslicesthis_Callback',0,0,1,0),0);
    rviconcell{1,end+1} = myicon('clearrvepi',gHnd.configiconholder,gIcn.clearrvepi,'Clear RV epi in selected slices according to mode',@() segmentation('clearslicesthis_Callback',0,0,0,1),0);
    rviconcell{1,end+1} = myicon('volumecurve',gHnd.configiconholder,gIcn.volumecurve,'Plot volume curve',@() lvpeter('plotvolumecurve'),0);
    
    g.Icons.rviconcell = rviconcell;
    
    %ROIFLOW
    
    roiflowiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'));
    roiflowiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'));
    roiflowiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'));
    %roiflowiconcell{1,end+1} = myicon('scaleROI',gHnd.configiconholder,gIcn.scaleROI,'Scale ROI',@() updatetool('scaleROI'));
    roiflowiconcell{1,end+1} = myicon('contrastbrightness',gHnd.configiconholder,gIcn.contrastbrightness,'Manually change contrast and brightness',@() buttondownfunctions('updatebuttondowns','Contrast'));
    roiflowiconcell{1,end+1} = myicon('autocontrast',gHnd.configiconholder,gIcn.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    roiflowiconcell{1,end+1} = myicon('flowstack',gHnd.configiconholder,gIcn.flowstack,'Go to flow stack',@() segment('viewspecial_Callback','flow'),0);
    roiflowiconcell{1,end+1} = lviconcell{1,7};

    roiflowiconcell{1,end+1} = myicon('putroi',gHnd.configiconholder,gIcn.putroi,'Place ROI',@() buttondownfunctions('updatebuttondowns','RoiPut'));
    roiflowiconcell{1,end+1} = myicon('roipen',gHnd.configiconholder,gIcn.roipen,'ROI pen',@() buttondownfunctions('updatebuttondowns','Roi'));
    roiflowiconcell{1,end+1} = myicon('balloonroi',gHnd.configiconholder,gIcn.balloonroi,'Semi automatic ROI tool',@() buttondownfunctions('updatebuttondowns','RoiBalloon'));
    roiflowiconcell{1,end+1} = myicon('trackingvessel',gHnd.configiconholder,gIcn.trackingvessel,'Track vessel over time [Alt-T]',@() vesselsnake_flowtrackroi('flowtrackroi'),0);
    roiflowiconcell{1,end+1} = myicon('refineroi',gHnd.configiconholder,gIcn.refineroi,'Refine Flow ROI [Alt-R]',@() flow('flowrefine_Callback'),0);
    roiflowiconcell{1,end+1} = myicon('refineroinext',gHnd.configiconholder,gIcn.refineroinext,'Propagate Flow ROI forward and refine [Alt-F]',@() flow('flowpropagate_Callback'),0);
    roiflowiconcell{1,end+1} = myicon('contractroi',gHnd.configiconholder,gIcn.contractroi,'Contract ROI',@() roi('expandcontract_Callback',-1),0);
    roiflowiconcell{1,end+1} = myicon('expandroi',gHnd.configiconholder,gIcn.expandroi,'Expand ROI',@() roi('expandcontract_Callback',1),0);
    roiflowiconcell{1,end+1} = myicon('unwrap',gHnd.configiconholder,gIcn.unwrap,'Unwrap flow',@() flowunwrap,0);
    roiflowiconcell{1,end+1} = myicon('palette',gHnd.configiconholder,gIcn.palette,'Set ROI color',@() roi('roisetcolor_Callback'),0);
    roiflowiconcell{1,end+1} = myicon('text',gHnd.configiconholder,gIcn.text,'Set ROI label',@() roi('roisetlabel_Callback'),0);
    roiflowiconcell{1,end+1} = myicon('plotflow',gHnd.configiconholder,gIcn.plotflow,'Plot flow [Ctrl-T]',@() reportflow,0);
    roiflowiconcell{1,end+1} = myicon('hideroi',gHnd.configiconholder,gIcn.hideroi,'Hide ROI',@() viewfunctions('hide_Callback'),2);
    roiflowiconcell{1,end+1} = myicon('hidetext',gHnd.configiconholder,gIcn.hidetext,'Hide text',@() viewfunctions('hide_Callback'),2);
    roiflowiconcell{1,end+1} = myicon('clearroi',gHnd.configiconholder,gIcn.clearroi,'Clear selected ROIs',@() roi('roidelete_Callback'),0);  
    roiflowiconcell{1,end+1} = myicon('clearallroi',gHnd.configiconholder,gIcn.clearallroi,'Clear all ROIs',@() roi('roiclearall_Callback'),0); 
    g.Icons.roiflowiconcell = roiflowiconcell;
    
    %Viablility
    viabilityiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'));
    viabilityiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'));
    viabilityiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'));
    viabilityiconcell{1,end+1} = myicon('contrastbrightness',gHnd.configiconholder,gIcn.contrastbrightness,'Manually change contrast and brightness',@() buttondownfunctions('updatebuttondowns','Contrast'));
    viabilityiconcell{1,end+1} = myicon('autocontrast',gHnd.configiconholder,gIcn.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    viabilityiconcell{1,end+1} = myicon('scarstack',gHnd.configiconholder,gIcn.scarstack,'Go to scar stack',@() segment('viewspecial_Callback','cinescar'),0);
    viabilityiconcell{1,end+1} = lviconcell{1,7};

    viabilityiconcell{1,end+1} = myicon('importfromother',gHnd.configiconholder,gIcn.importfromother,'Import LV segmentation from cine to scar image stack',@() segmentation('importfromcine2scar_Callback'),0);
    viabilityiconcell{1,end+1} = myicon('moveall',gHnd.configiconholder,gIcn.moveall,'Translate all contours',@() buttondownfunctions('updatebuttondowns','moveall'));    
    viabilityiconcell{1,end+1} = myicon('endopen',gHnd.configiconholder,gIcn.endopen,'Endo pen',@() buttondownfunctions('updatebuttondowns','Endo'));
    viabilityiconcell{1,end+1} = myicon('epipen',gHnd.configiconholder,gIcn.epipen,'Epi pen',@() buttondownfunctions('updatebuttondowns','Epi'));
    viabilityiconcell{1,end+1} = myicon('interpendo',gHnd.configiconholder,gIcn.interpendo,'Set interpolation points for Endo',@() buttondownfunctions('updatebuttondowns','EndoInterp'));
    viabilityiconcell{1,end+1} = myicon('interpepi',gHnd.configiconholder,gIcn.interpepi,'Set interpolation points for Epi',@() buttondownfunctions('updatebuttondowns','EpiInterp'));
    viabilityiconcell{1,end+1} = myicon('autoscar',gHnd.configiconholder,gIcn.autoscar,'Auto scar',@() viability('viabilityautoewa'),0);
    viabilityiconcell{1,end+1} = myicon('scarpen',gHnd.configiconholder,gIcn.scarpen,'Draw scar',@() buttondownfunctions('updatebuttondowns','Scar'));
    viabilityiconcell{1,end+1} = myicon('mopen',gHnd.configiconholder,gIcn.mopen,'Draw MO',@() buttondownfunctions('updatebuttondowns','MO'));
    viabilityiconcell{1,end+1} = myicon('rubberscar',gHnd.configiconholder,gIcn.rubberscar,'Manually remove scar segmentation',@() buttondownfunctions('updatebuttondowns','ScarRubber'));
    viabilityiconcell{1,end+1} = myicon('automar',gHnd.configiconholder,gIcn.automar,'Auto MaR',@() mar('auto_Callback'),0);
    viabilityiconcell{1,end+1} = myicon('marpen',gHnd.configiconholder,gIcn.marpen,'Draw MaR',@() buttondownfunctions('updatebuttondowns','MaR'));
    viabilityiconcell{1,end+1} = myicon('rubbermar',gHnd.configiconholder,gIcn.rubbermar,'Manually remove MaR segmentation',@() buttondownfunctions('updatebuttondowns','MaRRubber'));
    viabilityiconcell{1,end+1} = myicon('hidescar',gHnd.configiconholder,gIcn.hidescar,'Hide scar segmentation',@() viewfunctions('hide_Callback'),2);
    viabilityiconcell{1,end+1} = myicon('hidescarextent',gHnd.configiconholder,gIcn.hidescarextent,'Hide scar segmentation extent',@() viewfunctions('hide_Callback'),2);
    viabilityiconcell{1,end+1} = myicon('hidescarmanual',gHnd.configiconholder,gIcn.hidescarmanual,'Hide manual scar interaction',@() viewfunctions('viewhidemanualinteraction_Callback'),2);
    viabilityiconcell{1,end+1} = myicon('hidemar',gHnd.configiconholder,gIcn.hidemar,'Hide MaR segmentation',@() viewfunctions('hide_Callback'),2);
    viabilityiconcell{1,end+1} = myicon('clearscar',gHnd.configiconholder,gIcn.clearscar,'Clear scar segmentation',@() viability('viabilityclear_Callback'),0);  
    viabilityiconcell{1,end+1} = myicon('clearmar',gHnd.configiconholder,gIcn.clearmar,'Clear MaR segmentation',@() mar('clearall_Callback'),0);  
    g.Icons.viabilityiconcell = viabilityiconcell;
    
    %Strain
    strainiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'));
    strainiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'));
    strainiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'));
    strainiconcell{1,end+1} = myicon('contrastbrightness',gHnd.configiconholder,gIcn.contrastbrightness,'Manually change contrast and brightness',@() buttondownfunctions('updatebuttondowns','Contrast'));
    strainiconcell{1,end+1} = myicon('autocontrast',gHnd.configiconholder,gIcn.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    strainiconcell{1,end+1} = myicon('strainsaxstack',gHnd.configiconholder,gIcn.strainsaxstack,'Go to Strain SAX stack',@() segment('viewspecial_Callback','strainsax'),0);
    strainiconcell{1,end+1} = myicon('strainlaxstack',gHnd.configiconholder,gIcn.strainlaxstack,'Go to Strain LAX stack',@() segment('viewspecial_Callback','strainlax'),0);
    strainiconcell{1,end+1} = lviconcell{1,7};

    strainiconcell{1,end+1} = myicon('crop',gHnd.configiconholder,gIcn.crop,'Manual crop',@() buttondownfunctions('updatebuttondowns','crop'));
    strainiconcell{1,end+1} = myicon('drawingguide',gHnd.configiconholder,gIcn.drawingguide,'Drawing guidance for strain analysis',@() straintagging.straintagging('draw_guide_helper',[],0,1),0);
    strainiconcell{1,end+1} = myicon('endopen',gHnd.configiconholder,gIcn.endopen,'Endo pen',@() buttondownfunctions('updatebuttondowns','Endo'));
    strainiconcell{1,end+1} = myicon('epipen',gHnd.configiconholder,gIcn.epipen,'Epi pen',@() buttondownfunctions('updatebuttondowns','Epi'));
    strainiconcell{1,end+1} = myicon('interpendo',gHnd.configiconholder,gIcn.interpendo,'Set interpolation points for Endo',@() buttondownfunctions('updatebuttondowns','EndoInterp'));
    strainiconcell{1,end+1} = myicon('interpepi',gHnd.configiconholder,gIcn.interpepi,'Set interpolation points for Epi',@() buttondownfunctions('updatebuttondowns','EpiInterp'));
    strainiconcell{1,end+1} = myicon('balloonendo',gHnd.configiconholder,gIcn.balloonendo,'Semi automatic Endo tool',@() buttondownfunctions('updatebuttondowns','EndoBalloon'));
%     strainiconcell{1,end+1} = myicon('balloonepi',gHnd.configiconholder,gIcn.balloonepi,'Semi automatic Epi tool',@() lvsegmentation('smartepi_Callback'),0);
    strainiconcell{1,end+1} = myicon('rvendopen',gHnd.configiconholder,gIcn.rvendopen,'RV Endo pen',@() buttondownfunctions('updatebuttondowns','RVEndo'));
    strainiconcell{1,end+1} = myicon('interprvendo',gHnd.configiconholder,gIcn.interprvendo,'Set interpolation points for RV Endo',@() buttondownfunctions('updatebuttondowns','RVEndoInterp'));
    strainiconcell{1,end+1} = myicon('smooth',gHnd.configiconholder,gIcn.smooth,'Smooth current contour [O]',@() tools('smoothsegmentation_Callback'),0);
    strainiconcell{1,end+1} = myicon('point',gHnd.configiconholder,gIcn.point,'Place annotation point',@() buttondownfunctions('updatebuttondowns','Point'));
    strainiconcell{1,end+1} = myicon('strainftsax',gHnd.configiconholder,gIcn.strainftsaxanalysis,'Start FT Strain SAX analysis',@() straintagging.straintagging('init','cine','shortaxis'),0);
    strainiconcell{1,end+1} = myicon('strainftlax',gHnd.configiconholder,gIcn.strainftlaxanalysis,'Start FT Strain LAX analysis',@() straintagging.straintagging('init','cine','longaxis'),0);
    strainiconcell{1,end+1} = myicon('straintagsax',gHnd.configiconholder,gIcn.straintagsaxanalysis,'Start Tagging Strain SAX analysis',@() straintagging.straintagging('init','tagging','shortaxis'),0);
    strainiconcell{1,end+1} = myicon('straintaglax',gHnd.configiconholder,gIcn.straintaglaxanalysis,'Start Tagging Strain LAX analysis',@() straintagging.straintagging('init','tagging','longaxis'),0);
    strainiconcell{1,end+1} = myicon('clearalled',gHnd.permanenticonholder,gIcn.clearalled,'Clear all segmentation except in ED',@() callbackfunctions('segmentclearallbutdiastole_Callback'),0); 
    strainiconcell{1,end+1} = myicon('clearendo',gHnd.configiconholder,gIcn.clearendo,'Clear LV endo in selected slices according to mode',@() segmentation('clearslicesthis_Callback',1,0,0,0),0);
    strainiconcell{1,end+1} = myicon('clearepi',gHnd.configiconholder,gIcn.clearepi,'Clear LV epi in selected slices according to mode',@() segmentation('clearslicesthis_Callback',0,1,0,0),0);
    strainiconcell{1,end+1} = myicon('clearrv',gHnd.configiconholder,gIcn.clearrv,'Clear RV Endo in selected slices according to mode',@() segmentation('clearslicesthis_Callback',0,0,1,0),0);
    strainiconcell{1,end+1} = myicon('clearstrain',gHnd.configiconholder,gIcn.clearstrain,'Clear Strain in selected image stack',@() straintagging.straintagging('clearstrain_Callback'),0);
    strainiconcell{1,end+1} = myicon('clearstrainall',gHnd.configiconholder,gIcn.clearstrainall,'Clear Strain in all image stacks',@() straintagging.straintagging('clearstrainall_Callback'),0);
    g.Icons.strainiconcell = strainiconcell;
    
    %Analysis
    analysisiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'));
    analysisiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'));
    analysisiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'));
    analysisiconcell{1,end+1} = myicon('contrastbrightness',gHnd.configiconholder,gIcn.contrastbrightness,'Manually change contrast and brightness',@() buttondownfunctions('updatebuttondowns','Contrast'));
    analysisiconcell{1,end+1} = myicon('autocontrast',gHnd.configiconholder,gIcn.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    analysisiconcell{1,end+1} = myicon('perfusionstack',gHnd.configiconholder,gIcn.perfusionstack,'Go to perfusion stacks',@() segment('viewspecial_Callback','perfusion'),0);
    analysisiconcell{1,end+1} = lviconcell{1,7};

    analysisiconcell{1,end+1} = myicon('importfromother',gHnd.configiconholder,gIcn.importfromother,'Import LV segmentation from other image stack with snap',@() segmentation('importsegmentationwithsnap_Callback'),0);
    analysisiconcell{1,end+1} = myicon('measure',gHnd.configiconholder,gIcn.measure,'Place measurement',@() buttondownfunctions('updatebuttondowns','Measure'));
    analysisiconcell{1,end+1} = myicon('point',gHnd.configiconholder,gIcn.point,'Place annotation point',@() buttondownfunctions('updatebuttondowns','Point'));
    analysisiconcell{1,end+1} = myicon('roipen',gHnd.configiconholder,gIcn.roipen,'ROI pen',@() buttondownfunctions('updatebuttondowns','Roi'));
    analysisiconcell{1,end+1} = myicon('addroiinlv',gHnd.configiconholder,gIcn.addroiinlv,'Add ROIs to sector of LV wall in selected slices',@() roi('roiaddinsector_Callback'),0);
    analysisiconcell{1,end+1} = myicon('bullseye',gHnd.configiconholder,gIcn.bullseye,'Bullseye analysis [Ctrl-B]',@() reportbullseye,0);
    analysisiconcell{1,end+1} = myicon('AVPD',gHnd.configiconholder,gIcn.AVPD,'AV plane displacement',@() avplane,0);
    
    analysisiconcell{1,end+1} = myicon('perfusion',gHnd.configiconholder,gIcn.perfusion, 'Perfusion analysis',@() perfusion.perfusion,0);
    analysisiconcell{1,end+1} = myicon('perfusionscoring',gHnd.configiconholder,gIcn.perfusionscoring, 'Perfusion scoring',@() perfusion.perfusionscoring,0);
    analysisiconcell{1,end+1} = myicon('reportperslice',gHnd.configiconholder,gIcn.reportperslice, 'Report per slice',@() slicereport,0);
    analysisiconcell{1,end+1} = myicon('model3d',gHnd.configiconholder,gIcn.model3d, 'Show 3D model',@() report3dmodel,0);
    analysisiconcell{1,end+1} = myicon('hidemeasure',gHnd.configiconholder,gIcn.hidemeasure,'Hide measurements',@() viewfunctions('hide_Callback'),2);
    analysisiconcell{1,end+1} = myicon('hidepoint',gHnd.configiconholder,gIcn.hidepoint,'Hide annotation points',@() viewfunctions('hide_Callback'),2);
    analysisiconcell{1,end+1}= roiflowiconcell{1,end-2};%myicon('hideroi',gHnd.configiconholder,gIcn.hideroi,'Hide ROI',@() segment('viewhideroi_Callback'),2);
    analysisiconcell{1,end+1} = myicon('clearmeasure',gHnd.configiconholder,gIcn.clearmeasure,'Clear all measurements',@() callbackfunctions('measureclearall_Callback'),0);  
    analysisiconcell{1,end+1} = myicon('clearpoint',gHnd.configiconholder,gIcn.clearpoint,'Clear all annotation points',@() callbackfunctions('pointclearall_Callback'),0); 
    analysisiconcell{1,end+1} = myicon('clearroi',gHnd.configiconholder,gIcn.clearroi,'Clear selected ROIs',@() roi('roidelete_Callback'),0);  
    analysisiconcell{1,end+1} = myicon('clearallroi',gHnd.configiconholder,gIcn.clearallroi,'Clear all ROIs',@() roi('roiclearall_Callback'),0); 
    g.Icons.analysisiconcell = analysisiconcell;
    
    %Image
    imageiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'));
    imageiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'));
    imageiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'));
    imageiconcell{1,end+1} = myicon('contrastbrightness',gHnd.configiconholder,gIcn.contrastbrightness,'Manually change contrast and brightness',@() buttondownfunctions('updatebuttondowns','Contrast'));
    imageiconcell{1,end+1} = myicon('autocontrast',gHnd.configiconholder,gIcn.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    imageiconcell{1,end+1} = myicon('resetlight',gHnd.configiconholder,gIcn.resetlight,'Reset contrast and brightness',@() segment('resetlight_Callback'),0);
    imageiconcell{1,end+1} = myicon('autocontrastall',gHnd.configiconholder,gIcn.autocontrastall,'Set contrast and brightness to predefined values for all images',@() segment('autocontrastall_Callback'),0);
    imageiconcell{1,end+1} = myicon('resetlightall',gHnd.configiconholder,gIcn.resetlightall,'Reset contrast and brightness for all imagestacks',@() segment('resetlightall_Callback'),0);
    imageiconcell{1,end+1} = lviconcell{1,7};

    imageiconcell{1,end+1} = myicon('crop',gHnd.configiconholder,gIcn.crop,'Manual crop',@() buttondownfunctions('updatebuttondowns','crop'));
    imageiconcell{1,end+1} = myicon('cropall',gHnd.configiconholder,gIcn.cropall,'Auto crop all',@() autocropallgui('init'),0);
    imageiconcell{1,end+1} = myicon('croplv',gHnd.configiconholder,gIcn.croplv,'LV crop',@() lvsegmentation('croplvall',0),0);
    imageiconcell{1,end+1} = myicon('cineplay',gHnd.configiconholder,gIcn.cineplay,'Open cine tool',@() segment('cinetool_Callback'),2);
    imageiconcell{1,end+1} = myicon('movie',gHnd.configiconholder,gIcn.movie,'Open movie tool',@() export('exportmovierecorder_Callback'),0);
    imageiconcell{1,end+1} = myicon('click3d',gHnd.configiconholder,gIcn.click3d,'Set 3D point',@() buttondownfunctions('updatebuttondowns','click3d'));
    imageiconcell{1,end+1} = myicon('rotate90',gHnd.configiconholder,gIcn.rotate90,'Rotate 90 degrees clockwise',@() tools('rotate90right_Callback'),0);
    imageiconcell{1,end+1} = myicon('mpr',gHnd.configiconholder,gIcn.mpr,'Reconstruct image stack',@() reformater,0);
    imageiconcell{1,end+1} = myicon('mergestacks',gHnd.configiconholder,gIcn.mergestacks,'Merge stacks',@() mergestacks,0);    
    imageiconcell{1,end+1} = myicon('imageinfo',gHnd.configiconholder,gIcn.imageinfo,'View and adjust image info',@() tools('imageinfo_Callback'),0);
    imageiconcell{1,end+1} = myicon('patientinfo',gHnd.configiconholder,gIcn.patientinfo,'View and adjust patient info',@() tools('viewpatientinfo_Callback'),0);    
    imageiconcell{1,end+1}=roiflowiconcell{1,end-2};%myicon('hidetext',gHnd.configiconholder,gIcn.hidetext,'Hide text',@() viewfunctions('hide_Callback'),2);
    imageiconcell{1,end+1} = myicon('hideplus',gHnd.configiconholder,gIcn.hideplus,'Hide center cross',@() viewfunctions('hide_Callback'),2);
    imageiconcell{1,end+1} = myicon('hideintersections',gHnd.configiconholder,gIcn.hideintersections,'Hide intersection lines',@() viewfunctions('hide_Callback'),2);
    imageiconcell{1,end+1} = myicon('hideothercontour',gHnd.configiconholder,gIcn.hideothercontour,'Hide other contour points',@() viewfunctions('hide_Callback'),2);
    g.Icons.imageiconcell = imageiconcell;
    
    %TXmap
    txmapiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'));
    txmapiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'));
    txmapiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'));
    txmapiconcell{1,end+1} = myicon('contrastbrightness',gHnd.configiconholder,gIcn.contrastbrightness,'Manually change contrast and brightness',@() buttondownfunctions('updatebuttondowns','Contrast'));
    txmapiconcell{1,end+1} = myicon('autocontrast',gHnd.configiconholder,gIcn.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0);
    txmapiconcell{1,end+1} = myicon('resetlight',gHnd.configiconholder,gIcn.resetlight,'Reset contrast and brightness',@() segment('resetlight_Callback'),0);    
    txmapiconcell{1,end+1} = lviconcell{1,7};
 
    txmapiconcell{1,end+1} = myicon('importfromother',gHnd.configiconholder,gIcn.importfromother,'Import LV segmentation from cine to Tx map image stack',@() segmentation('importfromcine2txmap_Callback'),0);
    txmapiconcell{1,end+1} = myicon('roipen',gHnd.configiconholder,gIcn.roipen,'ROI pen',@() buttondownfunctions('updatebuttondowns','Roi'));
    txmapiconcell{1,end+1} = myicon('addroiinlv',gHnd.configiconholder,gIcn.addroiinlv,'Add ROIs to sector of LV wall in selected slices',@() roi('roiaddinsector_Callback'),0);
    txmapiconcell{1,end+1} = myicon('T1',gHnd.configiconholder,gIcn.T1, 'T1 analysis',@() txmap('init',1),0);
    txmapiconcell{1,end+1} = myicon('T2',gHnd.configiconholder,gIcn.T2, 'T2 analysis',@() txmap('init',2),0);
    txmapiconcell{1,end+1} = myicon('T2star',gHnd.configiconholder,gIcn.T2star, 'T2* analysis',@() t2star.t2star,0);
    txmapiconcell{1,end+1} = myicon('ecv',gHnd.configiconholder,gIcn.ecv, 'ECV analysis',@() ecv('init_Callback'),0);
    txmapiconcell{1,end+1} = myicon('T1precolormap',gHnd.configiconholder,gIcn.T1precolormap, 'T1 pre colormap',@() tools('setcolormap_Callback','t1pre'),0);
    txmapiconcell{1,end+1} = myicon('T1postcolormap',gHnd.configiconholder,gIcn.T1postcolormap, 'T1 post colormap',@() tools('setcolormap_Callback','t1post'),0);
    txmapiconcell{1,end+1} = myicon('T2colormap',gHnd.configiconholder,gIcn.T2colormap, 'T2 colormap',@() tools('setcolormap_Callback','t2'),0);
    txmapiconcell{1,end+1} = myicon('T2starcolormap',gHnd.configiconholder,gIcn.T2starcolormap, 'T2* colormap',@() tools('setcolormap_Callback','t2star'),0);
    txmapiconcell{1,end+1} = myicon('ecvcolormap',gHnd.configiconholder,gIcn.ecvcolormap, 'ECV colormap',@() tools('setcolormap_Callback','ecv'),0);
    txmapiconcell{1,end+1} = myicon('greycolormap',gHnd.configiconholder,gIcn.greycolormap, 'Grey colormap',@() tools('setcolormap_Callback','gray'),0);
    txmapiconcell{1,end+1}= roiflowiconcell{1,end-2};
    txmapiconcell{1,end+1} = myicon('clearroi',gHnd.configiconholder,gIcn.clearroi,'Clear selected ROIs',@() roi('roidelete_Callback'),0);
    txmapiconcell{1,end+1} = myicon('clearallroi',gHnd.configiconholder,gIcn.clearallroi,'Clear all ROIs',@() roi('roiclearall_Callback'),0);

    g.Icons.txmapiconcell = txmapiconcell;
    end
    
     %----------------------------------------
    function initribbonplaceholder(varargin)
    %--------------------------------------

    g = varargin{1};
    gHndtoggle = g.Handles.toggleiconholder;
    gIcntoggle = g.Icons.toggleicons;
    %initiate iconholder for toggle axes 
    iconCell = cell(1,7);
    iconCell{1} = myicon('ribbonlv',gHndtoggle,gIcntoggle.ribbonlvoff,...
      'LV', @() g.togglebuttonLV_Callback,1,1,gIcntoggle.ribbonlvon);
    iconCell{2} = myicon('ribbonrv',gHndtoggle,gIcntoggle.ribbonrvoff,...
      'RV',@() g.togglebuttonRV_Callback,1,1,gIcntoggle.ribbonrvon);
    iconCell{3} = myicon('ribbonflow',gHndtoggle,gIcntoggle.ribbonflowoff,...
      'ROI/FLOW',@() g.togglebuttonROIFLOW_Callback,1,1,gIcntoggle.ribbonflowon);
    iconCell{4} = myicon('ribbonviability',gHndtoggle,gIcntoggle.ribbonviabilityoff,...
      'Viability',@() g.togglebuttonVia_Callback,1,1,gIcntoggle.ribbonviabilityon);
    iconCell{5} = myicon('ribbontxmap',gHndtoggle,gIcntoggle.ribbontxmapoff,...
      'TX maps and ECV',@() g.togglebuttontxmap_Callback,1,1,gIcntoggle.ribbontxmapon);
    iconCell{6} = myicon('ribbonstrain',gHndtoggle,gIcntoggle.ribbonstrainoff,...
      'Strain',@() g.togglebuttonStrain_Callback,1,1,gIcntoggle.ribbonstrainon);
    iconCell{7} = myicon('ribbonanalysis',gHndtoggle,gIcntoggle.ribbonanalysisoff,...
      'Analysis',@() g.togglebuttonAnalysis_Callback,1,1,gIcntoggle.ribbonanalysison);
    iconCell{8} = myicon('ribbonimage',gHndtoggle,gIcntoggle.ribbonimageoff,...
      'Image',@() g.togglebuttonImage_Callback,1,1,gIcntoggle.ribbonimageon);
 
    gHndtoggle.add(iconCell);   
    pos = plotboxpos(gHndtoggle.axeshandle);
    currentpos = get(gHndtoggle.axeshandle,'position');
    set(gHndtoggle.axeshandle,'position',currentpos-[pos(1),0,0,0]);
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
        
    %--- Initialize all menu items
    segment('initmenu');
    
    %--- Add extra menu items (if applicable)
    addextramenuitems(g);

    end
        
    %------------------------------------------
    function updatelvreport(g)
    %------------------------------------------
    global SET NO
    %for segment we always use the values in the NO
    no = NO;    
    
    c = cell(1,8);
    if ~isempty(SET) && ~isempty(NO)
      lvmEDT = SET(no).LVM(SET(no).EDT)*1.05;
      lvmEST = SET(no).LVM(SET(no).EST)*1.05;
      edv = SET(no).EDV;
      esv = SET(no).ESV;
      sv = SET(no).SV;
      co = (SET(no).HeartRate*SET(no).SV)/100;
      if lvmEDT > 0 && lvmEDT < 1
        %convert to milligram 
        lvmEDT = 1000*lvmEDT;
        lvmEST = 1000*lvmEST;
        g.Handles.lvheadertext.String{1} = 'ED/ES LVM (mg)';
      else
        g.Handles.lvheadertext.String{1} = 'ED/ES LVM (g)';
      end
   
      if edv > 0 && edv < 1
        %convert to mikroliters
        edv = 1000*edv;
        g.Handles.lvheadertext.String{2} = ['EDV (',char(181),'l)'];
        esv = 1000*esv;
        g.Handles.lvheadertext.String{3} = ['ESV (',char(181),'l)'];
        sv = 1000*sv;
        g.Handles.lvheadertext.String{4} = ['SV (',char(181),'l)'];
        co = 1000*co;
        g.Handles.lvheadertext.String{6} = 'CO (ml/min)';
      else
        g.Handles.lvheadertext.String{2} = 'EDV (ml)';
        g.Handles.lvheadertext.String{3} = 'ESV (ml)';
        g.Handles.lvheadertext.String{4} = 'SV (ml)';
        g.Handles.lvheadertext.String{6} = 'CO (l/min)';
      end
      
      c{1} = round(lvmEDT);
      c{2} = round(lvmEST);
      c{3} = round(edv);
      c{4} = round(esv);
      c{5} = round(sv);
      c{6} = round(100*SET(no).EF);
      c{7} = round(co)/10;
      c{8} = round(SET(no).HeartRate);
    end
    
    c_inds = cellfun(@(x) isempty(x) || isnan(x) || x == 0,c);
    [c{c_inds}] = deal('---');
    c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
    g.Handles.lvreporttext.String = sprintf('%s / %s \n%s \n%s \n%s \n%s \n%s \n%s',c{:});
    
    %adjust fontsize according to the available space on screen
    pos = getpixelposition(g.Handles.lvreporttext);
    widthinpixel = pos(3);
    cmperwidth  = 0.053*widthinpixel/10;
    cmperlength = 0.0265*pos(4)/9.5;
    g.fontsizeincm = min(cmperwidth,cmperlength);
    set(g.Handles.lvreporttext,'FontSize',g.fontsizeincm);
    g.Handles.lvheadertext.FontSize = g.fontsizeincm;

    end
    
    %------------------------------------------
    function updatervreport(g)
    %------------------------------------------
    global SET NO
    %for segment cmr we always use the values in the NO
    no = NO;
    
    c = cell(1,6);
    
    if ~isempty(SET) && ~isempty(no)
      rvmEDT = SET(no).RVM(SET(no).EDT)*1.05;
      rvmEST = SET(no).RVM(SET(no).EST)*1.05;
      rvedv = SET(no).RVEDV;
      rvesv = SET(no).RVESV;
      rvsv = SET(no).RVSV;
      if rvmEDT > 0 && rvmEDT < 1
        %convert to milligram 
        rvmEDT = 1000*rvmEDT;
        rvmEST = 1000*rvmEST;
        g.Handles.rvheadertext.String{1} = 'ED/ES LVM (mg)';
      else
        g.Handles.rvheadertext.String{1} = 'ED/ES LVM (g)';
      end
   
      if rvedv > 0 && rvedv < 1
        %convert to mikroliters
        rvedv = 1000*rvedv;
        g.Handles.rvheadertext.String{2} = ['RVEDV (',char(181),'l)'];
        rvesv = 1000*rvesv;
        g.Handles.rvheadertext.String{3} = ['RVESV (',char(181),'l)'];
        rvsv = 1000*rvsv;
        g.Handles.rvheadertext.String{4} = ['RVSV (',char(181),'l)'];
      else
        g.Handles.rvheadertext.String{2} = 'RVEDV (ml)';
        g.Handles.rvheadertext.String{3} = 'RVESV (ml)';
        g.Handles.rvheadertext.String{4} = 'RVSV (ml)';
      end
      
      c{1} = round(rvmEDT);
      c{2} = round(rvmEST);
      c{3} = round(rvedv);
      c{4} = round(rvesv);
      c{5} = round(rvsv);
      c{6} = round(100*SET(no).RVEF);
    end
    
    c_inds = cellfun(@(x) isempty(x) || isnan(x) || x == 0,c);
    [c{c_inds}] = deal('---');
    c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
    g.Handles.rvreporttext.String = sprintf('%s / %s \n%s \n%s \n%s \n%s',c{:});
    set(g.Handles.rvreporttext,'FontSize',g.fontsizeincm);
    g.Handles.rvheadertext.FontSize = g.fontsizeincm;

    end
    
    %------------------------------------------
    function updateflowreport(g)
    %------------------------------------------
    global SET NO
    %for segment cmr we always use the values in the FlowNO
    no = NO;
    c = cell(1,6);
%     set(g.Handles.flowheaderroitext,'String','ROI-x');
     if ~isempty(SET) && ~isempty(no)
      if isfield(SET(no), 'Parent')
        if ~isempty(SET(no).Parent)
          %reset no to the magnitude image
          no = SET(no).Parent;
        end
      end
     end
    roinbr =[];
    if ~isempty(SET) && ~isempty(no)
      if ~isempty (SET(no).RoiCurrent)
        roinbr = SET(no).RoiCurrent;
      end
    end
    
    
    
  
%     if ~isempty(SET) && ~isempty(no) && ~isempty(SET(no).Flow)...
%         &&  ~isempty(SET(no).RoiCurrent) && ~isempty(SET(no).Flow.Result)...
    if ~isempty(SET)&& ~isempty(no) && ~isempty(SET(no).Flow) && ~isempty(roinbr) &&...
        ~isempty(SET(no).Flow.Result) && ~(length(SET(no).Flow.Result)< roinbr)...
        && ~strcmp(SET(no).Roi(roinbr).Name, 'Static tissue')
      
        set(g.Handles.flowheaderroitext,'String',SET(no).Roi(SET(no).RoiCurrent).Name);
        if SET(no).Flow.Result(roinbr).nettotvol > 100
          c{1} = round(SET(no).Flow.Result(roinbr).nettotvol);
        else
          c{1} = round(SET(no).Flow.Result(roinbr).nettotvol*100)/100;
        end
        c{2} = round(SET(no).Flow.Result(roinbr).netforwardvol*100)/100;
        c{3} = round(SET(no).Flow.Result(roinbr).netbackwardvol*100)/100;
        if SET(no).Flow.Result(roinbr).regfrac > 100
          c{4} = round(SET(no).Flow.Result(roinbr).regfrac);
        else
          c{4} = round(SET(no).Flow.Result(roinbr).regfrac*100)/100;
        end

        %define heart rate
        if not(isfield(SET(no).Flow,'HeartRate')) || isempty(SET(no).Flow.HeartRate)
          temphr = (60/(SET(no).TSize*SET(no).TIncr));
          if temphr < 35
            reportflow('defineheartrate',no);
          else
            SET(no).Flow.HeartRate = temphr;
          end
        end
        c{5} = round(SET(no).Flow.Result(roinbr).nettotvol*SET(no).Flow.HeartRate/1000*100)/100;
        c{6} = round(SET(no).Flow.HeartRate); %SET(no).HeartRate;
    else
%       %define heart rate
%         if not(isempty(no)) && (not(isfield(SET(no).Flow,'HeartRate')) || isempty(SET(no).Flow.HeartRate))
%           temphr = (60/(SET(no).TSize*SET(no).TIncr));
%           if temphr < 35
%             reportflow('defineheartrate',no);
%           else
%             SET(no).Flow.HeartRate = temphr;
%           end
%         end
     
    end
    
    if ~isempty(roinbr)
%       c{6} = round(SET(no).Flow.HeartRate);      
      set(g.Handles.flowheaderroitext,'String',SET(no).Roi(roinbr).Name);
    end
    
    c_inds = cellfun(@(x) isempty(x) || isnan(x) || x == 0,c);
    [c{c_inds}] = deal('---');
    
    c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
    g.Handles.flowreporttext.String = sprintf('%s \n%s \n%s \n%s \n%s \n%s \n%s',c{:});
    set(g.Handles.flowreporttext,'FontSize',g.fontsizeincm);
    set(g.Handles.flowheadertext,'FontSize',g.fontsizeincm);
    set(g.Handles.flowheaderroitext,'FontSize',g.fontsizeincm);
    end

    %-------------------------
    function updateflowaxes(g)
    %-------------------------
    %Update the flow report axes in main GUI
    global SET NO
    roinbr=[];
    no = NO;
    if ~isempty(SET) && ~isempty(no)
      if isfield(SET(no), 'Parent')
        if ~isempty(SET(no).Parent)
          % reset no to the magnitude image
          no = SET(no).Parent;
        end
      end
    end
    if ~isempty(SET) && ~isempty(no)
      roinbr = SET(no).RoiCurrent;
    end
    
    if isempty(SET) || isempty(no) || SET(no).TSize<2 || isempty(roinbr)|| isempty(SET(no).Flow) || isempty(SET(no).Flow.Result)  || length(SET(no).Flow.Result)< roinbr || isempty(SET(no).Flow.Result(roinbr)) || isempty(SET(no).Flow.Result(roinbr).netflow) 
      set([g.Handles.flowaxes;g.Handles.flowaxes.Children;...
        g.Handles.flowaxes.XLabel;g.Handles.flowaxes.YLabel],'Visible','off');
     
    else        
      %Time resolved
      t = SET(no).TimeVector*1000;
    
      set([g.Handles.flowaxes;g.Handles.flowaxes.Children;...
        g.Handles.flowaxes.XLabel;g.Handles.flowaxes.YLabel],'Visible','on');
            
      %set the xlim of the plot
      xlim(g.Handles.flowaxes,[t(1),t(end)])
      
      %flow curve
      netflow = SET(no).Flow.Result(roinbr).netflow; 
      
      %plot in color of roi
      curvecolor=SET(no).Roi(roinbr).LineSpec(1);
      
      g.Handles.flowcurve.XData = t;
      g.Handles.flowcurve.YData = netflow;
      g.Handles.flowcurve.Color = curvecolor;

      % set y-lim dependent on netflow values
      minflow = min(netflow);
      maxflow = max(netflow);
      set(g.Handles.flowaxes,'ylim',[floor((minflow-50)/100)*100 ceil((maxflow+50)/100)*100]);
      
      %outer time bars
      g.Handles.outerbarsflow.XData = [t(1) t(1) nan t(end) t(end)]; 
      g.Handles.outerbarsflow.YData = [ylim(g.Handles.flowaxes) nan ylim(g.Handles.flowaxes)]; 
      
      %current time
      g.Handles.timebarflow.XData = [t(SET(no).CurrentTimeFrame) t(SET(no).CurrentTimeFrame)];
      g.Handles.timebarflow.YData = ylim(g.Handles.flowaxes);
        
      %the flow axeshelpers
      temp = get(g.Handles.flowaxes,'ylim');
      ttemp = [t;t;nan(size(t))];
      temp = [repmat(temp(1),size(t));repmat(temp(1)+0.1*(temp(2)-temp(1)),size(t));nan(size(t))];
      g.Handles.flowaxeshelpers.XData = ttemp(:);
      g.Handles.flowaxeshelpers.YData = temp(:);
      
      %the default assumption
      g.Handles.flowtext.Position = [nan nan];
      if ~isempty(SET(SET(no).Flow.PhaseNo).Flow.PhaseCorr)
          if ~isfield(SET(SET(no).Flow.PhaseNo).Flow,'PhaseCorrTimeResolved')
              mywarning('Incompatible eddy current correction. Correction reset.',g.GUI.Segment);
          else
              g.Handles.flowtext.Position = [max(get(g.Handles.flowaxes,'XLim')),double(max(get(g.Handles.flowaxes,'YLim')))];
          end
      end
     end
    end
    
    %------------------------------------------
    function updatemeasurementreport(g)
    %------------------------------------------
      global SET NO
      %for segment we always use the values in the NO
      no = NO;
      
      if ~isempty(SET) && ~isempty(no)
      %different measurements we want to display
      hasscar = ~isempty(SET(no).Scar) && ~isempty(SET(no).Scar.Percentage);
      hasmar = ~isempty(SET(no).MaR) && ~isempty(SET(no).MaR.Result);
      hasatrialscar = isfield(SET(no),'AtrialScar') && ~isempty(SET(no).AtrialScar) && ~isempty(SET(no).AtrialScar.Percentage);
      else
        hasscar = false;
        hasmar = false;
        hasatrialscar = false;
      end
      
      if ~any([hasscar hasmar hasatrialscar])
        g.Handles.measurementuipanel.Visible = 'off';
        g.Handles.flowuipanel.Visible = 'on';
      else
        g.Handles.measurementuipanel.Visible = 'on';
        g.Handles.flowuipanel.Visible = 'off';
      end
      
      if hasscar
        scarpro2mlcoeff=sum(SET(no).Scar.MyocardMask(:))/100/1000*SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceGap+SET(no).SliceThickness);
        scarml=round(10*SET(no).Scar.Percentage*scarpro2mlcoeff)/10;
      end
      
      if hasmar && hasscar
        c = cell(1,4);
        c{1} = scarml;
        c{2} = round(10*SET(no).Scar.Percentage)/10;
        c{3} = round(10*SET(no).Scar.MOPercentage)/10;
        c{4} = (round(SET(no).MaR.Percentage(SET(no).EDT)));
        c{5} = (round(SET(no).MaR.Percentage(SET(no).EST)));
        
        c_inds = cellfun(@(x) isempty(x) || isnan(x),c);
        [c{c_inds}] = deal('---');
        c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
        g.Handles.measurementheadertext.String = sprintf('%s\n%s\n%s\n%s','Scar (ml)','Scar (%)','MO (%)','ED/ES MaR (%)');
        g.Handles.measurementreporttext.String = sprintf('%s\n%s\n%s\n%s/%s',c{:});
        
      elseif hasscar
        c = cell(1,3);
        c{1} = scarml;
        c{2} = round(10*SET(no).Scar.Percentage)/10;
        c{3} = round(10*SET(no).Scar.MOPercentage)/10;
        
        c_inds = cellfun(@(x) isempty(x) || isnan(x),c);
        [c{c_inds}] = deal('---');
        c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
        g.Handles.measurementheadertext.String = sprintf('%s\n%s\n%s','Scar (ml)','Scar (%)','MO (%)');
        g.Handles.measurementreporttext.String = sprintf('%s\n%s\n%s',c{:});
        
      elseif hasmar    
        c = cell(1,2);
        c{1} = round(SET(no).MaR.Percentage(SET(no).EDT));
        c{2} = round(SET(no).MaR.Percentage(SET(no).EST));
          
        c_inds = cellfun(@(x) isempty(x) || isnan(x),c);
        [c{c_inds}] = deal('---');
        c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
        g.Handles.measurementheadertext.String = sprintf('%s','ED/ES MaR (%)');
        g.Handles.measurementreporttext.String = sprintf('%s/%s',c{:});
        
      elseif hasatrialscar          
        c = cell(1,1);
        c{1} = round(10*SET(no).AtrialScar.Percentage)/10;
        
        c_inds = cellfun(@(x) isempty(x) || isnan(x),c);
        [c{c_inds}] = deal('---');
        c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
        g.Handles.measurementheadertext.String = sprintf('%s','Atrial Scar (%)');
        g.Handles.measurementreporttext.String = sprintf('%s',c{:});
      end
      viability('sliderupdate')
      set(g.Handles.measurementheadertext,'FontSize',g.fontsizeincm);
      set(g.Handles.measurementreporttext,'FontSize',g.fontsizeincm);
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
    
%     %-------------------------------------------------
%     function [stri,lstr] = measureasklabel(~,measureind)
%     %-------------------------------------------------
%     %Ask measurement
%     
%     s = myinputdlg({'Enter name'},'Name',1,{sprintf('Measure_%d',measureind)});
%     if isempty(s)
%       stri = '';
%       lstr = '';
%     else
%       stri = s{1};
%       lstr = s{1};
%     end
%     
%     end
    %-------------------------------------------------
    function [stri,lstr] = measureasklabel(g,measureind)
    %-------------------------------------------------
    %Ask measurement
    
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
        s = myinputdlg({'Enter name'},'Name',1,{sprintf('Measure_%d',measureind)});
        if isempty(s)
          stri = '';
          lstr = '';
        else
          stri = s{1};
          lstr = s{1};
        end
    end
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
    viewfunctions('switchtimeframe',1,true); %segment('setcurrenttimeframe',tf2);
    segment('updatevolume');
    end
  end
  
end