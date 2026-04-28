classdef segmentgui < maingui %skiptranslation
  %Contains methods and properties that are specific to (original) Segment.

  %#ok<*GVMIS>

  properties
    AxesTables = [];
    %     LVNO = [];
    %     RVNO = [];
    %     FlowNO = [];
    %     FlowROI = [];
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

     %------------------------
    function updatestackmode(g,isvisible)
      %----------------------
      % Init stack mode buttons (all/relevant stack)
      % Copied from segmentmr
      
      if isvisible
        visiblestate = 'on';
      else
        visiblestate = 'off';
      end
      set(g.Handles.stacksuibuttongroup,'Visible',visiblestate);
      if isvisible        
        % Update string for allstackstogglebutton
        allstackstr = upper(dprintf('All'));
        set(g.Handles.allstackstogglebutton,'String',allstackstr);

        % Update string for relevantstackstogglebutton
        relevantstackstr = upper(dprintf('Mode'));
        set(g.Handles.relevantstackstogglebutton,'String',relevantstackstr);

        selectedbtn = g.Handles.stacksuibuttongroup.SelectedObject;
        if (contains(selectedbtn.Tag,'all') && g.ShowRelevantStacksOnly) ...
            || (contains(selectedbtn.Tag,'relevant') && ~g.ShowRelevantStacksOnly)
          viewfunctions('togglestackmode')
        end
      end
    end

    %--------------------------
    function requestfocus(g)
      %--------------------------
      %Take away focus from any uicontrols
      myrequestfocus(g.fig); %Request focus to avoid problem with sliders
    end

    %---------------------------------
    function g = initIFUmenu(g)
      %---------------------------------
      %Initialise menu item for Instructions for Use
      parentmenu = g.Handles.usermanuals;
      
      g.Handles.researchmanual = uimenu(parentmenu,...
        'Label','User Manual Research', ...
        'Callback','segmenthelp(''usermanual_Callback'')');

      if g.isFDAVersion
        g.Handles.fdamanual = uimenu(parentmenu,...
          'Label','User Manual FDA Version', ...
          'Callback','segmenthelp(''usermanual_Callback'',true)');
      end
    end

    %-------------------------------
    function g = initcontextmenus(g)
      %-----------------------------
      %Initialise context menus
      g = g.initlaracontextmenu;
      g = g.initrvcontextmenu;
      g = g.initlvcontextmenu;
    end

    %----------------------
    function initbatchoperationsmenu(g)
      %----------------------
      % add menus only to segment
      if isa(g.GUI.Segment.handles.utilitymenu,'matlab.ui.container.Menu')
        parentmainmenu = g.GUI.Segment.handles.utilitymenu;
        parenth = findobj(parentmainmenu.Children, 'Text', 'Batch Operations on .mat Files');
        if isempty(parenth)
          parenth = parentmainmenu;
        end
        g.Handles.batchlaxsegstrainmenu = uimenu('Parent',parenth,...
          'Callback','utility(''batchlaxsegstrain_Callback'')',...
          'Label','Batch LAX segmentation and Strain MITT',...
          'Tag','batchlaxsegstrainmenu');
        if ~isdeployed || strcmpi(g.LicenseNumber,'p5jK3TcR')
          g.Handles.batchstrainmittmenu = uimenu('Parent',parenth,...
            'Callback','utility(''batchstrainmitt_Callback'')',...
            'Label','Batch Strain MITT',...
            'Tag','batchstrainmittmenu');
        end
        if ~isdeployed
          g.Handles.batchupdatestrainmittmenu = uimenu('Parent',parenth,...
            'Callback','utility(''batchupdatestrainmitt_Callback'')',...
            'Label','Batch Update Strain MITT',...
            'Tag','batchupdatestrainmittmenu');
        end
      end
    end

    %----------------------
    function initexportmenus(g)
      %----------------------
      % add menus only to segment
      if isa(g.GUI.Segment.handles.utilitymenu,'matlab.ui.container.Menu')
        parentmainmenu = g.GUI.Segment.handles.utilitymenu;
        parenth = findobj(parentmainmenu.Children, 'Text', 'Batch Export from .mat Files');
        if isempty(parenth)
          parenth = parentmainmenu;
        end

        % export Strain MITT
        if g.checkstrainlicense
          exportstrainmittmenu = uimenu('Parent',parenth,'Separator','on','Label','Export Strain MITT',...
            'Tag','exportstrainmittmenu');
          g.Handles.batchexportstrainmittsaxlvmenu = uimenu('Parent',exportstrainmittmenu,...
            'Callback','export(''exportmultiplestrainmitt_Callback'',''sax'',true)',...
            'Label','Export SAX',...
            'Tag','batchexportstrainmittsaxlvmenu');
          g.Handles.batchexportstrainmittlaxlvmenu = uimenu('Parent',exportstrainmittmenu,...
            'Callback','export(''exportmultiplestrainmitt_Callback'',''lax'',true)',...
            'Label','Export LAX',...
            'Tag','batchexportstrainmittlaxlvmenu');
          g.Handles.batchexportstrainmittlaxgraph = uimenu('Parent',exportstrainmittmenu,...
            'Label','Export global strain LAX graphs',...
            'Tag','batchexportstrainmittlaxgraph',...
            'Callback','export(''exportmultiplestrainmittgraphs_Callback'', ''lax'')');
          g.Handles.batchexportstrainmittsaxgraph = uimenu('Parent',exportstrainmittmenu,...
            'Label','Export global strain SAX graphs',...
            'Tag','batchexportstrainmittsaxgraph',...
            'Callback','export(''exportmultiplestrainmittgraphs_Callback'', ''sax'')');
          g.Handles.exportstrainmittmenu = exportstrainmittmenu;
          separatoroldstrain = 'off';
        else
          separatoroldstrain = 'on';
        end

        % export Strain 1st Generation
        labeloldstrain = 'Export Strain (1st Generation)';
        exportstrain3p0menu = uimenu('Parent',parenth,'Separator',separatoroldstrain,'Label',labeloldstrain,...
          'Tag','exportstrain3p0menu');
        g.Handles.batchexportstrainftmenu = uimenu('Parent',exportstrain3p0menu,...
          'Callback','export(''exportmultiplestrain_Callback'',''cine'', ''true'')',...
          'Label','Export Feature Tracking Strain',...
          'Tag','batchexportstrainftmenu');
        g.Handles.batchexportstrainftrvmenu = uimenu('Parent',exportstrain3p0menu,...
          'Callback','export(''exportmultiplestrainRV_Callback'',''cine'')',...
          'Label','Export Feature Tracking Strain RV',...
          'Tag','batchexportstrainftrvmenu' );
        g.Handles.batchexportstraingraphsftmenu = uimenu('Parent',exportstrain3p0menu,...
          'Callback','export(''exportmultiplestraingraphs_Callback'',''cine'')',...
          'Label','Export Feature Tracking Strain Graphs',...
          'Tag','batchexportstraingraphsftmenu');
        g.Handles.batchexportrvstraingraphsftmenu = uimenu('Parent',exportstrain3p0menu,...
          'Callback','export(''exportmultiplestrainRVgraphs_Callback'',''cine'')',...
          'Label','Export Feature Tracking Strain RV Graphs',...
          'Tag','batchexportrvstraingraphsftmenu');
        uimenu('Parent',exportstrain3p0menu,...
          'Callback','export(''exportmultiplestrain_Callback'',''tagging'')',...
          'Label','Export Tagging Strain',...
          'Tag','batchexportstraintaggingmenu' );
        g.Handles.exportstrain3p0menu = exportstrain3p0menu;

        % export AVPD
        exportavpdmenu = uimenu('Parent',parenth,'Separator','on','Label','Export longitudinal pumping variables',...
          'Tag','exportavpdmenu');
        g.Handles.exportleftavpdmenu = uimenu('Parent',exportavpdmenu,...
          'Callback','avplane(''exporttrackedlvavpd'', ''multiple'')',...
          'Label','Export Left AV-Plane Displacement');
        g.Handles.exportrightavpdmenu = uimenu('Parent',exportavpdmenu,...
          'Callback','avplane(''exporttrackedrvavpd'', ''multiple'')',...
          'Label','Export Right AV-Plane Displacement');
        g.Handles.exportallavpdmenu = uimenu('Parent',exportavpdmenu,...
          'Callback','avplane(''exportalldata'')',...
          'Label','Export all AVPD and longitudinal pumping');
        g.Handles.exportavpdmenu = exportavpdmenu;

      end
    end

      %---------------------------------
      function g = initAILVmenu(g)
      %---------------------------------
      %Initialise menu item for LV

        parentmenu = g.Handles.lvmenu;
        g.Handles.lvaisaxmenu = uimenu(parentmenu,...
          'Label','Legacy AI-based LV segmentation', ...
          'Callback','lvsegmentationml(''doaltlv2_Callback'')'); 

        g.Handles.lvaisaxedesmenu = uimenu(parentmenu,...
          'Label','Legacy AI-based LV segmentation for ED/ES', ...
          'Callback','lvsegmentationml(''doaltlv2_Callback'',true)');

        g.Handles.lvaisaxslicesmenu = uimenu(parentmenu,...
          'Label','Legacy AI-based LV segmentation for selected slices', ...
          'Callback','lvsegmentationml');
      end

       %---------------------------------
    function g = initautoloadermenu(g)
      %---------------------------------
      %Initialise menu item for autoloader in Segment
      parentmenu = g.Handles.autoloadermenu;
      % check license
      [~, ~, ~, enablevalue] = checkmodulelicense('J');
      g.Handles.autoloadfromfoldermenu = uimenu(parentmenu,...
        'Label',dprintf('Create Segment files from DICOM (default)'), ...
        'Tag','autoloadfromfoldermenu',...
        'Enable',enablevalue,...
        'Callback',['autoloader.autoloadfromfolder(false)' ...
        '']);
      g.Handles.autoloadfromsinglefoldermenu = uimenu(parentmenu,...
        'Label',dprintf('Create Segment files from DICOM (patients structured in subfolders)'), ...
        'Tag','autoloadfromsinglefoldermenu',...
        'Enable',enablevalue,...
        'Callback',['autoloader.autoloadfromfolder(true)' ...
        '']);
    end

    %----------------------
    function initmlnetworks(g)
    %----------------------

      %LV in LGE
      n = 1;
      g.Networks(n).Name = 'LV in LGE';
      g.Networks(n).NV = {'1338.nv'};
      g.Networks(n).FileParts = 10;
      g.Networks(n).LicenseCode = ''; %non existing code
      g.Networks(n).PreInstalled = false;
      g.Networks(n).Visible = true;
      g.Networks(n).Modality = 'MR';
      g.Networks(n).Description = 'Trained to segment left ventricle in LGE SAX images.';
      g.Networks(n).ImageFile = '';

      %LAX 2CH
      n = n + 1;
      g.Networks(n).Name = 'LAX 2CH';
      g.Networks(n).NV = {'5721.nv'};
      g.Networks(n).FileParts = 10;
      g.Networks(n).LicenseCode = ''; %non existing code
      g.Networks(n).PreInstalled = false;
      g.Networks(n).Visible = true;
      g.Networks(n).Modality = 'MR';
      g.Networks(n).Description = 'Trained to segment LV, and LA, in 2CH LAX images.';
      g.Networks(n).ImageFile = '';

      %LAX 3CH
      n = n + 1;
      g.Networks(n).Name = 'LAX 3CH';
      g.Networks(n).NV = {'5723.nv'};
      g.Networks(n).FileParts = 10;
      g.Networks(n).LicenseCode = ''; %non existing code
      g.Networks(n).PreInstalled = false;
      g.Networks(n).Visible = true;
      g.Networks(n).Modality = 'MR';
      g.Networks(n).Description = 'Trained to segment LV, and LA, in 3CH LAX images.';
      g.Networks(n).ImageFile = '';

      %LAX 4CH
      n = n + 1;
      g.Networks(n).Name = 'LAX 4CH';
      g.Networks(n).NV = {'5722.nv'};
      g.Networks(n).FileParts = 10;
      g.Networks(n).LicenseCode = ''; %non existing code
      g.Networks(n).PreInstalled = false;
      g.Networks(n).Visible = true;
      g.Networks(n).Modality = 'MR';
      g.Networks(n).Description = 'Trained to segment LV, RV, LA, RA in 4CH LAX images.';
      g.Networks(n).ImageFile = '';

    end
    %------------------------------------
    function checkinstallednetworks(g, waitforinput)
    %------------------------------------
    %Check if AI networks are downloaded

     if nargin < 2
        % waits for user to close the network wizard GUI manually, needed when starting directly from setupwizard
        waitforinput = false;
      end
    isadmin = iswindowsadmin; %User starts as admin first time after upgrade
    if isadmin || waitforinput
      type='all';
      networkwizard('checkavailability',type);
    end
    end

   %----------------------------------------
    function initconfigplaceholder(varargin)
   %---------------------------------------
      g = varargin{1};
      gHnd = g.Handles;
      gIcn = g.Icons.config;

      %Set flags to initialise myicon objects
      fdacleared = true;
      defaulttype = 1; %default
      defaultgroup = 1; %default
      defaultindent = []; %default

      %initcells
      lviconcell = cell(1,1);
      rviconcell = cell(1,1);
      laraiconcell = cell(1,1);
      analysisiconcell = cell(1,1);
      roiflowiconcell = cell(1,1);
      viabilityiconcell = cell(1,1);
      imageiconcell = cell(1,1);
      g.initfunctionicons();

      %Initialise icons that are reused in several tabs
      lvstackicon = myicon('lvstack',gHnd.configiconholder,gIcn.lvstack,'Go to LV stack',@() segment('viewspecial_Callback','lv'),0,defaultgroup,defaultindent,fdacleared);
      rvstackicon = myicon('rvstack',gHnd.configiconholder,gIcn.rvstack,'Go to RV stack',@() segment('viewspecial_Callback','rv'),0,defaultgroup,defaultindent,fdacleared);
      selectoneallicon = myicon('selectoneall',gHnd.configiconholder,gIcn.selectoneall,'Select one or all frames mode [1]',@()segment_main('singleframemode_Callback'),2,1,gIcn.selectallone,fdacleared);
      
      g.Icons.lviconcell = lviconcell;
      %RV
      g.Icons.rviconcell = rviconcell;
      %LA/RA 
       g.Icons.laraiconcell = laraiconcell;

%       %General pen
%       generalpeniconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'),defaulttype,defaultgroup,defaultindent,fdacleared);
%       generalpeniconcell{1,end+1} = myicon('contrastbrightness',gHnd.configiconholder,gIcn.contrastbrightness,'Manually change contrast and brightness',@() buttondownfunctions('updatebuttondowns','Contrast'),defaulttype,defaultgroup,defaultindent,fdacleared);
%       generalpeniconcell{1,end+1} = myicon('autocontrast',gHnd.configiconholder,gIcn.autocontrast,'Set contrast and brightness to predefined values',@() segment('autocontrast_Callback'),0,defaultgroup,defaultindent,fdacleared);
%       generalpeniconcell{1,end+1} = selectoneallicon;
% 
%       generalpeniconcell{1,end+1} = myicon('addgeneralobject',gHnd.configiconholder,gIcn.addgeneralobject,'Add a new object',@() generalpen.generalpenfunctions('createnewobject_Callback'),0,defaultgroup,defaultindent,fdacleared);
%       generalpeniconcell{1,end+1} = myicon('generalpen',gHnd.configiconholder,gIcn.generalpen,'General object pen',@() buttondownfunctions('updatebuttondowns','GeneralPen'),defaulttype,defaultgroup,defaultindent,fdacleared);
% %       generalpeniconcell{1,end+1} = myicon('interpgeneralpen',gHnd.configiconholder,gIcn.interpgeneralpen,'Set interpolation points for general object',@() buttondownfunctions('updatebuttondowns','GeneralPenInterp'),defaulttype,defaultgroup,defaultindent,fdacleared);
% %       generalpeniconcell{1,end+1} = myicon('hidegeneralpen',gHnd.configiconholder,gIcn.hidegeneralpen,'Hide selected object', @() viewfunctions('hide_Callback'),2,defaultgroup,defaultindent,fdacleared);
%       generalpeniconcell{1,end+1} = myicon('hidegeneralpenall',gHnd.configiconholder,gIcn.hidegeneralpenall,'Hide all objects', @() viewfunctions('hide_Callback'),2,defaultgroup,defaultindent,fdacleared);
% %       generalpeniconcell{1,end+1} = generaliconcell{1,1}; %hide interp
%       generalpeniconcell{1,end+1} = myicon('cleargeneralpen',gHnd.configiconholder,gIcn.cleargeneralpen,'Clear selected object',@() generalpen.generalpenfunctions('deleteobject_Callback'),0,defaultgroup,defaultindent,fdacleared);
%       generalpeniconcell{1,end+1} = myicon('cleargeneralpenall',gHnd.configiconholder,gIcn.cleargeneralpenall,'Clear all objects',@() generalpen.generalpenfunctions('deleteallobjects_Callback'),0,defaultgroup,defaultindent,fdacleared);
% %       generalpeniconcell{1,end+1} = myicon('volumecurve',gHnd.configiconholder,gIcn.volumecurve,'Plot volume curve',@() lvrvtools('plotvolumecurvegeneralpen'),0,defaultgroup,defaultindent,fdacleared);
%       generalpeniconcell{1,end+1} = myicon('label',gHnd.configiconholder,gIcn.text,'Set object''s label',@() generalpen.generalpenfunctions('setlabel_Callback'),0,defaultgroup,defaultindent,fdacleared);
% 
%       verticalline{1,1} = myicon('verticalline',gHnd.configiconholder,gIcn.verticalline,'',@()'',0);
%       generalpeniconcell{1,end+1} = verticalline{1,1};
% 
%       g.Icons.generalpeniconcell = generalpeniconcell;

      %ROIFLOW
      % check for Neusoft license
      isneusoft = g.isneusoft;

      roiflowiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'),defaulttype,defaultgroup,defaultindent,fdacleared);
      roiflowiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'),defaulttype,defaultgroup,defaultindent,fdacleared);
      roiflowiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'),defaulttype,defaultgroup,defaultindent,fdacleared);
      %roiflowiconcell{1,end+1} = myicon('scaleROI',gHnd.configiconholder,gIcn.scaleROI,'Scale ROI',@() updatetool('scaleROI'));
      
      [icontype0flow4d, ~, ~, enableflow4D] = checkmodulelicense('d');
      if strcmpi(enableflow4D,'on')
        roiflowiconcell{1,end+1} = myicon('flow4d',gHnd.configiconholder,gIcn.fourdflow,'4D Flow',@() flow4d.flow4d,icontype0flow4d);
      end
      flowstackicon = myicon('flowstack',gHnd.configiconholder,gIcn.flowstack,'Go to flow stack',@() segment('viewspecial_Callback','flow'),0,defaultgroup,defaultindent,fdacleared);
      roiflowiconcell{1,end+1} = flowstackicon;
      roiflowiconcell{1,end+1} = selectoneallicon;

      roiflowiconcell{1,end+1} = myicon('putroi',gHnd.configiconholder,gIcn.putroi,'Place ROI',@() buttondownfunctions('updatebuttondowns','RoiPut'),defaulttype,defaultgroup,defaultindent,fdacleared);
      roiflowiconcell{1,end+1} = myicon('roipen',gHnd.configiconholder,gIcn.roipen,'ROI pen',@() buttondownfunctions('updatebuttondowns','Roi'),defaulttype,defaultgroup,defaultindent,fdacleared);
      roiflowiconcell{1,end+1} = myicon('balloonroi',gHnd.configiconholder,gIcn.balloonroi,'Semi automatic ROI tool',@() buttondownfunctions('updatebuttondowns','RoiBalloon'),defaulttype,defaultgroup,defaultindent,fdacleared);
      if ~isneusoft
        roiflowiconcell{1,end+1} = myicon('trackingvessel',gHnd.configiconholder,gIcn.trackingvessel,'Track vessel over time [Alt-T]',@() vesselsnake_flowtrackroi('flowtrackroi'),0,defaultgroup,defaultindent,fdacleared);
      end
      roiflowiconcell{1,end+1} = myicon('refineroi',gHnd.configiconholder,gIcn.refineroi,'Refine Flow ROI [Alt-R]',@() flow('flowrefine_Callback'),0,defaultgroup,defaultindent,fdacleared);
      roiflowiconcell{1,end+1} = myicon('refineroinext',gHnd.configiconholder,gIcn.refineroinext,'Propagate Flow ROI forward and refine [Alt-F]',@() flow('flowpropagate_Callback'),0,defaultgroup,defaultindent,fdacleared);
      roiflowiconcell{1,end+1} = myicon('contractroi',gHnd.configiconholder,gIcn.contractroi,'Contract ROI',@() roi('expandcontract_Callback',-1),0,defaultgroup,defaultindent,fdacleared);
      roiflowiconcell{1,end+1} = myicon('expandroi',gHnd.configiconholder,gIcn.expandroi,'Expand ROI',@() roi('expandcontract_Callback',1),0,defaultgroup,defaultindent,fdacleared);
      if ~isneusoft
        roiflowiconcell{1,end+1} = myicon('unwrap',gHnd.configiconholder,gIcn.unwrap,'Unwrap flow',@() flowunwrap,0,defaultgroup,defaultindent,fdacleared);
      end
      roiflowiconcell{1,end+1} = myicon('palette',gHnd.configiconholder,gIcn.palette,'Set ROI color',@() roi('roisetcolor_Callback'),0,defaultgroup,defaultindent,fdacleared);
      roiflowiconcell{1,end+1} = myicon('text',gHnd.configiconholder,gIcn.text,'Set ROI label',@() roi('roisetlabel_Callback'),0,defaultgroup,defaultindent,fdacleared);
      if ~isneusoft
        roiflowiconcell{1,end+1} = myicon('plotflow',gHnd.configiconholder,gIcn.plotflow,'Plot flow [Ctrl-T]',@() reportflow,0,defaultgroup,defaultindent,fdacleared);
      end
      roiflowiconcell{1,end+1} = myicon('hideroi',gHnd.configiconholder,gIcn.hideroi,'Hide ROI',@() viewfunctions('hide_Callback'),2,defaultgroup,defaultindent,fdacleared);
      indhideroi = length(roiflowiconcell);
      roiflowiconcell{1,end+1} = myicon('hidetext',gHnd.configiconholder,gIcn.hidetext,'Hide text',@() viewfunctions('hide_Callback'),2,defaultgroup,defaultindent,fdacleared);
      indhidetext = length(roiflowiconcell);
      roiflowiconcell{1,end+1} = myicon('clearroi',gHnd.configiconholder,gIcn.clearroi,'Clear selected ROIs',@() roi('roidelete_Callback'),0,defaultgroup,defaultindent,fdacleared);
      roiflowiconcell{1,end+1} = myicon('clearallroi',gHnd.configiconholder,gIcn.clearallroi,'Clear all ROIs',@() roi('roiclearall_Callback'),0,defaultgroup,defaultindent,fdacleared);
      if strcmpi(enableflow4D,'on')
        roiflowiconcell{1,end+1} = myicon('clearfourdflow',gHnd.configiconholder,gIcn.clearfourdflow,'Clear 4D Flow',@() flow4d.flow4d('clearall_Callback'),icontype0flow4d);
      end
      g.Icons.roiflowiconcell = roiflowiconcell;

      %Viablility
      viabilityiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'));
      viabilityiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'));
      viabilityiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'));
      scarstackicon = myicon('scarstack',gHnd.configiconholder,gIcn.scarstack,'Go to scar stack',@() segment('viewspecial_Callback','cinescar'),0);
      viabilityiconcell{1,end+1} = scarstackicon;
      viabilityiconcell{1,end+1} = selectoneallicon;
      onlyEDES = false;
      onlycurrenttimeframe = true;
      iconautolvone = myicon('autolvone',gHnd.configiconholder,gIcn.autolvone,gettoolinfo('autolvone'),@() lvsegmentationml('doaltlv2_Callback',onlyEDES,onlycurrenttimeframe),0);
      
      viabilityiconcell{1,end+1} = iconautolvone;
      viabilityiconcell{1,end+1} = myicon('lgelv',gHnd.configiconholder,gIcn.autolvwand,gettoolinfo('autolvlge'),@() lgeml('lgesegmentationml'),0,defaultgroup,defaultindent,fdacleared);
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
      viabilityiconcell{1,end+1} = myicon('rubbermo',gHnd.configiconholder,gIcn.rubbermo,'Manually remove MO segmentation',@() buttondownfunctions('updatebuttondowns','MORubber'));      
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

      %% Strain 
      g.initstrainicons();
      
      %% Analysis
      analysisiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'),defaulttype,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'),defaulttype,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'),defaulttype,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('perfusionstack',gHnd.configiconholder,gIcn.perfusionstack,'Go to perfusion stacks',@() segment('viewspecial_Callback','perfusion'),0);
      analysisiconcell{1,end+1} = selectoneallicon;
      analysisiconcell{1,end+1} = iconautolvone;
      analysisiconcell{1,end+1} = myicon('importfromother',gHnd.configiconholder,gIcn.importfromother,'Import LV segmentation from other image stack with snap',@() segmentation('importsegmentationwithsnap_Callback'),0,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('measure',gHnd.configiconholder,gIcn.measure,'Place measurement',@() buttondownfunctions('updatebuttondowns','Measure'),defaulttype,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('point',gHnd.configiconholder,gIcn.point,'Place annotation point',@() buttondownfunctions('updatebuttondowns','Point'),defaulttype,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('roipen',gHnd.configiconholder,gIcn.roipen,'ROI pen',@() buttondownfunctions('updatebuttondowns','Roi'),defaulttype,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('addroiinlv',gHnd.configiconholder,gIcn.addroiinlv,'Add ROIs to sector of LV wall in selected slices',@() roi('roiaddinsector_Callback'),0,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('bullseye',gHnd.configiconholder,gIcn.bullseye,'Bullseye analysis [Ctrl-B]',@() reportbullseye,0,defaultgroup,defaultindent,fdacleared);
      if ~g.isneusoft
        analysisiconcell{1,end+1} = myicon('AVPD',gHnd.configiconholder,gIcn.AVPD,'AV plane displacement',@() avplane,0);
      end
      % check perfusion
      [icontype0perf, ~, ~, ~] = checkmodulelicense('p');
      if g.isneusoft
        % for Neusoft license pefusion scoring is the same as for perfusion
        icontype0perfscoring = icontype0perf;
      else
        % in Segment research perfusion scoring is available w/o license
        icontype0perfscoring = 0;
      end
      analysisiconcell{1,end+1} = myicon('perfusion',gHnd.configiconholder,gIcn.perfusion, 'Perfusion analysis',@() perfusion.perfusion,icontype0perf);
      analysisiconcell{1,end+1} = myicon('perfusionscoring',gHnd.configiconholder,gIcn.perfusionscoring, 'Perfusion scoring',@() perfusion.perfusionscoring,icontype0perfscoring);
      analysisiconcell{1,end+1} = myicon('reportperslice',gHnd.configiconholder,gIcn.reportperslice, 'Report per slice',@() slicereport,0,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('model3d',gHnd.configiconholder,gIcn.model3d, 'Show 3D model',@() report3dmodel,0,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('verticalline',gHnd.configiconholder,gIcn.verticalline,'',@()'',0);
      analysisiconcell{1,end+1} = myicon('hidemeasure',gHnd.configiconholder,gIcn.hidemeasure,'Hide measurements',@() viewfunctions('hide_Callback'),2,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('hidepoint',gHnd.configiconholder,gIcn.hidepoint,'Hide annotation points',@() viewfunctions('hide_Callback'),2,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = roiflowiconcell{1,indhideroi};
      analysisiconcell{1,end+1} = myicon('textbackground',gHnd.configiconholder,gIcn.textbackground,'Add text background',@() viewfunctions('addtextbackground'),2);
      ind=length(analysisiconcell);
      analysisiconcell{1,end+1} = myicon('verticalline',gHnd.configiconholder,gIcn.verticalline,'',@()'',0);
      analysisiconcell{1,end+1} = myicon('clearmeasure',gHnd.configiconholder,gIcn.clearmeasure,'Clear all measurements',@() callbackfunctions('measureclearall_Callback'),0,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('clearpoint',gHnd.configiconholder,gIcn.clearpoint,'Clear all annotation points',@() callbackfunctions('pointclearall_Callback'),0,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('clearroi',gHnd.configiconholder,gIcn.clearroi,'Clear selected ROIs',@() roi('roidelete_Callback'),0,defaultgroup,defaultindent,fdacleared);
      analysisiconcell{1,end+1} = myicon('clearallroi',gHnd.configiconholder,gIcn.clearallroi,'Clear all ROIs',@() roi('roiclearall_Callback'),0,defaultgroup,defaultindent,fdacleared);
      g.Icons.analysisiconcell = analysisiconcell;
      %update indent status for  textbackground icon
      g.Icons.analysisiconcell{ind}.isindented = g.Pref.BackgroundColor;
      %% Image
      imageiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'),defaulttype,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'),defaulttype,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'),defaulttype,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('autocontrastall',gHnd.configiconholder,gIcn.autocontrastall,'Set contrast and brightness to predefined values for all images',@() segment('autocontrastall_Callback'),0,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('resetlightall',gHnd.configiconholder,gIcn.resetlightall,'Reset contrast and brightness for all imagestacks',@() segment('resetlightall_Callback'),0,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = selectoneallicon;

      imageiconcell{1,end+1} = myicon('crop',gHnd.configiconholder,gIcn.crop,'Manual crop',@() buttondownfunctions('updatebuttondowns','crop'),defaulttype,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('cropall',gHnd.configiconholder,gIcn.cropall,'Auto crop all',@() autocropallgui('init'),0,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('croplv',gHnd.configiconholder,gIcn.croplv,'LV crop',@() lvsegmentation('croplvall',0),0,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('cineplay',gHnd.configiconholder,gIcn.cineplay,'Open cine tool',@() segment('cinetool_Callback'),2,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('movie',gHnd.configiconholder,gIcn.movie,'Open movie tool',@() export('exportmovierecorder_Callback'),0,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('click3d',gHnd.configiconholder,gIcn.click3d,'Set 3D point',@() buttondownfunctions('updatebuttondowns','click3d'),defaulttype,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('rotate90',gHnd.configiconholder,gIcn.rotate90,'Rotate 90 degrees clockwise',@() tools('rotate90right_Callback'),0,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('mpr',gHnd.configiconholder,gIcn.mpr,'Reconstruct image stack',@() reformater,0,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('mergestacks',gHnd.configiconholder,gIcn.mergestacks,'Merge stacks',@() mergestacks,0,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('imageinfo',gHnd.configiconholder,gIcn.imageinfo,'View image info',@() imageinfo('init'),0,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('patientinfo',gHnd.configiconholder,gIcn.patientinfo,'View and adjust patient info',@() tools('viewpatientinfo_Callback'),0,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('hideid',gHnd.configiconholder,gIcn.hideid,'Blind Subject Identity',@() tools('hideid_Callback'),2,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = roiflowiconcell{1,indhidetext};
      imageiconcell{1,end+1} = myicon('hideplus',gHnd.configiconholder,gIcn.hideplus,'Hide center cross',@() viewfunctions('hide_Callback'),2,defaultgroup,defaultindent,fdacleared);
      imageiconcell{1,end+1} = myicon('hideintersections',gHnd.configiconholder,gIcn.hideintersections,'Hide intersection lines',@() viewfunctions('hide_Callback'),2,defaultgroup,defaultindent,fdacleared);
      g.Icons.imageiconcell = imageiconcell;

      %% TXmap
      txmapiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'));
      txmapiconcell{1,end+1} = myicon('move',gHnd.configiconholder,gIcn.move,'Translate contour',@() buttondownfunctions('updatebuttondowns','move'));
      txmapiconcell{1,end+1} = myicon('scale',gHnd.configiconholder,gIcn.scale,'Scale object',@() buttondownfunctions('updatebuttondowns','scale'));
      txmapiconcell{1,end+1} = selectoneallicon;
      txmapiconcell{1,end+1} = iconautolvone;
      txmapiconcell{1,end+1} = myicon('importcorrect',gHnd.configiconholder,gIcn.importcorrect,'Import LV segmentation from cine to Tx map image stack with position correction',@() segmentation('importfromcine2txmap_Callback',true),0);
      txmapiconcell{1,end+1} = myicon('importfromother',gHnd.configiconholder,gIcn.importfromother,'Import LV segmentation from cine to Tx map image stack without postion correction',@() segmentation('importfromcine2txmap_Callback',false),0);
      txmapiconcell{1,end+1} = myicon('copybackandforward',gHnd.configiconholder,gIcn.copybackandforward,'Copy LV segmentation from current time frame to all time frames',@() tools('copytoalltimeframes_Callback','lv','true','true'),0);
      txmapiconcell{1,end+1} = myicon('moveall',gHnd.configiconholder,gIcn.moveall,'Translate all contours',@() buttondownfunctions('updatebuttondowns','moveall'));
      txmapiconcell{1,end+1} = myicon('endopen',gHnd.configiconholder,gIcn.endopen,'Endo pen',@() buttondownfunctions('updatebuttondowns','Endo'));
      txmapiconcell{1,end+1} = myicon('epipen',gHnd.configiconholder,gIcn.epipen,'Epi pen',@() buttondownfunctions('updatebuttondowns','Epi'));
      txmapiconcell{1,end+1} = myicon('interpendo',gHnd.configiconholder,gIcn.interpendo,'Set interpolation points for Endo',@() buttondownfunctions('updatebuttondowns','EndoInterp'));
      txmapiconcell{1,end+1} = myicon('interpepi',gHnd.configiconholder,gIcn.interpepi,'Set interpolation points for Epi',@() buttondownfunctions('updatebuttondowns','EpiInterp'));
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

      reviewiconcell{1,1} = myicon('select',gHnd.configiconholder,gIcn.select,'Select image stack or object',@() buttondownfunctions('updatebuttondowns','select'),1,defaultgroup,defaultindent,fdacleared);
      reviewiconcell{1,end+1} = lvstackicon;
      reviewiconcell{1,end+1} = rvstackicon;
      reviewiconcell{1,end+1} = flowstackicon;
      reviewiconcell{1,end+1} = scarstackicon;
      [icontype0mitt, ~, ~, ~] = checkmodulelicense('M');
      strainsaxstackicon = myicon('strainsaxstack',gHnd.configiconholder,gIcn.strainsaxstack,'Go to Strain SAX stack',@() segment('viewspecial_Callback','strainsax'),0);
      strainlaxstackicon = myicon('strainlaxstack',gHnd.configiconholder,gIcn.strainlaxstack,'Go to Strain LAX stack',@() segment('viewspecial_Callback','strainlax'),0);
      strainsaxanalysis = myicon('strainmittftsaxanalysis',gHnd.configiconholder,gIcn.strainmittftsaxanalysis,'Start FT Strain SAX analysis',@() strainmitt.strainmitt('init','cine','shortaxis'),icontype0mitt);
      strainlaxanalysis = myicon('strainmittftlaxanalysis',gHnd.configiconholder,gIcn.strainmittftlaxanalysis,'Start FT Strain LAX analysis',@() strainmitt.strainmitt('init','cine','longaxis'),icontype0mitt);
      reviewiconcell{1,end+1} = strainsaxstackicon;
      reviewiconcell{1,end+1} = strainsaxanalysis;
      reviewiconcell{1,end+1} = strainlaxstackicon;
      reviewiconcell{1,end+1} = strainlaxanalysis;
      reviewiconcell{1,end+1} = myicon('automateinfo',gHnd.configiconholder,gIcn.imageinfoai,'View AI AutoMate notifications',@() warningfunctions('showautoloaderwarnings'),0,defaultgroup,defaultindent,fdacleared);


      g.Icons.reviewiconcell = reviewiconcell;
    end

    

    %----------------------------------------
    function initribbonplaceholder(varargin)
      %--------------------------------------

      g = varargin{1};
      gHndtoggle = g.Handles.toggleiconholder;
      gIcntoggle = g.Icons.toggleicons;
      %initiate iconholder for toggle axes
      fdacleared = true;
      ind = 1;

      iconCell{ind} = myicon('ribbonfunction',gHndtoggle, gIcntoggle.ribbonfunctionoff,...
        'Function', @() g.togglebuttonFunction_Callback,1,1,gIcntoggle.ribbonfunctionon,fdacleared);

%       iconCell{end+1} = myicon('ribbongeneralpen',gHndtoggle,gIcntoggle.ribbongeneralpenoff,...
%         'GeneralPen',@() g.togglebuttonGeneralPen_Callback,1,1,gIcntoggle.ribbongeneralpenon,fdacleared);
      iconCell{end+1} = myicon('ribbonflow',gHndtoggle,gIcntoggle.ribbonflowoff,...
        'ROI/FLOW',@() g.togglebuttonROIFLOW_Callback,1,1,gIcntoggle.ribbonflowon,fdacleared);
      if ~g.isneusoft
        iconCell{end+1} = myicon('ribbonviability',gHndtoggle,gIcntoggle.ribbonviabilityoff,...
          'Viability',@() g.togglebuttonVia_Callback,1,1,gIcntoggle.ribbonviabilityon);
        iconCell{end+1} = myicon('ribbontxmap',gHndtoggle,gIcntoggle.ribbontxmapoff,...
          'TX maps and ECV',@() g.togglebuttontxmap_Callback,1,1,gIcntoggle.ribbontxmapon);
      else
        % ensure specific menus are invisible
        set([ ...
          g.Handles.mrmenu, ...
          g.Handles.ctmenu, ...
          g.Handles.spectmenu, ...
          g.Handles.avplanemenu, ...
          ], ...
          'Visible','off');
        
      end
      iconCell{end+1} = myicon('ribbonstrain',gHndtoggle,gIcntoggle.ribbonstrainoff,...
        'Strain',@() g.togglebuttonStrain_Callback,1,1,gIcntoggle.ribbonstrainon);
      iconCell{end+1} = myicon('ribbonanalysis',gHndtoggle,gIcntoggle.ribbonanalysisoff,...
        'Analysis',@() g.togglebuttonAnalysis_Callback,1,1,gIcntoggle.ribbonanalysison,fdacleared);
      iconCell{end+1} = myicon('ribbonimage',gHndtoggle,gIcntoggle.ribbonimageoff,...
        'Image',@() g.togglebuttonImage_Callback,1,1,gIcntoggle.ribbonimageon,fdacleared);
      iconCell{end+1} = myicon('ribbonreview',gHndtoggle,gIcntoggle.ribbonreviewoff,...
        'Review',@() g.togglebuttonReview_Callback,1,1,gIcntoggle.ribbonreviewon,fdacleared);

      
      gHndtoggle.add(iconCell);
      pos = plotboxpos(gHndtoggle.axeshandle);
      currentpos = get(gHndtoggle.axeshandle,'position');
      set(gHndtoggle.axeshandle,'position',currentpos-[pos(1),0,0,0]);
    end

    %----------------------
    function inittoolbar(g)
      %----------------------
      g.GUISettings.TopGapHeight = 0.135;
      g.GUISettings.BottomGapHeight = 0.133;
      g.GUISettings.RightGapWidth = 0.21;
      g.GUISettings.ShowColorbar =  false;

      %--- Initialize menu items under Utility and Plugins
      segment('initmenu');
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
          g.Handles.lvheadertext.String{6} = 'CO (L/min)';
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
      %check license
      if g.isneusoft
        return
      end
      no = NO;
      roitxt = '';
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

      roinbr = [];
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

        roitxt = SET(no).Roi(SET(no).RoiCurrent).Name;

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
          SET(no).Flow.HeartRate = temphr;
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
        roitxt = SET(no).Roi(roinbr).Name;
      end

      c_inds = cellfun(@(x) isempty(x) || isnan(x) || x == 0,c);
      [c{c_inds}] = deal('---');

      c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
      set(g.Handles.flowheaderroitext,'String',roitxt);
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
      if g.isneusoft
        set(g.Handles.flowuipanel, 'Visible', 'off');
        return
      end
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
        hasla = ismember(no,findfunctions('findlaxnowithla')) || ismember(no,findfunctions('findlaxnowithlaest'));
        hasra = ismember(no,findfunctions('findlaxnowithra')) || ismember(no,findfunctions('findlaxnowithraest'));
      else
        hasscar = false;
        hasmar = false;
        hasatrialscar = false;
        hasla = false;
        hasra = false;
      end

      if any([hasscar hasmar hasatrialscar hasla hasra])
        g.Handles.measurementuipanel.Visible = 'on';
        g.Handles.flowuipanel.Visible = 'off';
      else
        g.Handles.measurementuipanel.Visible = 'off';
        if ~g.isneusoft
          g.Handles.flowuipanel.Visible = 'on';
        end
      end

      if hasscar
        scarml = calcfunctions('calcscarml', no);
      end

      if hasmar && hasscar
        c = cell(1,7);
        c{1} = round(SET(no).LVM *1.05); %g
        c{2} = round(scarml, 1);
        c{3} = round(scarml * 1.05, 1); %g
        c{4} = round(SET(no).Scar.Percentage, 1);
        c{5} = round(SET(no).Scar.MOPercentage, 1);
        c{6} = (round(SET(no).MaR.Percentage(SET(no).EDT)));
        c{7} = (round(SET(no).MaR.Percentage(SET(no).EST)));

        c_inds = cellfun(@(x) isempty(x) || isnan(x),c);
        [c{c_inds}] = deal('---');
        c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
        g.Handles.measurementheadertext.String = sprintf('%s\n%s\n%s\n%s\n%s\n%s','LVM (g)', 'Scar (ml)', 'Scar (g)','Scar (%)','MO (%)','ED/ES MaR (%)');
        g.Handles.measurementreporttext.String = sprintf('%s\n%s\n%s\n%s\n%s/%s\n%s',c{:});

      elseif hasscar
        c = cell(1,5);
        c{1} = round(SET(no).LVM *1.05); %g
        c{2} = round(scarml, 1);
        c{3} = round(scarml * 1.05, 1); %g
        c{4} = round(SET(no).Scar.Percentage, 1);
        c{5} = round(SET(no).Scar.MOPercentage, 1);

        c_inds = cellfun(@(x) isempty(x) || isnan(x),c);
        [c{c_inds}] = deal('---');
        c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
        g.Handles.measurementheadertext.String = sprintf('%s\n%s\n%s\n%s\n%s','LVM (g)', 'Scar (ml)', 'Scar (g)', 'Scar (%)','MO (%)');
        g.Handles.measurementreporttext.String = sprintf('%s\n%s\n%s\n%s\n%s',c{:});

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
      
      elseif hasla || hasra
        c = cell(1,2);
        [laedv,laesv,laef,laev,no2ch,no4ch] = calcfunctions('calclavalues');

        c{1} = round(laedv);
        c{2} = round(laesv);
        c{3} = round(laev);
        c{4} = round(laef);
        c{5} = no2ch;
        c{6} = no4ch;

        [raeda,raesa,nora] = calcfunctions('calcravalues');

        c{7} = round(raeda);
        c{8} = round(raesa);
        c{9} = nora;
        c_inds = cellfun(@(x) isempty(x) || isnan(x),c);
        [c{c_inds}] = deal('---');
        c = cellfun(@(x) num2str(x),c,'UniformOutput',false);
        raedastring = ['RA EDA (cm',char(178),')'];
        raesastring = ['RA ESA (cm',char(178),')'];
        g.Handles.measurementheadertext.String = sprintf('%s\n%s\n%s\n%s\n%s / %s\n\n%s\n%s\n%s', ...
          'LA EDV (ml)','LA ESV (ml)', ...
          'LA EV (ml)','LA EF (%)', ...
          '#2CH','#4CH', ...
          raedastring,raesastring, ...
          '#4CH');
        g.Handles.measurementreporttext.String = sprintf('%s\n%s\n%s\n%s\n%s / %s\n\n%s\n%s\n%s',c{:});

        %adjust fontsize according to the available space on screen
        pos = getpixelposition(g.Handles.measurementheadertext);
        widthinpixel = pos(3);
        cmperwidth  = 0.040*widthinpixel/10.5;
        cmperlength = 0.0265*pos(4)/9.75;
        g.fontsizeincm = min(cmperwidth,cmperlength);
      end

      if hasla || hasra
        title = 'LA/RA';
      else
        title = 'Viability';
      end

      viability('sliderupdate')    
      set(g.Handles.measurementtitletext,'String',title);
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

    %----------------------
    function initfdaversion(g)
      %----------------------
      % Function to check if the license and the user preferences allow for
      % FDA version or not. Calls callbacks functions to set the menus' 
      % visiblity accordingly. edit segmentmakemat

      if g.isFDAVersion
        %Make buttons visible
        fdabuttongroup = g.Handles.fdabuttongroup;
        set(fdabuttongroup,'Visible','on'); %default is off

        %Initialise the radio buttons
        fdabutton = g.Handles.fdaradiobutton;
        researchbutton = g.Handles.researchradiobutton;
        set(fdabutton,'String','FDA cleared version','Callback','callbackfunctions(''fdabutton_Callback'')');
        set(researchbutton,'String','Research version','Callback','callbackfunctions(''researchbutton_Callback'')');

        %Adapt radiobuttons length
        g.resizebuttongroup(fdabuttongroup);

        %Check user's preferences for running FDA or research version
        if g.Pref.RunFDAVersion
          set(fdabutton,'Value',true);
          callbackfunctions('fdabutton_Callback');
        else
          set(researchbutton,'Value',true);
          callbackfunctions('researchbutton_Callback');
        end
      end
    end

    %----------------------------
    function updatetimebaraxes(g)
      %----------------------------
      %overloads maingui method to show/hide play buttons
      global SET NO

      updatetimebaraxes@maingui(g);
      no = NO;

      if isempty(no) || SET(no).TSize < 2
        set(g.Handles.playuipanel, 'Visible', 'off');
      else
        set(g.Handles.playuipanel,'Visible','on');
        g.updatefps(no);
      end
    end

  end

  methods(Static)

    %--------------------
    function resizebuttongroup(buttongroup)
      %--------------------
      %Helper function to resize the radiobuttons and buttongroup based on
      %their extent

      buttons = num2cell(flipud(buttongroup.Children));
      set([buttons{:} buttongroup],'units','pixels');
      bufferspace = 20; %pixels
      nbuttons = numel(buttons);
      leftbufferspace = 5;

      %Calculate total width of buttons
      buttonwidths = cellfun(@(buttons) get(buttons,'Extent'),buttons,'UniformOutput',false);
      buttonwidths = cellfun(@(extent) extent(3), buttonwidths);
      totalwidth = sum(buttonwidths) + leftbufferspace + bufferspace*nbuttons;

      %Adjust buttons position
      buttongroupposition = get(buttongroup, 'Position');
      set(buttongroup,'Position',...
        [buttongroupposition(1) + buttongroupposition(3) - totalwidth, ...
        buttongroupposition(2), totalwidth, buttongroupposition(4)]);
      currentxpos = leftbufferspace;
      for loop = 1:nbuttons
        buttonposition = get(buttons{loop},'Position');
        set(buttons{loop},'Position',...
          [currentxpos, buttonposition(2), buttonwidths(loop) + bufferspace, buttonposition(4)]);
        currentxpos = currentxpos + buttonwidths(loop) + bufferspace;
      end

      %Revert to original units
      set([buttons{:} buttongroup],'units','normalized');
    end

    %--------------------
    function setupwelcomewindow
      %--------------------
      %Method to display the software welcome window
      versionhello('init');
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

    %------------------------------------------------
    function menuoptions = getmultiplepatientsoptions
      %------------------------------------------------
      %Return options for a user trying to load data from patient
      %with different ID.
      menuoptions = {'Abort loading of new image stacks',...
        'Close previously loaded image stacks',...
        'Load all image stacks (not recommended)'};
    end

    %------------------------------------------------
    function [ok, msg] = loadingerrorgadgetron(ignoreddata)
      %------------------------------------------------
      msg = sprintf('Warning for series %s: %s',ignoreddata(1).series,ignoreddata(1).msg);
      ok = true;
    end

  end

end