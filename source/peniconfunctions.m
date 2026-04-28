classdef peniconfunctions %skiptranslation
  % class containing functions to handle pen icons (LV/RV/LA/RA)
  methods (Static)

    %-----------------------
    function penbuttondown()
      %---------------------
      % draw with currently chosen pen
      pen = peniconfunctions.getcurrentpen();
      buttondownfunctions('updatebuttondowns',pen);
    end

    %-----------------------------
    function peninterpbuttondown()
      %---------------------------
      % set interpolation points with currently chosen pen
      pen = peniconfunctions.getcurrentpen();
      currenttool = [pen,'Interp'];
      buttondownfunctions('updatebuttondowns',currenttool)
    end

    %------------------------------
    function penballoonbuttondown()
      %----------------------------
      % draw balloon with currently chosen pen (only for LV endo/RV endo)
      pen = peniconfunctions.getcurrentpen();
      currenttool = [pen,'Balloon'];
      buttondownfunctions('updatebuttondowns',currenttool)
    end

    %------------------
    function refine()
      %----------------
      % contract contour based on currently chosen pen
      pen = peniconfunctions.getcurrentpen();
      switch pen
        case 'Endo'
          lv('segmentrefineendo_Callback');
        case 'Epi'
          lv('segmentrefineepi_Callback');
        case 'RVEndo'
          rv('segmentrefinervendo_Callback');
        case 'RVEpi'
          % not implemented
        case 'LA'
          % not implemented
        case 'RA'
          % not implemented
      end
    end

    %------------------
    function contract()
      %----------------
      % contract contour based on currently chosen pen
      pen = peniconfunctions.getcurrentpen();
      switch pen
        case 'Endo'
          lv('segmentexpandcontract_Callback',-1,'endo')
        case 'Epi'
          lv('segmentexpandcontract_Callback',-1,'epi')
        case 'RVEndo'
          rv('expandcontract_Callback',-1,'rvendo')
        case 'RVEpi'
          rv('expandcontract_Callback',-1,'rvepi')
        case 'LA'
          % not implemented yet
        case 'RA'
          % not implemented yet
      end
    end

    %----------------
    function expand()
      %--------------
      % expand contour based on currently chosen pen
      pen = peniconfunctions.getcurrentpen();
      switch pen
        case 'Endo'
          lv('segmentexpandcontract_Callback',1,'endo')
        case 'Epi'
          lv('segmentexpandcontract_Callback',1,'epi')
        case 'RVEndo'
          rv('expandcontract_Callback',1,'rvendo')
        case 'RVEpi'
          rv('expandcontract_Callback',1,'rvepi')
        case 'LA'
          % not implemented yet
        case 'RA'
          % not implemented yet
      end
    end

    %---------------
    function clear()
      %-------------
      % clear contour based on currently chosen pen (LV and RV clear all
      % corresponding contours)
      pen = peniconfunctions.getcurrentpen();
      switch pen
        case {'Endo','Epi'}
          callbackfunctions('segmentclearalllv_Callback');
        case {'RVEndo','RVEpi'}
          callbackfunctions('segmentclearallrv_Callback');
        case 'LA'
          generalpen.atriumpenfunctions('deleteobject_Callback','la');
        case 'RA'
          generalpen.atriumpenfunctions('deleteobject_Callback','ra');
      end
    end

    %--------------------
    function clearslice()
      %------------------
      % clear contour in current slice based on currently chosen pen
      pen = peniconfunctions.getcurrentpen();
      switch pen
        case 'Endo'
          segmentation('clearslicesthis_Callback',1,0,0,0);
        case 'Epi'
          segmentation('clearslicesthis_Callback',0,1,0,0);
        case 'RVEndo'
          segmentation('clearslicesthis_Callback',0,0,1,0);
        case 'RVEpi'
          segmentation('clearslicesthis_Callback',0,0,0,1);
        case 'LA'
          generalpen.atriumpenfunctions('clearslice_Callback','la');
        case 'RA'
          generalpen.atriumpenfunctions('clearslice_Callback','ra');
      end
    end

    %--------------------------
    function selectpen(penname)
      %------------------------
      % select pen with the provided pen name
      global DATA

      peniconfunctions.updateballoonicon(penname);
      peniconfunctions.updatesmoothicon(penname);
      peniconfunctions.updaterefineicon(penname);
      peniconfunctions.updateexpandicon(penname);
      peniconfunctions.updatecontracticon(penname);

      penmode = peniconfunctions.getcurrentpenmode;
      penname = [penname,penmode];

      DATA.CurrentTool = penname;
      buttondownfunctions('updatebuttondowns',penname);
    end

    %--------------
    function hide(heartpart)
      %------------
      % hide/show contour based on currently chosen pen and corresponding
      % hide status
      global DATA
      DATA.HideState.(upper(heartpart)) = ~ DATA.HideState.(upper(heartpart));
      viewfunctions('hide_Callback')
    end

    %-------------------------------------
    function [pen,penicon] = getcurrentpen
      %-----------------------------------
      % get currenlty chosen pen
      global DATA %#ok<*GVMIS>
      switch DATA.CurrentTheme
        case 'strain'
          currenticons = [DATA.Icons.strainiconcell{:}];
          penicons = {
            'endopen', 'Endo';
            'epipen', 'Epi';
            'rvendopen','RVEndo';
            'lapen', 'LA';
            'rapen','RA'
            };
        case 'function'
          currenticons = [DATA.Icons.functioniconcell{:}];
          penicons = {
            'endopen', 'Endo';
            'epipen', 'Epi';
            'rvendopen','RVEndo';
            'rvepipen', 'RVEpi';
            'lapen', 'LA';
            'rapen','RA'
            };
      end
      iconnames = {currenticons.name};

      [ispresent, indpen] = ismember(penicons(:,1),iconnames);
      indpen = indpen(ispresent);
      isindented = [currenticons(indpen).isindented];
      ind = find(isindented);
      if ~isempty(ind)
        pen = penicons{ind(1),2};
        penicon = penicons{ind(1),1};
      else
        pen = '';
        penicon = '';
      end
    end

    %--------------------------------------------------
    function [penmode, penmodeicon] = getcurrentpenmode
      %------------------------------------------------
      % get currently chosen pen mode (pen/interpolation/balloon)
      global DATA
      switch DATA.CurrentTheme
        case 'strain'
          currenticons = [DATA.Icons.strainiconcell{:}];
          penmodeicons = {
            'pen','';
            'interpolate','Interp';
            };
        case 'function'
          currenticons = [DATA.Icons.functioniconcell{:}];
          penmodeicons = {
            'pen','';
            'interpolate','Interp';
            'balloon','Balloon';
            };
      end
      iconnames = {currenticons.name};

      [ispresent, indpen] = ismember(penmodeicons(:,1),iconnames);
      indpen = indpen(ispresent);
      isindented = [currenticons(indpen).isindented];
      ind = find(isindented);
      if ~isempty(ind)
        penmode = penmodeicons{ind(1),2};
        penmodeicon = penmodeicons{ind(1),1};
      else
        penmode = '';
        penmodeicon = '';
      end
    end

    %-------------------------------
    function updatehideicon(penname)
      %-----------------------------
      % indents/undents hide icon basen on current pen hide status
      global DATA
      iconholder = DATA.Handles.configiconholder;
      switch penname
        case {'Endo','Epi'}
          hidestate = DATA.HideState.LV;
        case {'RVEndo','RVEpi'}
          hidestate = DATA.HideState.RV;
        case 'LA'
          hidestate = DATA.HideState.LA;
        case 'RA'
          hidestate = DATA.HideState.RA;
      end
      runicon = false;
      if hidestate
        iconholder.indent('hidesegmentation',runicon);
      else
        iconholder.undent('hidesegmentation',runicon);
      end
    end

    %----------------------------------
    function updateballoonicon(penname)
      %--------------------------------
      % enable/disable balloon icon based on pen
      global DATA
      iconholder = DATA.Handles.configiconholder;
      functionicons = [DATA.Icons.functioniconcell{:}];
      iconnames = {functionicons.name};
      [ispresent, indicon] = ismember('balloon',iconnames);
      indicon = indicon(ispresent);
      if ~isempty(indicon)
        balloonicon = functionicons(indicon);
        if strcmpi(penname,'endo') || strcmpi(penname,'rvendo')
          % check if chosen
          if ~balloonicon.enabled
            balloonicon.enable(DATA.Pref.RunFDAVersion);
          end
        else
          % check if chosen
          if balloonicon.isindented
            % indent normal pen
            runicon = true;
            iconholder.indent('pen',runicon);
          end
          balloonicon.disable();
        end
      end
    end

    %----------------------------------
    function updatesmoothicon(penname)
      %--------------------------------
      % enable/disable smooth icon based on pen
      global DATA

      if strcmp(DATA.CurrentTheme,'strain')
        functionicons = [DATA.Icons.strainiconcell{:}];
      else
        functionicons = [DATA.Icons.functioniconcell{:}];
      end
      iconnames = {functionicons.name};
      [ispresent, indicon] = ismember('smooth',iconnames);
      indicon = indicon(ispresent);
      if ~isempty(indicon)
        smoothicon = functionicons(indicon);
        if strcmpi(penname,'la') || strcmpi(penname,'ra')
          % check if chosen
          smoothicon.disable();
        else
          smoothicon.enable(DATA.Pref.RunFDAVersion);
        end
      end
    end

    %----------------------------------
    function updaterefineicon(penname)
      %--------------------------------
      % enable/disable refine icon based on pen
      global DATA
      functionicons = [DATA.Icons.functioniconcell{:}];
      iconnames = {functionicons.name};
      [ispresent, indicon] = ismember('refine',iconnames);
      indicon = indicon(ispresent);
      if ~isempty(indicon)
        refineicon = functionicons(indicon);
        if strcmpi(penname,'la') || strcmpi(penname,'ra') || strcmpi(penname,'RVEpi')
          % check if chosen
          refineicon.disable();
        else
          refineicon.enable(DATA.Pref.RunFDAVersion);
        end
      end
    end

    %----------------------------------
    function updateexpandicon(penname)
      %--------------------------------
      % enable/disable expand icon based on pen
      global DATA

      functionicons = [DATA.Icons.functioniconcell{:}];
      iconnames = {functionicons.name};
      [ispresent, indicon] = ismember('expand',iconnames);
      indicon = indicon(ispresent);
      if ~isempty(indicon)
        refineicon = functionicons(indicon);
        if strcmpi(penname,'la') || strcmpi(penname,'ra')
          % check if chosen
          refineicon.disable();
        else
          refineicon.enable(DATA.Pref.RunFDAVersion);
        end
      end
    end
    %----------------------------------
    function updatecontracticon(penname)
      %--------------------------------
      % enable/disable contract icon based on pen
      global DATA

      functionicons = [DATA.Icons.functioniconcell{:}];
      iconnames = {functionicons.name};
      [ispresent, indicon] = ismember('contract',iconnames);
      indicon = indicon(ispresent);
      if ~isempty(indicon)
        refineicon = functionicons(indicon);
        if strcmpi(penname,'la') || strcmpi(penname,'ra')
          % check if chosen
          refineicon.disable();
        else
          refineicon.enable(DATA.Pref.RunFDAVersion);
        end
      end
    end
  end
end