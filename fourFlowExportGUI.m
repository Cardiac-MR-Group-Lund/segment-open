function varargout = fourFlowExportGUI(varargin)
%FourFlow export GUI

%Christoffer Green

if nargin==0;
  init;
else
  macro_helper(varargin{:});  
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
end;

%------------
function init
%------------
%Initializes GUI
global DATA SET

if isempty(SET)
  myfailed('There are no image stacks to export');
  return
end

%Initialize GUI
gui = mygui('fourFlowExportGUI.fig');
DATA.GUI.FourFlowExport = gui; %Store to global variable.

%Find what stacks that are possible to export and create suitable names for them.
stacks = prepareforexport;

gui.stacks = stacks;

%gui.stacks(n).data.OrgTSize == time steps
%gui.stacks(n).data.TIncr == time increment

%Fill listbox
%cellstri = cat(1, stacks(:).name); %Convert the struct to a cell array. 1 is along first dimension.
cells = {};
for i = 1:length(stacks(:))
    cells{i} = stacks(i).name;
end;

set(gui.handles.stacksListbox,'value',1:length(stacks(:)));
set(gui.handles.stacksListbox,'string',cells);
set(gui.handles.stacksListbox,'Max',length(cells));

%---------------------------------------------
function stacks = prepareforexport
%---------------------------------------------
global SET

%%% Loop all stacks
velc = 1;
saxc = 1;
CH2c = 1;
CH3c = 1;
CH4c = 1;
RVOTc = 1;
QFlowc = 1;
anatc = 1;
novector = 1:length(SET);

stacks = [];
stacks.type = '';
stacks.name = '';

loopsize = length(SET);
pos = 1;
currentstack = 1;
cineno = findfunctions('findno');
while pos <= loopsize
  if isempty(SET(pos).Flow) && (SET(pos).TSize > 1) && (SET(pos).ZSize == 1)
    stacks(currentstack).type = 'anatomy_time_resolved.';
    stacks(currentstack).name = sprintf('%s_%d',stacks(currentstack).type,RVOTc); %#ok<AGROW>
    stacks(currentstack).nos = pos; %#ok<AGROW>
    stacks(currentstack).data = SET(pos); %#ok<AGROW>
    %Fix for 3D_volume RVOT Acquisition: OBS Does not support Segmentation
    %---------------------------------------------------
%     if stacks(currentstack).data.ZSize > 1
%       %1. Change Imageposition to current slice:
%       zdir = cross(SET(pos).ImageOrientation(1:3),SET(pos).ImageOrientation(4:6));
%       stacks(currentstack).data.ImagePosition = stacks(currentstack).data.ImagePosition(:)'-(SET(pos).CurrentSlice-1)*(SET(pos).SliceThickness+SET(pos).SliceGap)*zdir(:)';
%       
%       %2. Set slice parameters to a single slice
%       stacks(currentstack).data.ZSize = 1;
%       stacks(currentstack).data.StartSlice = 1;
%       stacks(currentstack).data.EndSlice = 1;
%       stacks(currentstack).data.SliceGap = 0;
%       stacks(currentstack).data.OrgZSize = 1;
%       
%       %3. Take current slice of volume
%       stacks(currentstack).data.IM = single(squeeze(SET(pos).IM(:,:,:,SET(pos).CurrentSlice)));
%       
%       %4. Set final slice parameter correctly
%       stacks(currentstack).data.CurrentSlice = 1;
%     end
    %---------------------------------------------------
	currentstackName = stacks(currentstack).name;
	RVOTc = RVOTc+1;
    currentstack = currentstack+1;
	
	%Segmentation/Mesh support for standard views in fourflow: 
	if not(isempty(SET(pos).EndoX))
        stacks(currentstack).type = 'mesh';
        stacks(currentstack).subtype = 'Endo';
        stacks(currentstack).name = sprintf('%s_Endo',currentstackName);
        stacks(currentstack).nos = pos;
        stacks(currentstack).data = SET(pos);
        stacks(currentstack).x = SET(pos).EndoX;
        stacks(currentstack).y = SET(pos).EndoY;
        stacks(currentstack).z = (0:size(stacks(currentstack).x,3)-1);
        stacks(currentstack).timesteps = SET(pos).TIncr*(0:(SET(pos).TSize-1));
        currentstack = currentstack+1;
	end
  elseif ((SET(pos).TSize == 1) && (SET(pos).ZSize==1))
	 stacks(currentstack).type = 'anatomy';
	 stacks(currentstack).name = sprintf('%s_%d',stacks(currentstack).type,anatc);
	 stacks(currentstack).nos = pos;
	 stacks(currentstack).data = SET(pos);
	 anatc = anatc+1;
	 currentstack = currentstack+1;
  elseif strcmp(SET(pos).ImageType,'Vessel_3D_Probability')
    stacks(currentstack).type = 'Vessel_3D_Probability';
    stacks(currentstack).name = sprintf('%s',stacks(currentstack).type);
    stacks(currentstack).nos = pos;
    stacks(currentstack).data = SET(pos);
    stacks(currentstack).timesteps = SET(pos).TIncr*(0:(SET(pos).TSize-1));
    currentstack = currentstack+1;
  
  %QFlow magnitude: SB 24/4
  elseif strcmp(SET(pos).ImageType, 'Flow (magnitude)')%&& (SET(pos).ZSize<2)%strcmp(SET(1).SeriesDescription(1:2),'QF')
%     stacks(currentstack).type = ['QFlow', SET(pos).SeriesDescription(3:end)]; 
stacks(currentstack).type = 'QFlow'; 
    stacks(currentstack).name = sprintf('%s_%d',stacks(currentstack).type,QFlowc); %#ok<AGROW>
    stacks(currentstack).nos = pos; %#ok<AGROW>
    stacks(currentstack).data = SET(pos); %#ok<AGROW>
    %Reset Flow structure to only include magnitude image as an anatomical
    %plane:
    stacks(currentstack).data.Flow = [];
    QFlowc = QFlowc+1;
    %Fix for multi-slice 2D-flow:
    %Fix for 3D_volume RVOT Acquisition: OBS Does not support Segmentation
    %---------------------------------------------------
    if stacks(currentstack).data.ZSize > 1
      %1. Change Imageposition to current slice:
      zdir = cross(SET(pos).ImageOrientation(1:3),SET(pos).ImageOrientation(4:6));
      stacks(currentstack).data.ImagePosition = stacks(currentstack).data.ImagePosition(:)'-(SET(pos).CurrentSlice-1)*(SET(pos).SliceThickness+SET(pos).SliceGap)*zdir(:)';
      
      %2. Set slice parameters to a single slice
      stacks(currentstack).data.ZSize = 1;
      stacks(currentstack).data.StartSlice = 1;
      stacks(currentstack).data.EndSlice = 1;
      stacks(currentstack).data.SliceGap = 0;
      stacks(currentstack).data.OrgZSize = 1;
      
      %3. Take current slice of volume
      stacks(currentstack).data.IM = single(squeeze(SET(pos).IM(:,:,:,SET(pos).CurrentSlice)));
      
      %4. Set final slice parameter correctly
      stacks(currentstack).data.CurrentSlice = 1;
    end
    %---------------------------------------------------
    currentstack = currentstack+1;
  elseif ~isempty(findstr('vortexprob',SET(pos).ImageType))   
    numstr = plugin_exporttofourflow('numfromstr', SET(pos).ImageType); %pos, SET(pos).ImageType);
    
    if length(numstr) == 2
      stacks(currentstack).type = 'vortexprob'; %#ok<*AGROW>
      stacks(currentstack).name = sprintf('%s_%s_%s',stacks(currentstack).type,numstr{1},numstr{2});
      stacks(currentstack).nos = pos;
      stacks(currentstack).data = SET(pos);
      currentstack = currentstack+1;
    else
      error('VortexProb have to have both radius and sigma integers.');
    end
    % To be implemented (ftle, volumetracking)
    %   elseif ~isempty(findstr('ftle',SET(pos).ImageType))
    %     result = 'ftle';
    %   elseif ~isempty(findstr('volumetracking',SET(pos).ImageType))
    %     result = 'volumetracking';
%   elseif 
% 	%Anatomy Volume
% 	stacks(currentstack).type = 'Anatomy Volume';
%     currentstackName = sprintf('%s_%d',stacks(currentstack).type,saxc);
%     stacks(currentstack).name = currentstackName;
%     stacks(currentstack).nos = pos;
%     stacks(currentstack).data = SET(pos);
%     saxc = saxc+1;
%     currentstack = currentstack+1;
  elseif (ismember(pos,cineno) && SET(pos).ZSize>1) || ((SET(pos).TSize == 1) && (SET(pos).ZSize>1))
    stacks(currentstack).type = 'sax';
    currentstackName = sprintf('%s_%d',stacks(currentstack).type,saxc);
    stacks(currentstack).name = currentstackName;
    stacks(currentstack).nos = pos;
    stacks(currentstack).data = SET(pos);
    saxc = saxc+1;
    currentstack = currentstack+1;
    if not(isempty(SET(pos).EndoX))
        stacks(currentstack).type = 'mesh';
        stacks(currentstack).subtype = 'Endo';
        stacks(currentstack).name = sprintf('%s_Endo',currentstackName);
        stacks(currentstack).nos = pos;
        stacks(currentstack).data = SET(pos);
        stacks(currentstack).x = SET(pos).EndoX;
        stacks(currentstack).y = SET(pos).EndoY;
        stacks(currentstack).z = (0:size(stacks(currentstack).x,3)-1);
        stacks(currentstack).timesteps = SET(pos).TIncr*(0:(SET(pos).TSize-1));
        currentstack = currentstack+1;
    end
    if not(isempty(SET(pos).EpiX))
        stacks(currentstack).type = 'mesh';
        stacks(currentstack).subtype = 'Epi';
        stacks(currentstack).name = sprintf('%s_Epi',currentstackName);
        stacks(currentstack).nos = pos;
        stacks(currentstack).data = SET(pos);
        stacks(currentstack).x = SET(pos).EpiX;
        stacks(currentstack).y = SET(pos).EpiY;
        stacks(currentstack).z = (0:size(stacks(currentstack).x,3)-1);
        stacks(currentstack).timesteps = SET(pos).TIncr*(0:(SET(pos).TSize-1));
        currentstack = currentstack+1;
    end
    if not(isempty(SET(pos).RVEndoX))
        stacks(currentstack).type = 'mesh';
        stacks(currentstack).subtype = 'RVEndo';
        stacks(currentstack).name = sprintf('%s_RVEndo',currentstackName);
        stacks(currentstack).nos = pos;
        stacks(currentstack).data = SET(pos);
        stacks(currentstack).x = SET(pos).RVEndoX;
        stacks(currentstack).y = SET(pos).RVEndoY;
        stacks(currentstack).z = (0:size(stacks(currentstack).x,3)-1);
        stacks(currentstack).timesteps = SET(pos).TIncr*(0:(SET(pos).TSize-1));
        currentstack = currentstack+1;
    end
    if not(isempty(SET(pos).RVEpiX))
        stacks(currentstack).type = 'mesh';
        stacks(currentstack).subtype = 'RVEpi';
        stacks(currentstack).name = sprintf('%s_RVEpi',currentstackName);
        stacks(currentstack).nos = pos;
        stacks(currentstack).data = SET(pos);
        stacks(currentstack).x = SET(pos).RVEpiX;
        stacks(currentstack).y = SET(pos).RVEpiY;
        stacks(currentstack).z = (0:size(stacks(currentstack).x,3)-1);
        stacks(currentstack).timesteps = SET(pos).TIncr*(0:(SET(pos).TSize-1));
        currentstack = currentstack+1;
    end
  elseif strcmp(SET(pos).ImageViewPlane,'2CH') && SET(pos).TSize > 1
    stacks(currentstack).type = '2CH';
    stacks(currentstack).name = sprintf('%s_%d',stacks(currentstack).type,CH2c);
    stacks(currentstack).nos = pos;
    stacks(currentstack).data = SET(pos);
    CH2c = CH2c+1;
    currentstack = currentstack+1
  elseif strcmp(SET(pos).ImageViewPlane,'3CH') && SET(pos).TSize > 1
    stacks(currentstack).type = 'lvoutflow';%'3CH';
    stacks(currentstack).name = sprintf('%s_%d',stacks(currentstack).type,CH3c);
    stacks(currentstack).nos = pos;
    stacks(currentstack).data = SET(pos);
    CH3c = CH3c+1;
	currentstackName = stacks(currentstack).name;
    currentstack = currentstack+1;
	
	%Segmentation/Mesh support for standard views in fourflow: 
	if not(isempty(SET(pos).EndoX))
        stacks(currentstack).type = 'mesh';
        stacks(currentstack).subtype = 'Endo';
        stacks(currentstack).name = sprintf('%s_Endo',currentstackName);
        stacks(currentstack).nos = pos;
        stacks(currentstack).data = SET(pos);
        stacks(currentstack).x = SET(pos).EndoX;
        stacks(currentstack).y = SET(pos).EndoY;
        stacks(currentstack).z = (0:size(stacks(currentstack).x,3)-1);
        stacks(currentstack).timesteps = SET(pos).TIncr*(0:(SET(pos).TSize-1));
        currentstack = currentstack+1;
	end
	if not(isempty(SET(pos).EpiX))
		stacks(currentstack).type = 'mesh';
		stacks(currentstack).subtype = 'Epi';
		stacks(currentstack).name = sprintf('%s_Epi',currentstackName);
		stacks(currentstack).nos = pos;
		stacks(currentstack).data = SET(pos);
		stacks(currentstack).x = SET(pos).EpiX;
		stacks(currentstack).y = SET(pos).EpiY;
		stacks(currentstack).z = (0:size(stacks(currentstack).x,3)-1);
		stacks(currentstack).timesteps = SET(pos).TIncr*(0:(SET(pos).TSize-1));
		currentstack = currentstack+1;
	end
	if not(isempty(SET(pos).RVEndoX))
		stacks(currentstack).type = 'mesh';
		stacks(currentstack).subtype = 'RVEndo';
		stacks(currentstack).name = sprintf('%s_RVEndo',currentstackName);
		stacks(currentstack).nos = pos;
		stacks(currentstack).data = SET(pos);
		stacks(currentstack).x = SET(pos).RVEndoX;
		stacks(currentstack).y = SET(pos).RVEndoY;
		stacks(currentstack).z = (0:size(stacks(currentstack).x,3)-1);
		stacks(currentstack).timesteps = SET(pos).TIncr*(0:(SET(pos).TSize-1));
		currentstack = currentstack+1;
	end
	if not(isempty(SET(pos).RVEpiX))
		stacks(currentstack).type = 'mesh';
		stacks(currentstack).subtype = 'RVEpi';
		stacks(currentstack).name = sprintf('%s_RVEpi',currentstackName);
		stacks(currentstack).nos = pos;
		stacks(currentstack).data = SET(pos);
		stacks(currentstack).x = SET(pos).RVEpiX;
		stacks(currentstack).y = SET(pos).RVEpiY;
		stacks(currentstack).z = (0:size(stacks(currentstack).x,3)-1);
		stacks(currentstack).timesteps = SET(pos).TIncr*(0:(SET(pos).TSize-1));
		currentstack = currentstack+1;
	end
  elseif strcmp(SET(pos).ImageViewPlane,'4CH') && SET(pos).TSize > 1
    stacks(currentstack).type = '4CH'; %#ok<AGROW>
    stacks(currentstack).name = sprintf('%s_%d',stacks(currentstack).type,CH4c); %#ok<AGROW>
    stacks(currentstack).nos = pos; %#ok<AGROW>
    stacks(currentstack).data = SET(pos); %#ok<AGROW>
    CH4c = CH4c+1;
    currentstack = currentstack+1
  elseif ~isempty(SET(pos).Flow) && SET(pos).ZSize > 1
    flowstruct = SET(pos).Flow;
    if length([flowstruct.PhaseX flowstruct.PhaseY flowstruct.PhaseNo]) == 3
      stacks(currentstack).type = 'velocity'; %#ok<AGROW>
      stacks(currentstack).name = sprintf('%s_%d',stacks(currentstack).type,velc); %#ok<AGROW>
      stacks(currentstack).nos = pos; %#ok<AGROW>
      stacks(currentstack).data = SET(pos); %#ok<AGROW>
      velc = velc+1;
      pos = pos+3;
      currentstack = currentstack+1
    %else
    %  stacks(currentstack).type = 'notsupported';
    end
  %else
  %  stacks(currentstack).type = 'notsupported';
  end
  pos = pos+1
end

%----------------------
function close_Callback %#ok<DEFNU>
%----------------------
global DATA

disp('Closing.')

%This saves the position
DATA.GUI.FourFlowExport = close(DATA.GUI.FourFlowExport);
 
%------------------------------
function export_Callback
%------------------------------
disp('export_Callback');
global DATA NO

gui = DATA.GUI.FourFlowExport; %get from the global variable
%Ask for filename
temp = pwd;
if exist(DATA.Pref.exportpath,'dir')
  cd(DATA.Pref.exportpath);
else
  mywarning('Export path does not exist, please check preferences.');
end;

pathname = myuigetdir('DATA.Pref.datapath', 'Choose empty directory for export');
choice = 2;
while length(dir(pathname))>2 && choice == 2
  stri = sprintf(['You have chosen a non-empty directory. ' ...
    'Automatically generated files with same name will replace current files if such exists.']);
  choice = mymenu(stri,{'Cancel','Choose new directory','Continue'});
  if choice == 1
    myfailed('Export aborted by user.');
    return
  elseif choice == 2
    pathname = myuigetdir('DATA.Pref.datapath', 'Choose empty directory for export');
  else
    continue
  end
end

%Export all stacks as extended with two periods?
cd(temp);

smallestTime = -1;
largestNumberOfFrames = -1;
stack_export_list = get(gui.handles.stacksListbox,'String');

%Pre-process each selected stacks: 
for n = 1:length(stack_export_list)
    data = gui.stacks(n).data;
    %if smallestTime == -1
    %    smallestTime = data.TIncr;
    %elseif smallestTime < data.TIncr
    %    smallestTime = data.TIncr;    
    %end;
    if largestNumberOfFrames == -1
        largestNumberOfFrames = data.TSize;
        smallestTime = data.TIncr;
    elseif largestNumberOfFrames < data.TSize
        largestNumberOfFrames = data.TSize;    
        smallestTime = data.TIncr;
    end;
end;

if largestNumberOfFrames < 2
	largestNumberOfFrames = 10;
	smallestTime = 0.1;
end

h = waitbar(0,'Please wait exporting data...');
selectedForExportindices = get(gui.handles.stacksListbox,'value');
for n = 1:length(selectedForExportindices)
    currentIndex = selectedForExportindices(n);
    stacktype = gui.stacks(currentIndex).type;
    if(strcmp(stacktype, 'mesh'))
        stack = gui.stacks(currentIndex);
        data = stack.data;
        factor = largestNumberOfFrames/data.TSize;
        data.TIncr = smallestTime;
        data.TSize = largestNumberOfFrames;
        if strcmp(stack.subtype, 'Endo')
            stack.x = tools('upsampletemporal', factor, data.EndoX,'segmentation','nearest');
            stack.y = tools('upsampletemporal', factor, data.EndoY,'segmentation','nearest');
            stack.z = (0:size(stack.x,3)-1);
            stack.timesteps = data.TIncr*(0:(data.TSize-1));
        end
        if strcmp(stack.subtype, 'Epi')
            stack.x = tools('upsampletemporal', factor, data.EpiX,'segmentation','nearest');
            stack.y = tools('upsampletemporal', factor, data.EpiY,'segmentation','nearest');
            stack.z = (0:size(stack.x,3)-1);
            stack.timesteps = data.TIncr*(0:(data.TSize-1));
        end
        if strcmp(stack.subtype, 'RVEndo')
            stack.x = tools('upsampletemporal', factor, data.RVEndoX,'segmentation','nearest');
            stack.y = tools('upsampletemporal', factor, data.RVEndoY,'segmentation','nearest');
            stack.z = (0:size(stack.x,3)-1);
            stack.timesteps = data.TIncr*(0:(data.TSize-1));
		end
		if strcmp(stack.subtype, 'RVEpi')
            stack.x = tools('upsampletemporal', factor, data.RVEpiX,'segmentation','nearest');
            stack.y = tools('upsampletemporal', factor, data.RVEpiY,'segmentation','nearest');
            stack.z = (0:size(stack.x,3)-1);
            stack.timesteps = data.TIncr*(0:(data.TSize-1));
		end  
		tmpno=NO;
		NO=stack.nos;
        matlab2ensightgolddeformingmesh(pathname, stack.name, '', stack.x, stack.y, stack.z, stack.timesteps);
		NO=tmpno;
	else
        data = gui.stacks(currentIndex).data;
        tmpno = NO;
        NO = gui.stacks(currentIndex).nos;
        %factor = smallestTime/data.TIncr;
        factor = largestNumberOfFrames/data.TSize;
        disp('factor selected as: ')
        disp(factor)
		if data.TSize==1
			imdata = squeeze(repmat(data.IM, [1 1 largestNumberOfFrames 1]));
			factor = 1;
		else
			imdata = data.IM;
		end
        data.IM = tools('upsampletemporal', factor, imdata, 'image');
        data.TIncr = smallestTime;
        data.TSize = largestNumberOfFrames;
        [datastruct, timesteps] = preparedataforexport(pathname, stack_export_list{currentIndex}, false, data);
        myworkon;
        matlab2ensightgold_multiple_parts(pathname, stack_export_list{currentIndex}, timesteps, datastruct);
        if(strcmp(stacktype, 'velocity')||strcmp(stacktype, 'Vessel_3D_Probability'))
          matlab2vtk(pathname, stack_export_list{currentIndex}, datastruct);
        end;
        myworkoff;
        NO = tmpno;
    end;
    waitbar(n / length(selectedForExportindices))
end;
close(h) 

function [datastruct, timesteps] = preparedataforexport(~, ~, writetwoperiods, currentData)
global SET
%Set up the data structure
datastruct = [];
datastruct.Data = {};
datastruct.VarName = {};
datastruct.Size = [];
datastruct.Desc = '';
datastruct.Origin = [];
datastruct.Delta = [];
datastruct.DirDim = [];
datastruct.TimeSteps = [];
datastruct.NumberOfTimeSteps = 0;

%--- create datastructure for export
%--- find type of image stack
type = 'anatomical'; %default

if ~isempty(currentData.Flow)
  %Some kind of flow...
  if isempty(currentData.Flow.PhaseNo)&&...
      (~isempty(currentData.Flow.PhaseX))&&...
      (~isempty(currentData.Flow.PhaseY))
    %Through plane missing => in plane
    type = 'flow2dinplane';
  end;
  if isempty(currentData.Flow.PhaseX)&&...
      isempty(currentData.Flow.PhaseY)&&...
      (~isempty(currentData.Flow.PhaseNo))
    type = 'flowthroughplane';
  end;
  if (~isempty(currentData.Flow.PhaseNo))&&...
      (~isempty(currentData.Flow.PhaseX))&&...
      (~isempty(currentData.Flow.PhaseY))
    type = 'flow3d';
  end;
  if isequal(type,'anatomical')
    mywarning('Could not determine type of flow.');
  end;
end; %flow

%--- fix datastruct (general)
datastruct.Size = [currentData.XSize currentData.YSize currentData.ZSize currentData.TSize]; 
datastruct.Origin = currentData.ImagePosition/1000; % %Convert to [m]
datastruct.Delta = [...
  currentData.ResolutionX ...
  currentData.ResolutionY ...
  -(currentData.SliceThickness+currentData.SliceGap)]/1000; % %Convert to [m]
zdir = cross(currentData.ImageOrientation(1:3),currentData.ImageOrientation(4:6));
datastruct.DirDim = [...
  currentData.ImageOrientation(4:6)' currentData.ImageOrientation(1:3)' zdir']; 

%--- image type specific ...
switch type
  case 'anatomical'
    disp('anatomical');
    datastruct.Desc = sprintf('Anatomical %s',currentData.ImageType); 
    datastruct.VarName = sprintf('MagnitudeP%03dN',1); 

    %Ensight uses x y z t
    datastruct.Data = permute(...
      plugin_exporttofourflow('resampletimeframes', currentData.IM,currentData.TSize),...
      [1 2 4 3]);  
  case 'flowthroughplane'
    disp('flowthroughplane');
    %--- Through plane velocity
    datastruct.Desc = 'Through plane vel (mag)'; 
    datastruct.VarName = cell(1,2);
    datastruct.VarName{1} = sprintf('MagnitudeP%03dN',1);
    datastruct.VarName{2} = sprintf('V_through_plane_%02d',1); 
    nop = currentData.Flow.PhaseNo;
    nom = currentData.Flow.MagnitudeNo;
    datastruct.Data = cell(1,2); 
    datastruct.Data{1} = ...
      permute(plugin_exporttofourflow('resampletimeframes', SET(nom).IM ,currentData.TSize),[1 2 4 3]);  
    datastruct.Data{2} = ...
      permute(plugin_exporttofourflow('resampletimeframes', velsegment2real(SET(nop).IM, SET(nop).VENC), currentData.TSize),[1 2 4 3]);
  case 'flow2dinplane'
    disp('flow2dinplane');
    %--- In plane velocity encoded
    datastruct.Desc = sprintf('In-plane vel %s',currentData.ImageType); 
    datastruct.VarName = sprintf('V_in_plane_%02d',1); 
    nox = currentData.Flow.PhaseX;
    noy = currentData.Flow.PhaseY;
    nom = currentData.Flow.MagnitudeNo;
    datastruct.Data = cell(1,2); 
    datastruct.Data{1} = ...
      permute(plugin_exporttofourflow('resampletimeframes', velsegment2real(SET(nom).IM,SET(nom).VENC),currentData.TSize),[1 2 4 3]); 

    %Set third component to zero.
    datastruct.Data{2} = cat(5,...
      permute(plugin_exporttofourflow('resampletimeframes', velsegment2real(SET(noy).IM,SET(noy).VENC),currentData.TSize),[1 2 4 3]),...
      permute(plugin_exporttofourflow('resampletimeframes', velsegment2real(SET(nox).IM,SET(nox).VENC),currentData.TSize),[1 2 4 3]),...
      zeros(size(datastruct.Data{1})));        

  case 'flow3d'
    disp('flow3d');
    %--- 3D flow
    datastruct.Desc = sprintf('3Dflow%02d',1); 
    datastruct.VarName = sprintf('V_3D_%02d',1); 
    nox = currentData.Flow.PhaseX;
    noy = currentData.Flow.PhaseY;
    noz = currentData.Flow.PhaseNo;

    datastruct.Data = cat(5,...
      permute(plugin_exporttofourflow('resampletimeframes', velsegment2real(SET(nox).IM,SET(nox).VENC),currentData.TSize),[1 2 4 3]),...
      permute(plugin_exporttofourflow('resampletimeframes', velsegment2real(SET(noy).IM,SET(noy).VENC),currentData.TSize),[1 2 4 3]),...
      permute(plugin_exporttofourflow('resampletimeframes', velsegment2real(SET(noz).IM,SET(noz).VENC),currentData.TSize),[1 2 4 3]));           
  otherwise
    myfailed(dprintf('Type %s not recognized.',type));
    return;
end;
currentData.TIncr;

% Generate timestep vector
timesteps = currentData.TIncr*(0:(currentData.TSize-1));
datastruct.TimeSteps = timesteps;
datastruct.NumberOfTimeSteps = length(timesteps);

% Extend to two periods if user wants to
if writetwoperiods
  % Duplicated timestep vector
  timesteps = currentData.TIncr*(0:(2*(currentData.TSize)-1));

    % Change time size
    datastruct.Size(4) = 2*datastruct.Size(4);
    
    % Handle single Data and cell Data
    if isa(datastruct.Data,'cell')
      % JT: This branch is untested.
      for jj = 1:length(datastruct.Data)
        datastruct.Data{jj} = cat(4, datastruct.Data{jj}, datastruct.Data{jj});
      end
    else
      datastruct.Data = cat(4, datastruct.Data, datastruct.Data);
    end 
    
    length(timesteps)
end;

%----------------------------------------------------
function [vreal] = velsegment2real(vsegment, venc)
%----------------------------------------------------
% Helper function.
% Converts segment velocities to 'real' velocities.
%
% real           segment
%  -venc           0
%  0               1/2
%  +venc           1
%
% venc is assumed to be in cm/s, and v_real is in m/s.

vreal = (vsegment - 1/2)*2*venc/100;