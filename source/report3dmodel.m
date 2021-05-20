function report3dmodel(varargin)
% REPORT3DMODEL
% Functions for doing 3D model reports

% Moved out from segment_main by Nisse Lundahl

%GUI to show 3D display of segmentation.

% TODO
% - get longaxismotion corrected
% - rvendo may be improperly sampled radially. This should be fixed in
%   checkconsistancy.
% - -clean up persistants...
% Check if possible not to use grid

global DATA SET 

if nargin > 0 && isequal(varargin{1},'rotated')
  report3dmodelrotated;
  return
end

if nargin==0
  %------------------------------------------ No input arguments => setup
  init;
else
  gui = DATA.GUI.Report3DModel;
  %Input arguments start different callbacks
  switch varargin{1}
    case 'play'
      DATA.StartFrame = gui.t;
      DATA.StartTime = now;
      while ~gui.closing && get(gui.handles.playtogglebutton,'value')
        gui.t = segment('getframenumber');
        pause(0.5*SET(gui.no).BeatTime/SET(gui.no).TSize);
        report3dmodel('update');
      end
    case 'radialvel'
      if get(gui.handles.radialvelradiobutton,'value')
        gui.handles.radvel = calcfunctions('calcradialvelocity',gui.no);
      else
        set(DATA.Handles.endosurf,...
          'Cdata',repmat(97,size(EndoX,1),size(EndoX,3)));
      end
      report3dmodel('update');
    case 'prev'
      i=find(gui.t==gui.indTn);
      if (i==1)
        gui.t=gui.indTn(end);
      else
        gui.t=gui.indTn(i-1);
      end
      report3dmodel('update');
    case 'next'
      i=find(gui.t==gui.indTn);
      if (i==length(gui.indTn))
        gui.t=gui.indTn(1);
      else
        gui.t=gui.indTn(i+1);
      end
      report3dmodel('update');
    case 'dias'
      gui.t=SET(gui.no).EDT;
      report3dmodel('update');
    case 'sys'
      gui.t=SET(gui.no).EST;
      report3dmodel('update');
    case 'close'
      gui.closing=true;
      pause(0.5)
      %pause(0.5*SET(gui.no).BeatTime/SET(gui.no).TSize);
      %Close the windows
      try
        delete(4);
      catch %#ok<CTCH>
        %Do nothing if failed
      end
      try
        delete(gui.fig);
      catch %#ok<CTCH>
        %Do nothing if failed
      end
      DATA.GUI.Report3DModel = [];
    case 'update'
      %--- Refresh the display
      if (gui.closing)
        try
          set(gui.handles.playtogglebutton,'value',0);
        catch %#ok<CTCH>
        end
        return;
      end
      %try
      %Update endosurf
      if get(gui.handles.endocheckbox,'value')
        set(DATA.Handles.endosurf,...
          'visible','on',...
          'Xdata',squeeze(gui.EndoX(:,gui.t,:)),... %EH:
          'Ydata',squeeze(gui.EndoY(:,gui.t,:)),... %EH:
          'Zdata',gui.EndoZ,...
          'Alphadata',mygetvalue(gui.handles.endoslider)*ones(size(gui.EndoZ)));
        if get(gui.handles.radialvelradiobutton,'value')
          %97 = 64+1+32
          set(DATA.Handles.endosurf,...
            'Cdata',97+(32/DATA.MaxRadialVel)*squeeze(gui.handles.radvel(:,t,ind)));
        end
      else
        try
          set(DATA.Handles.endosurf,'visible','off');
        catch %#ok<CTCH>
        end
        
      end

      %Update episurf
      if get(gui.handles.epicheckbox,'value')
        set(DATA.Handles.episurf,...
          'visible','on',...
          'Xdata',squeeze(gui.EpiX(:,gui.t,:)),... %EH:
          'Ydata',squeeze(gui.EpiY(:,gui.t,:)),... %EH:
          'Zdata',gui.EpiZ,...
          'Alphadata',mygetvalue(gui.handles.epislider)*ones(size(gui.EpiZ)));
      else
        try
          set(DATA.Handles.episurf,'visible','off');
        catch %#ok<CTCH>
        end
      end

      %Update RVendo surf
      if get(gui.handles.rvendocheckbox,'value')
        set(DATA.Handles.rvendosurf,...
          'visible','on',...
          'Xdata',squeeze(gui.RVEndoX(:,gui.t,:)),... %EH:
          'Ydata',squeeze(gui.RVEndoY(:,gui.t,:)),... %EH:
          'Zdata',gui.RVEndoZ,...
          'Alphadata',mygetvalue(gui.handles.rvendoslider)*ones(size(gui.RVEndoZ)));
      else
        try
          set(DATA.Handles.rvendosurf,'visible','off');
        catch %#ok<CTCH>
        end
      end
      
      %Update RVepi surf
      if get(gui.handles.rvepicheckbox,'value')
        set(DATA.Handles.rvepisurf,...
          'visible','on',...
          'Xdata',squeeze(gui.RVEpiX(:,gui.t,:)),... %EH:
          'Ydata',squeeze(gui.RVEpiY(:,gui.t,:)),... %EH:
          'Zdata',gui.RVEpiZ,...
          'Alphadata',mygetvalue(gui.handles.rvepislider)*ones(size(gui.RVEpiZ)));
      else
        try
          set(DATA.Handles.rvepisurf,'visible','off');
        catch %#ok<CTCH>
        end
      end
      
      %Update annotation points
      if mygetvalue(gui.handles.annotationpointscheckbox)
        set(DATA.Handles.annotationpoints3d, ...
          'visible','on', ...
          'XData',gui.PointX, ...
          'YData',gui.PointY, ...
          'ZData',gui.PointZ);
      else
        try
          set(DATA.Handles.annotationpoints3d,'visible','off');
        catch %#ok<CTCH>
        end
      end
      
      
      %Update measurements
      if mygetvalue(gui.handles.measurementscheckbox)
        for mloop = 1:numel(gui.MeasureX)
        set(DATA.Handles.measurements3d(mloop), ...
          'visible','on', ...
          'XData',gui.MeasureX{mloop}, ...
          'YData',gui.MeasureY{mloop}, ...
          'ZData',gui.MeasureZ{mloop});
        end
      else
        try
          set(DATA.Handles.measurements3d,'visible','off');
        catch %#ok<CTCH>
        end
      end
      
      %Update scar
      if mygetvalue(gui.handles.scarcheckbox)
        set(DATA.Handles.scar3d, ...
          'visible','on', ...
          'XData',gui.ScarX, ...
          'YData',gui.ScarY, ...
          'ZData',gui.ScarZ);
      else
        try
          set(DATA.Handles.scar3d,'visible','off');
        catch %#ok<CTCH>
        end
      end

      %Update short axis slice
      if get(gui.handles.shortslicescheckbox,'value')
        slice = 1+SET(gui.no).ZSize-round(mygetvalue(gui.handles.shortsliceslider));
        %the most basal should give z=0
        z = repmat(...
          (SET(gui.no).SliceThickness+SET(gui.no).SliceGap)*(slice-gui.basalslice),...
          SET(gui.no).XSize,SET(gui.no).YSize);
        set(DATA.Handles.lvcutshort,'Zdata',-z);
        set(DATA.Handles.lvcutshort,'visible','on',...
          'CDATA',1+(63*double(squeeze(SET(gui.no).IM(:,:,gui.t,slice)))));
      else
        set(DATA.Handles.lvcutshort,'visible','off');
      end

      %Update long-axis slice
      if get(gui.handles.longslicescheckbox,'value')
        z = max(min(SET(gui.no).YSize*mygetvalue(gui.handles.longsliceslider),SET(gui.no).YSize),1);
        set(DATA.Handles.lvcutlong,...
          'visible','on',...
          'Ydata',repmat(SET(gui.no).ResolutionY*z,[SET(gui.no).XSize SET(gui.no).ZSize]),...
          'Cdata',1+(63*double(squeeze(SET(gui.no).IM(:,round(z),gui.t,:)))));
      else
        set(DATA.Handles.lvcutlong,'visible','off');
      end

      if DATA.Record
        drawnow;
        DATA.MovieFrame = mygetframe(4);
        export('exportmovierecorder_Callback','newframe');
      end

      %catch
      %Did not manage to update correctly, give it up!
      %  report3dmodel('close');
      %end
  end
end


%----------------------------
function report3dmodelrotated
%----------------------------
%Display 3D model of rotated image stacks.
global SET NO

fig = figure(88);
clf;
hold on;
stateandicon = viewfunctions('iconson',{'hidelv','hiderv'});
state = [stateandicon{:,1}];
if state(1)
  if ~isempty(SET(NO).EndoX)
    rotatedplothelper(...
      SET(NO).EndoX,...
      SET(NO).EndoY,...
      'r-');
  end
  
  if ~isempty(SET(NO).EpiX)
    rotatedplothelper(...
      SET(NO).EpiX,...
      SET(NO).EpiY,...
      'g-');
  end
end

if state(2)  
  if ~isempty(SET(NO).RVEndoX)
    rotatedplothelper(...
      SET(NO).RVEndoX,...
      SET(NO).RVEndoY,...
      'm-');
  end

  if ~isempty(SET(NO).RVEpiX)
    rotatedplothelper(...
      SET(NO).RVEpiX,...
      SET(NO).RVEpiY,...
      'c-');
  end
end

ind = find(SET(NO).Point.T==SET(NO).CurrentTimeFrame);
for loop=1:length(ind)
  [x,y,z] = calcfunctions('cyl2cart',...
    SET(NO).Point.X(ind(loop)),...
    SET(NO).Point.Y(ind(loop)),...
    SET(NO).Point.Z(ind(loop)),...
    SET(NO).ZSize,...
    SET(NO).RotationCenter);
  plot3(x,y,z,'r*');
end

axis image off;
hold off;

set(fig,'Numbertitle','off','Name','3D visualization of rotated stacks.');
cameratoolbar(fig,'show');

%------------------------------------------
function rotatedplothelper(xc,yc,linecolor)
%------------------------------------------
%Helper function to display segmentation in rotated image stacks.
global SET NO

tf = SET(NO).CurrentTimeFrame;

xmin = 1e10;
xmax = -1e10;
for zloop=1:SET(NO).ZSize
  xmin = min(min(xc(:,tf,zloop)),xmin);
  xmax = max(max(xc(:,tf,zloop)),xmax);
end

xrange = linspace(xmin,xmax,20);

for zloop=1:SET(NO).ZSize
  %drawcontour
  [x,y,z] = calcfunctions('cyl2cart',...
    xc(:,tf,zloop),...
    yc(:,tf,zloop),...
    zloop,...
    SET(NO).ZSize,... %Number of slices
    SET(NO).RotationCenter); %Rotation center
  plot3(x,y,z,linecolor);   
end

%---Draw_lines

for loop=1:length(xrange)
  %Loop over rows

  r = [];
  slice = [];
  for zloop=1:SET(NO).ZSize
    
    xt = xc(:,tf,zloop);
    yt = yc(:,tf,zloop);
    
    %Find indices that crosses
    ind = find(...
      ((xt(2:end)>xrange(loop))&(xt(1:(end-1))<xrange(loop)))|...
      ((xt(2:end)<xrange(loop))&(xt(1:(end-1))>xrange(loop))));
    for iloop=1:length(ind)
      if yt(ind(iloop))>SET(NO).RotationCenter
        r = [r yt(ind(iloop))-SET(NO).RotationCenter]; %#ok<AGROW>
        slice = [slice zloop]; %#ok<AGROW>
      else
        r = [r SET(NO).RotationCenter-yt(ind(iloop))]; %#ok<AGROW>
        slice = [slice zloop+SET(NO).ZSize]; %#ok<AGROW>
      end
    end
  end
  
  %--- "Sort it"
  while ~isempty(slice)
    drawr = [];
    drawslice = [];

    nextslice = min(slice); %Start with first slice
    while any(slice==nextslice)
      %Find maximum radius at nextslice
      [~,ind] = max(r.*(slice==nextslice));

      %Add to list of lines to draw
      drawr = [drawr r(ind)]; %#ok<AGROW>
      drawslice = [drawslice slice(ind)]; %#ok<AGROW>

      %Remove taken point from list
      logim = true(size(r));
      logim(ind) = false;
      r = r(logim);
      slice = slice(logim);

      %Continue along the slices
      nextslice = nextslice+1;
    end

    %--- Plot the line section
    if length(drawr)>1
      %Upsample it
      drawr = interp1(drawr,linspace(1,length(drawr),length(drawr)*10),'pchip');
      drawslice = interp1(drawslice,linspace(1,length(drawslice),length(drawslice)*10),'linear');

      %Make it go around
      if 1==1
        if ...
            isequal(drawslice(1),1)&&...
            isequal(drawslice(end),SET(NO).ZSize*2)
          drawslice = [drawslice linspace(0,1,10)]; %#ok<AGROW>
          drawr = [drawr interp1([drawr(end) drawr(1)],linspace(1,2,10))]; %#ok<AGROW>
        end
      end
      
      angle = (drawslice-1)*pi/SET(NO).ZSize;
      x = drawr.*sin(angle)+SET(NO).RotationCenter;
      y = drawr.*cos(angle)+SET(NO).RotationCenter;
      z = repmat(xrange(loop),size(x));
      
      plot3(x,y,z,linecolor);
    end
  end
end %Loop over rows

%------------
function init
%------------
%initiating the 3D model interface
global DATA SET NO

if SET(NO).Rotated
  report3dmodelrotated;
  return;
end

if SET(NO).ZSize<2
  myfailed('Need to have more than one slice for viewing 3D-model.',DATA.GUI.Segment);
  return;
end

if isempty(SET(NO).EndoX) || all(isnan(SET(NO).EndoX(:)))
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end

tempnos=NO;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

%Initialize
  if SET(NO).TSize>1
    if isequal(SET(NO).EDT,SET(NO).EST)
      mywarning('Warning, end-diastole occurs at the same time as end-systole. Use autodetect under edit menu.',DATA.GUI.Segment);
    end
  end

%   segment('checkconsistency',1:SET(NO).TSize,1:SET(NO).ZSize);

  doepi = ~isempty(SET(NO).EpiX) && ~all(isnan(SET(NO).EpiX(:)));
  dorvendo = ~isempty(SET(NO).RVEndoX) && ~all(isnan(SET(NO).RVEndoX(:)));
  dorvepi = ~isempty(SET(NO).RVEpiX) && ~all(isnan(SET(NO).RVEpiX(:)));
  dopoints = ~isempty(SET(NO).Point.X);
  domeasures = ~isempty(SET(NO).Measure);
  doscar = ~isempty(SET(NO).Scar);

  %Check data
  ind = any(~isnan(SET(NO).EndoX(1,:,:)),2);
  if doepi
    ind = ind|any(~isnan(SET(NO).EpiX(1,:,:)),2);
  end
  if dorvendo
    ind = ind|any(~isnan(SET(NO).RVEndoX(1,:,:)),2);
  end
  if dorvepi
    ind = ind|any(~isnan(SET(NO).RVEpiX(1,:,:)),2);
  end
  nslices = sum(ind);
  
  if (size(SET(NO).EndoX,2)>1)
    indT=(sum(~isnan(SET(NO).EndoX(1,:,:)),3)>0);
  else
    indT=1;
  end
  indTn=find(indT);
  t=indTn(1);

  fendo=find(any(~isnan(SET(NO).EndoX(1,:,:)),2));
  if doepi
    fepi = find(any(~isnan(SET(NO).EpiX(1,:,:)),2));
  end
  if dorvendo
    frvendo = find(any(~isnan(SET(NO).RVEndoX(1,:,:)),2));
  end
  if dorvepi
    frvepi = find(any(~isnan(SET(NO).RVEpiX(1,:,:)),2));
  end

  if ((1+fendo(end)-fendo(1))~=length(fendo)) ||...
     (doepi&&((1+fepi(end)-fepi(1))~=length(fepi))) ||...
     (dorvendo &&((1+frvendo(end)-frvendo(1))~=length(frvendo))) || ...
     (dorvepi &&((1+frvepi(end)-frvepi(1))~=length(frvepi)))
      myfailed('There are gaps in the segmentation, continuous segmentation required.',DATA.GUI.Segment);
    return
  end
  
  %longaxis = SET(NO).Longaxis-1; %=>mm
    
  %Create plot figure
  plotfig = figure(4);
  setupicon(plotfig);
  set(plotfig,'units','pixels','position',[193 200 500 500]);
  set(plotfig,'renderer','opengl',...
    'Name','3D model',...
    'numbertitle','off',...
    'menubar','none',...
    'CloseRequestFcn',sprintf('report3dmodel(''close'')'));
  %This CloseRequestFcn is set within guide for plotmodel.fig
  
  %Open GUI
  gui = mygui('plotmodel.fig');
  DATA.GUI.Report3DModel = gui;
    
  % Use system color scheme for figure:
  set(gui.fig,'Color',DATA.GUISettings.BackgroundColor); %get(0,'defaultUicontrolBackgroundColor'));

  gui.closing=false;
    
  %Which one is the most basal slice
%   if nslices>1
%     gui.basalslice = find(ind);
%     gui.basalslice = gui.basalslice(1);
%     gui.numslice = nslices;
%   else
%     gui.basalslice = 1;
%     gui.numslice = 0;
%   end
  gui.numslices=nslices;
  gui.basalslice=1;
  gui.t = t;
  gui.indTn = indTn;
  gui.no = NO;
  
  %Modify GUI components
  temp = get(gui.fig,'units');
  set(gui.fig,'units','pixels');
  set(gui.fig,'position',[700 200 300 500]);
  set(gui.fig,'units',temp);

  set(gui.handles.shortsliceslider,...
    'min',1,...
    'max',SET(NO).ZSize,...
    'value',1, ...
    'SliderStep',1/(SET(NO).ZSize-1)*[1 1]);
  set(gui.handles.longsliceslider,...
    'min',1/(SET(NO).XSize-1),...
    'max',1,...
    'value',1/2, ...
    'SliderStep',1/(SET(NO).XSize-1)*[1 1]);
  set(gui.handles.shortslicescheckbox,'value',1);
  set(gui.handles.longslicescheckbox,'value',1);
  set(gui.handles.endocheckbox,'value',1);
  
  load Icons.mat
  set(gui.handles.previcon,'CData',icon.prev);
  set(gui.handles.playtogglebutton,'CData',icon.play);
  set(gui.handles.nexticon,'CData',icon.next);

  set(gui.handles.playtogglebutton,'value',0);
  set(gui.handles.nexticon,'value',0);
  set(gui.handles.previcon,'value',0);
  if ~all(indT)
    set(gui.handles.playtogglebutton,'enable','off');  
  end
  if isempty(find(gui.indTn==SET(NO).EST, 1))
    set(gui.handles.systolebutton,'enable','off');
  end    
  if isempty(find(gui.indTn==SET(NO).EDT, 1))
    set(gui.handles.diastolebutton,'enable','off');
  end    
  if (size(gui.indTn)==1)
    set(gui.handles.playtogglebutton,'enable','off');
    set(gui.handles.nexticon,'enable','off');
    set(gui.handles.previcon,'enable','off');
    set(gui.handles.diastolebutton,'enable','off');
    set(gui.handles.systolebutton,'enable','off');
  end
  
  if ~doepi
    set(gui.handles.epicheckbox,'enable','off');
  else
    set(gui.handles.epicheckbox,'value',1);
  end
  if ~dorvendo
    set(gui.handles.rvendocheckbox,'enable','off');
  else
    set(gui.handles.rvendocheckbox,'value',1);
  end
  if ~dorvepi
    set(gui.handles.rvepicheckbox,'enable','off');
  else
    set(gui.handles.rvepicheckbox,'value',1);
  end
  if ~dopoints
    set(gui.handles.annotationpointscheckbox,'enable','off');
  else
    set(gui.handles.annotationpointscheckbox,'value',1);
  end
  if ~domeasures
    set(gui.handles.measurementscheckbox,'enable','off');
  else
    set(gui.handles.measurementscheckbox,'value',1);
  end
  if ~doscar
    set(gui.handles.scarcheckbox,'enable','off');
  else
    set(gui.handles.scarcheckbox,'value',1);
  end

  if nslices<=2
    %No slices => disable them
    set(gui.handles.endocheckbox,'value',0);
    set(gui.handles.epicheckbox,'value',0);
    set(gui.handles.rvendocheckbox,'value',0);
    set(gui.handles.rvepicheckbox,'value',0);
    set([...
      gui.handles.endocheckbox ...
      gui.handles.epicheckbox...
      gui.handles.rvendocheckbox...
      gui.handles.rvepicheckbox],'enable','off');
  end

  if nslices>2
    gui.EndoX = SET(NO).ResolutionX*SET(NO).EndoX;
    gui.EndoY = SET(NO).ResolutionY*SET(NO).EndoY;
    gui.EndoZ = repmat((1-(1:SET(NO).ZSize))*(SET(NO).SliceThickness+SET(NO).SliceGap),...
      [DATA.NumPoints 1]);
    
    if doepi
      gui.EpiX = SET(NO).ResolutionX*SET(NO).EpiX;
      gui.EpiY = SET(NO).ResolutionY*SET(NO).EpiY;
      gui.EpiZ = repmat((1-(1:SET(NO).ZSize))*(SET(NO).SliceThickness+SET(NO).SliceGap),...
        [DATA.NumPoints 1]);
    end
    
    
%--- Code to upsample contours in Z. Prettier when working, but causes
%many problems, notably when aligning LV with RV. To work, one needs to
%interpolate RV as well, and keep track of the zslices that are releveant
%for each. /JU

%     %Prepare for endocardium
%     tempxendo = SET(NO).ResolutionX*SET(NO).EndoX(:,:,ind);
%     tempyendo = SET(NO).ResolutionY*SET(NO).EndoY(:,:,ind);
% 
%     if doepi
%       tempxepi  = SET(NO).ResolutionX*SET(NO).EpiX(:,:,ind);
%       tempyepi  = SET(NO).ResolutionY*SET(NO).EpiY(:,:,ind);
%     end    
%     if dorv
%       tempxrvendo  = SET(NO).ResolutionX*SET(NO).RVEndoX(:,:,ind);
%       tempyrvendo  = SET(NO).ResolutionY*SET(NO).RVEndoY(:,:,ind);
%     end
% 
%     newzsize = max(16,nslices); %Number of planes to resample it to.
%     EndoX = zeros(...
%       size(SET(NO).EndoX,1),...
%       size(SET(NO).EndoX,2),...
%       newzsize);
%     EndoY = EndoX;
%     EndoZ = zeros(SET(NO).TSize,newzsize); %sum(ind) is number of planes segmented.
%     if doepi
%       EpiX = zeros(...
%         size(SET(NO).EpiX,1),...
%         size(SET(NO).EpiX,2),...
%         newzsize);
%      EpiY = EpiX;
%      EpiZ = EndoZ;
%     end
%     if dorv
%       RVEndoX = zeros(...
%         size(SET(NO).RVEndoX,1),...
%         size(SET(NO).RVEndoX,2),...
%         newzsize);
%      RVEndoY = EndoX;
%      RVEndoZ = EndoZ;
%     end    
%     for tloop=gui.indTn
%       %Find correct factor
% 
%       %temp2 below was NaN'ing when not all timeframes have been segmented
%       %ESV/EDV are not defined unless segmentation has Sys/Dias time
%       %calculations. Going [0..1] works.
%       if (SET(NO).TSize==1) || ...
%          any(isnan([SET(NO).EDV SET(NO).ESV SET(NO).LVV(tloop)]))
%         newind = linspace(0,1,newzsize);
%       else
%         temp = (longaxis/(SET(NO).SliceThickness+SET(NO).SliceGap))/nslices;
%         %temp=>0..1 normalized over nslices number of slices to remove
%         temp2 = (1-(SET(NO).LVV(tloop)-SET(NO).ESV)/(SET(NO).EDV-SET(NO).ESV));
%         %temp2=>0 at EDT, temp2=>1 at EST
%         temp = temp*temp2;
%         newind = linspace(temp,1,newzsize);
%       end
%       %Resample endocardium
%       EndoX(:,tloop,:) = interp1(...
%         linspace(0,1,nslices)',...
%         squeeze(tempxendo(:,tloop,:))',...
%         newind')';
%       EndoY(:,tloop,:) = interp1(...
%         linspace(0,1,nslices)',...
%         squeeze(tempyendo(:,tloop,:))',...
%         newind')';
%       if doepi
%         %Resample epicardium
%         EpiX(:,tloop,:) = interp1(...
%           linspace(0,1,nslices)',...
%           squeeze(tempxepi(:,tloop,:))',...
%           newind')';
%         EpiY(:,tloop,:) = interp1(...
%           linspace(0,1,nslices)',...
%           squeeze(tempyepi(:,tloop,:))',...
%           newind')';
%       end
%       if dorv
%         %Resample RV endocardium
%         RVEndoX(:,tloop,:) = interp1(...
%           linspace(0,1,nslices)',...
%           squeeze(tempxrvendo(:,tloop,:))',...
%           newind')';
%         RVEndoY(:,tloop,:) = interp1(...
%           linspace(0,1,nslices)',...
%           squeeze(tempyrvendo(:,tloop,:))',...
%           newind')';
%       end
%       EndoZ(tloop,:) = newind*...
%         sum(ind)*(SET(NO).SliceThickness+SET(NO).SliceGap);
%       EpiZ(tloop,:) = newind*...
%         sum(ind)*(SET(NO).SliceThickness+SET(NO).SliceGap);
%       RVEndoZ(tloop,:) = newind*...
%         sum(ind)*(SET(NO).SliceThickness+SET(NO).SliceGap);
%     end
%     clear tempxendo tempyendo tempxepi tempyepi;

  %--------------
  %--- Prepare RV
  %--------------
  if dorvendo
    %RVEndoX RVEndoY RVEndoZ
    
    %--- Create RV EndoZ, EpiZ
    %gui.RVEndoZ = repmat(-(1:SET(NO).ZSize)*(SET(NO).SliceThickness+SET(NO).SliceGap),...
    %  [DATA.NumPoints 1]); %Changed 2020-07-26 after bug-report

    gui.RVEndoZ = repmat((1-(1:SET(NO).ZSize))*(SET(NO).SliceThickness+SET(NO).SliceGap),...
      [DATA.NumPoints 1]);

    %--- Create RV EndoX, EndoY, EpiX, EpiY
    gui.RVEndoX = SET(NO).ResolutionX*SET(NO).RVEndoX;
    gui.RVEndoY = SET(NO).ResolutionY*SET(NO).RVEndoY;

    if 1==1
    %--- Fix with contours to re-allign them with a distant point.
    for tloop=1:SET(NO).TSize
      %Find center.
      tempx = gui.RVEndoX(:,tloop,:);
      tempx = tempx(:);
      tempx = tempx(~isnan(tempx));
      meanx = mean(tempx);

      tempy = gui.RVEndoY(:,tloop,:);
      tempy = tempy(:);      
      tempy = tempy(~isnan(tempy));      
      meany = mean(tempy);
      
      %find point furthest away.
      dist = sqrt((tempx-meanx).^2+(tempy-meany).^2);
      [~,ind]=max(dist);
      if isempty(tempx)
        px = NaN;
        py = NaN;
      else
        px = tempx(ind);
        py = tempy(ind);
      end
      
      for zloop=1:SET(NO).ZSize
        xr = gui.RVEndoX(:,tloop,zloop);
        yr = gui.RVEndoY(:,tloop,zloop);        
        
        %Find point closest to cursor
        [~,inda] = min(abs(complex(px-xr,py-yr)));
        
        gui.RVEndoX(1:(DATA.NumPoints-inda+1),tloop,zloop) = xr(inda:end);
        gui.RVEndoY(1:(DATA.NumPoints-inda+1),tloop,zloop) = yr(inda:end);
        gui.RVEndoX((DATA.NumPoints+1-inda):end,tloop,zloop) = xr(1:inda);
        gui.RVEndoY((DATA.NumPoints+1-inda):end,tloop,zloop) = yr(1:inda);
      end
    end
    end 
  end
  
  if dorvepi
    %RVEndoX RVEndoY RVEndoZ
    
    %--- Create RVEndoZ
    %gui.RVEpiZ = repmat(-(1:SET(NO).ZSize)*(SET(NO).SliceThickness+SET(NO).SliceGap),...
    %  [DATA.NumPoints 1]); %Fixed 2020-07-26 after bug report.
    gui.RVEpiZ = repmat((1-(1:SET(NO).ZSize))*(SET(NO).SliceThickness+SET(NO).SliceGap),...
      [DATA.NumPoints 1]);


    %--- Create EndoX, EndoY
    gui.RVEpiX = SET(NO).ResolutionX*SET(NO).RVEpiX;
    gui.RVEpiY = SET(NO).ResolutionY*SET(NO).RVEpiY;

    if 1==1
    %--- Fix with contours to re-allign them with a distant point.
    for tloop=1:SET(NO).TSize
      %Find center.
      tempx = gui.RVEpiX(:,tloop,:);
      tempx = tempx(:);
      tempx = tempx(~isnan(tempx));
      meanx = mean(tempx);

      tempy = gui.RVEpiY(:,tloop,:);
      tempy = tempy(:);      
      tempy = tempy(~isnan(tempy));      
      meany = mean(tempy);
      
      %find point furthest away.
      dist = sqrt((tempx-meanx).^2+(tempy-meany).^2);
      [~,ind]=max(dist);
      if isempty(tempx)
        px = NaN;
        py = NaN;
      else
        px = tempx(ind);
        py = tempy(ind);
      end
      
      for zloop=1:SET(NO).ZSize
        xr = gui.RVEpiX(:,tloop,zloop);
        yr = gui.RVEpiY(:,tloop,zloop);        
        
        %Find point closest to cursor
        [~,inda] = min(abs(complex(px-xr,py-yr)));
        
        gui.RVEpiX(1:(DATA.NumPoints-inda+1),tloop,zloop) = xr(inda:end);
        gui.RVEpiY(1:(DATA.NumPoints-inda+1),tloop,zloop) = yr(inda:end);
        gui.RVEpiX((DATA.NumPoints+1-inda):end,tloop,zloop) = xr(1:inda);
        gui.RVEpiY((DATA.NumPoints+1-inda):end,tloop,zloop) = yr(1:inda);
      end
    end
    end 
  end
  
  %Prepare annotation points
  if dopoints
    gui.PointX = SET(NO).ResolutionX*SET(NO).Point.X;
    gui.PointY = SET(NO).ResolutionY*SET(NO).Point.Y;
    gui.PointZ = (SET(NO).SliceThickness + SET(NO).SliceGap) * ...
      (-SET(NO).Point.Z + 1);
  end
  
  %Prepare measurements
  if domeasures
    gui.MeasureX = cellfun(@(x)(SET(NO).ResolutionX*x), ...
      {SET(NO).Measure.X},'UniformOutput',false);
    gui.MeasureY = cellfun(@(x)(SET(NO).ResolutionY*x), ...
      {SET(NO).Measure.Y},'UniformOutput',false);
    gui.MeasureZ = cellfun(@(x)((SET(NO).SliceThickness + SET(NO).SliceGap) * ...
      (-x + 1)),{SET(NO).Measure.Z},'UniformOutput',false);
  end
  
  %Prepare scar data
  if doscar
    ScarX = [];
    ScarY = [];
    ScarZ = [];
    for z = 1:SET(NO).ZSize
      scarcontour = contourc(double(SET(NO).Scar.Result(:,:,z)),1);
      defcol = 1;
      while defcol < size(scarcontour,2)
        numpoints = scarcontour(2,defcol);
        ScarX = [ScarX scarcontour(2,defcol+(1:numpoints))]; %#ok<AGROW>
        ScarY = [ScarY scarcontour(1,defcol+(1:numpoints))]; %#ok<AGROW>
        ScarZ = [ScarZ z*ones(1,numpoints)]; %#ok<AGROW>
        defcol = defcol + numpoints + 1;
      end
    end
    gui.ScarX = SET(NO).ResolutionX*ScarX;
    gui.ScarY = SET(NO).ResolutionY*ScarY;
    gui.ScarZ = (SET(NO).SliceThickness + SET(NO).SliceGap) * (-ScarZ+1);
  end
  
%-----  Check if need to add apex, currently disabled.
    if nslices==-1    %was>3
      gui.EndoX(:,:,end+1) = mean(mean(gui.EndoX(:,:,1)));
      gui.EndoY(:,:,end+1) = mean(mean(gui.EndoY(:,:,1)));
      gui.EndoZ(:,end+1) = gui.EndoZ(:,end)+(SET(NO).SliceThickness+SET(NO).SliceGap);
      if doepi
        gui.EpiX(:,:,end+1) = mean(mean(gui.EpiX(:,:,1)));
        gui.EpiY(:,:,end+1) = mean(mean(gui.EpiY(:,:,1)));
        gui.EpiZ(:,end+1) = gui.EpiZ(:,end)+(SET(NO).SliceThickness+SET(NO).SliceGap);
      end
      if dorvendo
        gui.RVEndoX(:,:,end+1) = mean(mean(gui.RVEndoX(:,:,1)));
        gui.RVEndoY(:,:,end+1) = mean(mean(gui.RVEndoY(:,:,1)));
        gui.RVEndoZ(:,end+1) = gui.RVEndoZ(:,end)+(SET(NO).SliceThickness+SET(NO).SliceGap);
      end
      if dorvepi
        gui.RVEpiX(:,:,end+1) = mean(mean(gui.RVEpiX(:,:,1)));
        gui.RVEpiY(:,:,end+1) = mean(mean(gui.RVEpiY(:,:,1)));
        gui.RVEpiZ(:,end+1) = gui.RVEpiZ(:,end)+(SET(NO).SliceThickness+SET(NO).SliceGap);
      end
    end
  end

%------- Draw
  %Create endocardium
  DATA.Handles.endosurf = surf(...
    squeeze(gui.EndoX(:,t,:)),... %EH:
    squeeze(gui.EndoY(:,t,:)),... %EH:
    gui.EndoZ,... %ZData
    'FaceAlpha','flat',...
    'AlphaDataMapping','none',...
    'AlphaData',ones([size(gui.EndoX,1) size(gui.EndoX,3)]),...
    'FaceColor','red');
  %32*ones(size(EndoX,1),size(EndoX,3)));
    %'facecolor','red');
  hold on; %Add more things :-)

  %Create epicardium
  if doepi
    DATA.Handles.episurf = surf(...
      squeeze(gui.EpiX(:,t,:)),... %EH:
      squeeze(gui.EpiY(:,t,:)),... %EH:
      gui.EpiZ,... %ZData
      'FaceAlpha','flat',...
      'AlphaDataMapping','none',...
      'AlphaData',repmat(0.5,[size(gui.EpiX,1) size(gui.EpiX,3)]),...
      'FaceColor','green');
  end  
  if dorvendo
    DATA.Handles.rvendosurf = surf(...
      squeeze(gui.RVEndoX(:,t,:)),... %EH:
      squeeze(gui.RVEndoY(:,t,:)),... %EH:
      gui.RVEndoZ,...
      'FaceAlpha','flat',...
      'AlphaDataMapping','none',...
      'AlphaData',repmat(0.5,size(gui.RVEndoZ)),...
      'FaceColor','magenta');
    hold on;
  end
  if dorvepi
    DATA.Handles.rvepisurf = surf(...
      squeeze(gui.RVEpiX(:,t,:)),... %EH:
      squeeze(gui.RVEpiY(:,t,:)),... %EH:
      gui.RVEpiZ,...
      'FaceAlpha','flat',...
      'AlphaDataMapping','none',...
      'AlphaData',repmat(0.5,size(gui.RVEpiZ)),...
      'FaceColor','cyan');
    hold on;
  end
  
  %Create annotation points
  if dopoints
    DATA.Handles.annotationpoints3d = plot3(...
      gui.PointX, ...
      gui.PointY, ...
      gui.PointZ, ...
      'bo');
    hold on
  end
  
  %Create measurements
  if domeasures
    DATA.Handles.measurements3d = zeros(1,numel(gui.MeasureX));
    for mloop = 1:size(gui.MeasureX,2)
      DATA.Handles.measurements3d(mloop) = plot3(...
        gui.MeasureX{mloop}, ...
        gui.MeasureY{mloop}, ...
        gui.MeasureZ{mloop}, ...
        'w');
    end
    hold on
  end
  
  %Create scar
  if doscar
    DATA.Handles.scar3d = plot3(...
        gui.ScarX, ...
        gui.ScarY, ...
        gui.ScarZ, ...
        'y.');
%     DATA.Handles.scarsurf = surf(...
%       gui.ScarX,... 
%       gui.ScarY,... 
%       gui.ScarZ,... 
%       'FaceAlpha','flat',...
%       'AlphaDataMapping','none',...
%       'AlphaData',repmat(0.5,[size(gui.EpiX,1) size(gui.EpiX,3)]),...
%       'FaceColor','yellow');

  end
  
  %------------------------------------
  %--- Prepare for short axis slice cut
  %------------------------------------  
  [x,y] = ndgrid(SET(NO).ResolutionX*(1:SET(NO).XSize),...
                 SET(NO).ResolutionY*(1:SET(NO).YSize));
  z = repmat((SET(NO).SliceThickness+SET(NO).SliceGap)*(SET(NO).CurrentSlice-1),SET(NO).XSize,SET(NO).YSize);
  DATA.Handles.lvcutshort = surf(...
    x,...
    y,...
    -z,...
    32*double(SET(NO).IM(:,:,t,SET(NO).CurrentSlice)),...
    'facecolor','interp',...
    'facelighting','none',...
    'edgecolor','none');
  set(DATA.Handles.lvcutshort,'cdata',...
    repmat(double(SET(NO).IM(:,:,t,SET(NO).CurrentSlice)),[1 1 3]));
  set(DATA.Handles.lvcutshort,'CDataMapping','direct');
  hold on;
  
  %Prepare for long axis slice cut
  [x,y] = ndgrid(SET(NO).ResolutionX*(1:SET(NO).XSize),...
                (SET(NO).SliceThickness+SET(NO).SliceGap)*((1:SET(NO).ZSize)-gui.basalslice));
  z = round(SET(NO).XSize/2);
  DATA.Handles.lvcutlong = surf(...
    x,...
    repmat(SET(NO).ResolutionX*z,size(x)),...
    -y,...
    32*double(squeeze(SET(NO).IM(:,z,t,:))),...
    'facecolor','interp',...
    'facelighting','none',...
    'edgecolor','none');
  set(DATA.Handles.lvcutlong,'cdata',...
    repmat(squeeze(double(SET(NO).IM(:,z,t,:))),[1 1 3]));
  set(DATA.Handles.lvcutlong,'CDataMapping','direct');
  hold off;

  %Fix colors
  %colormap([gray(64);jet(64)]);
  colormap(gray(64));
  %set(4,'mincolormap',128);
  %49 = 33+16, was 97, 65+32=97
  %set(DATA.Handles.endosurf,...
  %  'Cdata',repmat(97,size(EndoX,1),size(EndoX,3)));

  axis image vis3d off;

  report3dmodel('update');
  cameratoolbar(plotfig,'setmode','orbit');
  %rotate3d on;
%--------------------------------------------------------- End Initialize
