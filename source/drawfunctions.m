function [varargout] = drawfunctions(varargin)
% Functions for drawing in panels

% Klas

%#ok<*GVMIS>

%Invoke subfunction

global DATA 
if DATA.Autoloader || DATA.Batch
  return;
end
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%-----------------------
function drawthumbnails(calculatepreview,sliderupdated) 
%-----------------------
%Draw all thumbnails. Calculatepreview is a boolean
%indicating if thumbnails needs to be redrawn.

global DATA SET
persistent setlength

if DATA.Silent || (~DATA.DataLoaded) || isempty(DATA.Handles.datasetaxes)
  return;
end

if isempty(DATA.VisibleThumbnails)
  DATA.VisibleThumbnails = 1:min(DATA.Pref.NumberVisibleThumbnails,length(SET));
  setlength = length(SET);
end

thumbsize = DATA.GUISettings.ThumbnailSize;
if nargin<2 ||(nargin==2 && not(sliderupdated))
  segment('updateslider');
end

if nargin<1
  calculatepreview=0;
end

if setlength~=length(SET)
  calculatepreview=1;
end
setlength=length(SET);

if calculatepreview
  delete(DATA.Handles.datasetaxes.Children)
  calcfunctions('calcdatasetpreview');
end

DATA.Handles.datasetpreviewimage = image(DATA.DATASETPREVIEW,'parent',DATA.Handles.datasetaxes);
set(DATA.Handles.datasetpreviewimage,'ButtondownFcn',...
  sprintf('%s(''thumbnail_Buttondown'')','segment'));
%colormap(DATA.Colormap); %no longer needed, now truecolor /JU
axis(DATA.Handles.datasetaxes,'image','off','ij');
% axis(DATA.Handles.datasetaxes,)
if DATA.ShowRelevantStacksOnly
  relevantstacks = DATA.RelevantStacks;
  firstind = find(relevantstacks == DATA.VisibleThumbnails(1));
  firstind = firstind-1;
  lastind = find(relevantstacks == DATA.VisibleThumbnails(end));
else
  firstind = DATA.VisibleThumbnails(1)-1;
  lastind = DATA.VisibleThumbnails(end);
end
ylim =[firstind lastind]*thumbsize+1;
set(DATA.Handles.datasetaxes,'ylim',ylim ,'dataaspectratio', [1 1 1]);
hold(DATA.Handles.datasetaxes,'on');

%draw frames
drawthumbnailframes;
%print number of thumbnail in the upper left corner of the image
DATA.printthumbnailnumber(thumbsize);

hold(DATA.Handles.datasetaxes,'off');

%---------------------------
function drawthumbnailframes
%---------------------------
%Draw frames around thumbnail images
global DATA SET NO

if isempty(DATA.Handles.datasetaxes)
  return
end

thumbsize=DATA.GUISettings.ThumbnailSize;
try
  delete(DATA.Handles.datasetpreviewline);
  delete(DATA.Handles.datasetflowline);
catch %#ok<CTCH>
end

hold(DATA.Handles.datasetaxes,'on');
%Make frames used for linked images
DATA.Handles.datasetflowline = [];
linkedtono = SET(NO).Linked(SET(NO).Linked ~= NO);
if ~isempty(linkedtono)
  alllinked = unique([SET(DATA.VisibleThumbnails).Linked]);
  ind = ismember(alllinked,linkedtono);
  visibleandlinked = alllinked(ind);
  if ~isempty(visibleandlinked)
    for loop =  visibleandlinked
      if DATA.ShowRelevantStacksOnly
        numpos = find(DATA.RelevantStacks == loop);
      else
        numpos = loop;
      end
      if ~isempty(numpos)

        ypos = (numpos-1)*thumbsize+1;
        DATA.Handles.datasetflowline =  ...
          [DATA.Handles.datasetflowline ...
          plot(DATA.Handles.datasetaxes,...
          [1    1                thumbsize        thumbsize 1   ],...
          [ypos ypos+thumbsize-1 ypos+thumbsize-1 ypos      ypos],...
          'color',DATA.GUISettings.ThumbFlowLineColor)...
          ];
      end
    end
  end
end

%draw frame around current image
if DATA.ShowRelevantStacksOnly
  ind = find(DATA.RelevantStacks == NO);
else
  ind = NO;
end
ypos = (ind-1)*thumbsize+1;
DATA.Handles.datasetpreviewline =  plot(DATA.Handles.datasetaxes,...
  [1    1                thumbsize        thumbsize 1   ],...
  [ypos ypos+thumbsize-1 ypos+thumbsize-1 ypos      ypos],...
  'color',DATA.GUISettings.ThumbLineColor);

%---------------------
function killhandles 
%-----------------------
%kill the graphical handles
global DATA

if ~isfield(DATA.Handles,'imageaxes')
  return
end

graphics = [];
for i = 1:length(DATA.Handles.imageaxes)
  graphics = [graphics;DATA.Handles.imageaxes(i).Children];
end

delete(graphics)
delete(DATA.Handles.boxaxes)
delete(DATA.Handles.imageaxes)

%------------------------
function setxynan(panels)
%------------------------
%this function hides all overlay and axes by nan setting all objects.
global DATA

%Get global graphics handles
handles = DATA.Handles;
haxes = handles.imageaxes;
%all imageaxes needs to be checked.
if nargin == 0
  panels = 1:length(haxes);
end

%If a panel contains any image information then it needs to be cleared as
%it is likely that there should be a new stack in that panel after the
%configuration of the panels.
graphics = [];
for i = panels
  if ~isempty(handles.imagehandles(i).CData)
    graphics = [graphics;haxes(i).Children]; %#ok<AGROW>
  end
end

%nan set all objects in panel image is always last so fish by setting it as
%empty
for i = 1:length(graphics)
  if strcmp(graphics(i).Type,'contour')
    set(graphics(i),'ZData',[1,0;0,0],'Visible','off'); %Visible needed in 3DPrint due to Matlab bu
  elseif strcmp(graphics(i).Type,'image')
    set(graphics(i),'CData',[]);
  elseif isprop(graphics(i),'XData')
    set(graphics(i),'XData',nan,'YData',nan);
  elseif isprop(graphics(i),'Position')
    set(graphics(i),'Position',[nan nan]);
  end
end
        
%---------------------
function drawtext(panel) 
%-----------------------
%display the text in the main GUI

global DATA SET

if isempty(panel)
  return
end

if ~DATA.issegment3dp
  %--- Not 3DP
  no = DATA.ViewPanels(panel);
  if no == 0
    % this can happen when panel is empty/ does not have an image (yet)
    return
  end
  slice = [];
  
  if ismember(DATA.ViewPanelsType{panel},{'one','montage','montagerow'})
    slice = SET(no).CurrentSlice;
  end

  %we do not need to do the full update here it is sufficient if we update
  %the current image and corresponding lines in the others

  dicomimtype = SET(no).DICOMImageType;
  seriesdesc = SET(no).SeriesDescription;
  ctf = SET(no).TimeVector(SET(no).CurrentTimeFrame);
  stri = sprintf('%s\n',dprintf('Stack #%d',no));

  if ~isempty(dicomimtype)
    stri = [stri sprintf('%s\n',dicomimtype)];
  end
  if ~isempty(seriesdesc)
    stri = [stri sprintf('%s\n',seriesdesc)];
  end
  if ~isnan(ctf)
    stri = [stri sprintf('%s: %.0f ms (%d / %d)\n',dprintf('Time'),1000*ctf,SET(no).CurrentTimeFrame,SET(no).TSize)];
  end

  if (SET(no).ZSize>1) && ~isempty(slice)
    stri = [stri sprintf('%s: %d / %d\n',dprintf('Slice'),slice,SET(no).ZSize)];
  end

  venc = SET(no).VENC;
  if ~isempty(venc) && not(venc==0)
    stri = [stri sprintf('VENC: %d cm/s',venc)];
  end

  commentlist = SET(no).Comment;
  if ~isempty(commentlist)
    stri = [stri dprintf('Image comment: ')];
    for nbrcomment = 1 : length(commentlist)
      stri = [stri sprintf('%s (%s)\n', commentlist(nbrcomment).Text, commentlist(nbrcomment).Username)]; %#ok<AGROW>
    end
  end

  set(DATA.Handles.text(panel),'String',stri);

else
  %--- 3DP
  doshow = DATA.Handles.tabholder.isiconindented('showtexticon');

  if doshow && (panel <= length(DATA.ViewPanels))
    %--- We should show it, compute the string
    no = DATA.ViewPanels(panel);
    [rmax,gmax,bmax] = segment3dp.tools('xyz2rgb',SET(no).XSize,SET(no).YSize,SET(no).ZSize);
    slice = [];
    view = '';
    numslices = SET(no).ZSize;
    switch DATA.ViewPanelsType{panel}
      case 'trans3DP'
        slice = SET(no).LevelSet.View.RSlice;
        numslices = rmax;
        view = 'Transversal';
      case 'sag3DP'
        slice = SET(no).LevelSet.View.GSlice;
        numslices = gmax;
        view = 'Sagittal';
      case 'cor3DP'
        slice = SET(no).LevelSet.View.BSlice;
        numslices = bmax;
        view = 'Coronal';
      case 'speedim'
        color = SET(no).LevelSet.Pen.Color;
        switch color
          case 'r'
            slice = SET(no).LevelSet.View.RSlice;
            numslices = rmax;
            view = 'Transversal / Speed';
          case 'g'
            slice = SET(no).LevelSet.View.GSlice;
            numslices = gmax;
            view = 'Sagittal / Speed';
          case 'b'
            slice = SET(no).LevelSet.View.BSlice;
            numslices = bmax;
            view = 'Coronal / Speed';
        end
    end

    %Get data
    seriesdesc = SET(no).SeriesDescription;

    if isfield(SET(no),'IMBackup') && ~isempty(SET(no).IMBackup)
      enhanced = true;
    else
      enhanced = false;
    end

    %Start
    if enhanced
      stri = [dprintf('Enhanced') newline];
    else
      stri = '';
    end

    if ~isempty(view)
      stri = [stri sprintf('View:%s\n',view)];
    end

    if ~isempty(seriesdesc)
      stri = [stri sprintf('Ser:%s\n',seriesdesc)];
    end

    if (numslices>1) && ~isempty(slice)
      stri = [stri sprintf('%s: %d / %d\n',dprintf('Slice'),slice,numslices)];
    end

    commentlist = SET(no).Comment;
    if ~isempty(commentlist)
      stri = [stri dprintf('Image comment: ')];
      for nbrcomment = 1 : length(commentlist)
        stri = [stri sprintf('%s (%s)\n', commentlist(nbrcomment).Text, commentlist(nbrcomment).Username)]; %#ok<AGROW>
      end
    end
  else
    stri = '';
  end

  %--- Prepare to update
  if doshow
    visstate = 'on';
  else
    visstate = 'off';
  end

  %Update
  set(DATA.Handles.text(panel),'String',stri,'Visible',visstate);

end


%---------------------
function drawmeasures(panel) 
%-----------------------
%draw measurements in the main GUI
global DATA SET

no = DATA.ViewPanels(panel);

scale = viewfunctions('getscale',panel);
slices = viewfunctions('slicesinpanel',panel);
if not(isempty(SET(no).Measure))
  [measure,slice] = viewfunctions('getmeasurecoords',panel);
  x = nan;
  y = nan;
  xo = nan; %used to set other points not in this slice (measurementoutsideplane)
  yo = nan;
  
  %This makes measurement text rendering compact and index will depend on
  %number of measurements in the current timeframe and slice
  if DATA.Pref.BackgroundColor
    bgcolor = 'k';
  else
    bgcolor = 'none';
  end
  textcounter = 1;
  nummeasures = length(SET(no).Measure);
  for loop = 1: nummeasures
    ziv = round(measure(loop).Z);
    ziv = min(ziv):max(ziv);
    if any(ismember(slices,ziv)) && ...
        (SET(no).Measure(loop).T==SET(no).CurrentTimeFrame)||(isnan(SET(no).Measure(loop).T))
      
      %sometimes the column format is not used
      if ~iscolumn(measure(loop).Y)
        measure(loop).Y = measure(loop).Y';
        measure(loop).X = measure(loop).X';
      end
      
      if ~all(ismember(ziv , slices))
        DATA.Handles.measurement(panel).LineStyle = '--';
      else
        DATA.Handles.measurement(panel).LineStyle = '-';        
      end
      
      %if the current view is montage we need to translate the
      %contour to the correct slice if the entire      

      if contains(DATA.ViewPanelsType{panel},'montage')   
        for i = 1:length(ziv)
          [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(ziv(i)==slices,1));
          imdim = zoomfunctions.getxysize(no,panel);
          xt = (xl-1)*imdim.XSize - imdim.XStart +1;
          yt = (yl-1)*imdim.YSize - imdim.YStart +1;
          contoury = measure(loop).Y;
          contourx = measure(loop).X;
          [contourx,contoury] = zoomfunctions.maskoutzoomedcontour(no,panel,contourx,contoury);
          if ~any(isnan(contourx)) && ~any(isnan(contoury))
            x = [x;nan;contoury+yt];
            y = [y;nan;contourx+xt];
          else
          end
        end
      else
        imdim = zoomfunctions.getxysize(no,panel);
        contourx = scale*(measure(loop).X);
        contoury = scale*(measure(loop).Y);
        x = [x;nan;contoury];
        y = [y;nan;contourx];
        curslice = SET(no).CurrentSlice; %get current slice
        logind = (round(measure(loop).Z) ~= curslice); %logical index to not in this slice
        if sum(logind)>0
          xo = [xo;nan;scale*(measure(loop).Y(logind))];
          yo = [yo;nan;scale*(measure(loop).X(logind))];          
        end
        xt = 0;
        yt = 0;
      end
      
      if ~DATA.Run
        if ~all(isnan(x)) && ~all(isnan(y))
          if ~isfield(SET(no).Measure(loop),'Offset') || isempty(SET(no).Measure(loop).Offset)
            offset = [10 0];
          else
            offset = SET(no).Measure(loop).Offset;
          end

          pos1 = measure(loop).Y(end)+offset(1);
          pos2 = measure(loop).X(end)+offset(2);
          % Ensure position inside the image boundaries
          lengthboder = 25;
          heightborder = 10;
          pos1 = max(imdim.YStart+2, min(imdim.YEnd-lengthboder, pos1)); 
          pos2 = max(imdim.XStart+heightborder, min(imdim.XEnd-heightborder, pos2));

          SET(no).Measure(loop).Offset(1) = pos1 - measure(loop).Y(end);
          SET(no).Measure(loop).Offset(2) = pos2 - measure(loop).X(end);

          textpos = scale*[pos1+yt  pos2+xt];
          
          DATA.Handles.measurementtext(panel,textcounter).Position = textpos;
          DATA.Handles.measurementtext(panel,textcounter).String = sprintf('%s\n%0.1f [mm]',SET(no).Measure(loop).Name,SET(no).Measure(loop).Length);
          DATA.Handles.measurementtext(panel,textcounter).BackgroundColor = bgcolor;
          set(DATA.Handles.measurementtext(panel,textcounter), ...
            'ButtonDownFcn',@(hObject, event)buttondownfunctions('measuretextbuttondown',panel,hObject), ...
            'UserData',textcounter,'Tag','measurementtext');
          
          set(DATA.Handles.measurementtextline(panel,textcounter),'XData',[scale*(measure(loop).Y(end)+yt), textpos(1)],'YData',[scale*(measure(loop).X(end)+xt), textpos(2)]);
          textcounter = textcounter + 1;
        end
      end
    end
  end
  
  if ~DATA.Run
    set(DATA.Handles.measurementtext(panel,textcounter:end),...      
      'Position',[nan nan]);
     set(DATA.Handles.measurementtextline(panel,textcounter:end),...
      'XData',nan,'YData',nan);
  end
  DATA.Handles.measurement(panel).XData = x;
  DATA.Handles.measurement(panel).YData = y;
  DATA.Handles.measurementoutsideplane(panel).XData = xo;
  DATA.Handles.measurementoutsideplane(panel).YData = yo;
  
end


%---------------------
function draworthoanglehandle(panel) 
%-----------------------
%This function adds the orthoangle handle to the ortho panel
global DATA SET

scale = viewfunctions('getscale',panel);
no = DATA.ViewPanels(panel);

if isnan(DATA.Handles.orthoanglehandle.XData)
  %default is quarter from the image border
  [x,y] = calcfunctions('calcplaneintersections',no,no,'orth','gla');
  
  %this is the intersection point of all planeintersections
  xc = SET(no).HLA.slice;
  yc = SET(no).VLA.slice;
  
  %then we find in which direction on the gla line there is the most
  %space and place the anglehandle point here.
  val1 = norm([xc,yc]-[x(1),y(1)]);
  val2 = norm([xc,yc]-[x(2),y(2)]);
  
  if val1<val2
    x = xc + (xc-x(1))/2;
    y = yc + (yc-y(1))/2;
  else
    x = xc + (xc-x(end))/2;
    y = yc + (yc-y(end))/2;
  end
  
  DATA.Handles.orthoanglehandle(panel).Parent = DATA.Handles.imageaxes(panel);
  DATA.Handles.orthoanglehandle(panel).YData = scale*x;%[nan y(i([1 end]))*[0.3;0.7]];
  DATA.Handles.orthoanglehandle(panel).XData = scale*y;%[nan x(i([1 end]))*[0.3;0.7]];
end

%------------------------
function drawpanel(panel)
%------------------------
%does what it's told from drawlist
global DATA

list = DATA.drawlist;
if isempty(list)
  return
end

if (~DATA.Silent) && (~isempty(panel))
  for i = 1:length(list{panel})
    fcn = list{panel}{i};
    if isa(fcn,'function_handle')
      fcn();
    else
      eval(fcn);
    end
  end
end

%---------------------
function drawroi(panel,colortypes) 
%-----------------------
%draw the ROI in the main GUI
global DATA SET

no = DATA.ViewPanels(panel);

if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
  originalno = no;% no actually shown in panel
  no = SET(no).Flow.MagnitudeNo;
  % resetting that current time frame is always the one that is shown in
  % the panel
  SET(no).CurrentTimeFrame = SET(originalno).CurrentTimeFrame;
end

scale = viewfunctions('getscale',panel);
slicestoinclude = viewfunctions('slicesinpanel',panel);

%Define colortypes. These are then used to concatenate into the correct
%index in the x,y variables which finally are dealt into the different
%color roi graphics objects.'
% % % textreset = false;
if nargin == 1%colortypes
  colortypes = 'cgbrwkym';
% % %   set(DATA.Handles.roitext(panel,:),'Position',[nan nan]);
% % %   textreset = true;
end

x = cell(1,length(colortypes));
y = cell(1,length(colortypes));
% % % maxnumroitxt = 12;
%only work with rois that exist in the included slices
numroisinslice = nnz(cellfun(@(x,y) any(~isempty(x)) && any(x==slicestoinclude),{SET(no).Roi.Z}));
roistodo = find(cellfun(@(x,y,z) numroisinslice ~= 0 && any(~isempty(x)) && any(x==slicestoinclude) && ismember(y(1),colortypes) && any(~isempty(z)),{SET(no).Roi.Z},{SET(no).Roi.LineSpec},{SET(no).Roi.Area}));
% % % if ~textreset
% % %    if numroisinslice <= maxnumroitxt 
% % %     if length(DATA.Handles.roitext(panel,:))>= numroisinslice
% % %       %reset all roitexts that are left after using ctrl-z
% % % %       set(DATA.Handles.roitext(panel,:),'Position',[nan nan]);
% % %     else
% % %       set(DATA.Handles.roitext(panel,1:numroisinslice),'Position',[nan nan]);
% % %     end 
% % %   else
% % %     set(DATA.Handles.roitext(panel,:),'Position',[nan nan]);
% % %   end
% % % end
linewidth = DATA.Pref.LineWidth;
if isempty(linewidth)
  linewidth = 1;
end
% % % if DATA.Pref.BackgroundColor
% % %   bgcolor = 'k';
% % % else
% % %   bgcolor = 'none';
% % % end
%KG: 
for loop = 1:length(roistodo)
  %KG:
  [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slicestoinclude==SET(no).Roi(roistodo(loop)).Z));
  imdim = zoomfunctions.getxysize(no,panel);

  xt = (xl-1)*imdim.YSize - imdim.YStart +1;
  yt = (yl-1)*imdim.XSize - imdim.XStart +1;

  origxdata = SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame);
  margin = 0.5;
  origxdata(origxdata > imdim.YEnd+margin) = nan;
  origxdata(origxdata < imdim.YStart-margin) = nan;

  origydata = SET(no).Roi(roistodo(loop)).X(:,SET(no).CurrentTimeFrame);
  origydata(origydata > imdim.XEnd+margin) = nan;
  origydata(origydata < imdim.XStart-margin) = nan;

  if roistodo(loop) == SET(no).RoiCurrent    
    DATA.Handles.roicurrent(panel).XData = scale*(xt+origxdata);
    DATA.Handles.roicurrent(panel).YData = scale*(yt+origydata);
    if strcmp(SET(no).Roi(roistodo(loop)).LineSpec(1),'y')
      color = get(DATA.Handles.yroi,'color'); %yellow is actually orange
      color = color{1};
    else
      color = SET(no).Roi(roistodo(loop)).LineSpec(1);
    end
    DATA.Handles.roicurrent(panel).Color = color;
    DATA.Handles.roicurrent(panel).LineWidth = linewidth + 1;
  else 
    cind = regexp(colortypes,SET(no).Roi(roistodo(loop)).LineSpec(1));
    x{cind} = cat(1,x{cind},nan,xt+origxdata);
    y{cind} = cat(1,y{cind},nan,yt+origydata);
    DATA.Handles.([colortypes(cind),'roi'])(panel).LineWidth = linewidth;
  end
end
if nargin  == 1
  drawroitext(panel);
end
% % % if ~DATA.Run
% % %   %We only display text in one and orth viewmode
% % %   if any(strcmp(DATA.ViewPanelsType{panel},{'one','orth'})) || any(strcmp(DATA.ViewPanelsType{panel},{'montage','montagerow'}))
% % % 
% % %     for loop = 1:length(roistodo)
% % %       %if ~all(isnan(SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame)))
% % %       %get sign for text
% % %       %KG:
% % %       if SET(no).Roi(roistodo(loop)).Sign > 0
% % %         roisign = '';
% % %       else
% % %         roisign = ' (-)';
% % %       end
% % %       % % %       cind = regexp(colortypes,SET(no).Roi((loop)).LineSpec(1));
% % % 
% % %       %ROI position
% % % 
% % %       % we also need to check if the roitext is within the panel limits
% % %       % otherwise we dont plot the text
% % %       %[ymin,ix] = min(SET(no).Roi((loop)).Y(:,SET(no).CurrentTimeFrame));
% % %       if numroisinslice <= maxnumroitxt
% % %         [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slicestoinclude==SET(no).Roi(roistodo(loop)).Z));
% % %         imdim = zoomfunctions.getxysize(no,panel);
% % % 
% % %         xt = (xl-1)*imdim.YSize - imdim.YStart +1;
% % %         yt = (yl-1)*imdim.XSize - imdim.XStart +1;
% % %         %This checks if it is a flow roi we are dealing with (do we have the output generated from calcflow in all timeframes). If so
% % %         %do not plot the Std and Mean intensities.
% % %         if ~isempty(SET(no).Flow)
% % %           [ymin,ix] = min(SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame));
% % %           ypos = ymin-1;
% % %           xpos = SET(no).Roi(roistodo(loop)).X(ix,SET(no).CurrentTimeFrame);
% % %           halign = 'right';
% % %           if isempty(SET(no).Roi(roistodo(loop)).Area)
% % %             labelstr = sprintf('%s%s',SET(no).Roi(roistodo(loop)).Name,roisign);
% % %           else
% % %             labelstr = sprintf('%s%s\n%3.1f [cm^2]',...
% % %               SET(no).Roi(roistodo(loop)).Name,roisign, ...
% % %               SET(no).Roi(roistodo(loop)).Area(SET(no).CurrentTimeFrame));
% % %           end
% % %         else
% % %           ypos = mean(SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame));
% % %           xpos = mean(SET(no).Roi(roistodo(loop)).X(:,SET(no).CurrentTimeFrame));
% % %           halign = 'left';
% % %           labelstr = sprintf('%s%s\n%3.1f [cm^2]\n%3.1f +/- %3.1f',...
% % %             SET(no).Roi(roistodo(loop)).Name,roisign, ...
% % %             SET(no).Roi(roistodo(loop)).Area(SET(no).CurrentTimeFrame), ...
% % %             SET(no).Roi(roistodo(loop)).Mean(SET(no).CurrentTimeFrame), ...
% % %             SET(no).Roi(roistodo(loop)).StD(SET(no).CurrentTimeFrame));
% % %         end
% % %         texthandles = DATA.Handles.roitext(panel,:);
% % %         pos = {texthandles.Position};
% % %         textind = find(cellfun(@(x) any(isnan(x(:))), pos));
% % %         if ~isempty(textind)
% % %           textloop = textind(1);
% % %           DATA.Handles.roitext(panel,(textloop)).Position = scale*[ypos+xt xpos+yt];
% % %           DATA.Handles.roitext(panel,(textloop)).HorizontalAlignment = halign;
% % %           DATA.Handles.roitext(panel,(textloop)).String = labelstr;
% % %           DATA.Handles.roitext(panel,(textloop)).BackgroundColor = bgcolor;
% % %         end        
% % %       end
% % %     end
% % %   else
% % %     %set the different colortypes
% % %     set(DATA.Handles.roitext(panel,:),'Position',[nan nan]);
% % %   end
% % % end

%set the different colortypes
for i = 1:length(colortypes)
  DATA.Handles.([colortypes(i),'roi'])(panel).XData = scale*x{i};
  DATA.Handles.([colortypes(i),'roi'])(panel).YData = scale*y{i};
end

%--------------------------
function drawroitext(panel)
%--------------------------
% function to show roi text
global DATA SET
texthandles = DATA.Handles.roitext(panel,:);
set(texthandles,'Position',[nan nan]);
if ~DATA.Run % display when not playing
  if DATA.Pref.BackgroundColor
    bgcolor = 'k';
  else
    bgcolor = 'none';
  end
  no = DATA.ViewPanels(panel);
  if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
    originalno = no;% no actually shown in panel
    no = SET(no).Flow.MagnitudeNo;
    % resetting that current time frame is always the one that is shown in
    % the panel
    SET(no).CurrentTimeFrame = SET(originalno).CurrentTimeFrame;
  end


  scale = viewfunctions('getscale',panel);
  slicestoinclude = viewfunctions('slicesinpanel',panel);
  maxnumroitext = viewfunctions('getnumberofmaxroitext',panel);
  %only work with rois that exist in the included slices
  numroisinslice = nnz(cellfun(@(x,y) any(~isempty(x)) && any(x==slicestoinclude),{SET(no).Roi.Z}));
  roistodo = find(cellfun(@(x,z) numroisinslice ~= 0 && any(~isempty(x)) && any(x==slicestoinclude) && any(~isempty(z)),{SET(no).Roi.Z},{SET(no).Roi.Area}));
  if numroisinslice <= maxnumroitext
    imageaxeshandles = DATA.Handles.imageaxes(panel);
    for loop = 1:length(roistodo)
      if SET(no).Roi(roistodo(loop)).Sign > 0
        roisign = '';
      else
        roisign = ' (-)';
      end
      [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slicestoinclude==SET(no).Roi(roistodo(loop)).Z));
      imdim = zoomfunctions.getxysize(no,panel);

      yt = (yl-1)*imdim.YSize - imdim.YStart +1;
      xt = (xl-1)*imdim.XSize - imdim.XStart +1;

      %This checks if it is a flow roi we are dealing with (do we have the output generated from calcflow in all timeframes). If so
      %do not plot the Std and Mean intensities.
      if ~isempty(SET(no).Flow)
        [ymin,ix] = min(SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame));
        ypos = ymin-1;
        xpos = SET(no).Roi(roistodo(loop)).X(ix,SET(no).CurrentTimeFrame);
        halign = 'right';
        if isempty(SET(no).Roi(roistodo(loop)).Area)
          labelstr = sprintf('%s%s',SET(no).Roi(roistodo(loop)).Name,roisign);
        else
          labelstr = sprintf('%s%s\n%3.1f [cm^2]',...
            SET(no).Roi(roistodo(loop)).Name,roisign, ...
            SET(no).Roi(roistodo(loop)).Area(SET(no).CurrentTimeFrame));
        end
      else
        ypos = mean(SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame));
        xpos = mean(SET(no).Roi(roistodo(loop)).X(:,SET(no).CurrentTimeFrame));
        halign = 'left';
        labelstr = sprintf('%s%s\n%3.1f [cm^2]\n%3.1f +/- %3.1f',...
          SET(no).Roi(roistodo(loop)).Name,roisign, ...
          SET(no).Roi(roistodo(loop)).Area(SET(no).CurrentTimeFrame), ...
          SET(no).Roi(roistodo(loop)).Mean(SET(no).CurrentTimeFrame), ...
          SET(no).Roi(roistodo(loop)).StD(SET(no).CurrentTimeFrame));
      end

      pos = {texthandles.Position};
      textind = find(cellfun(@(x) any(isnan(x(:))), pos));
      
      if ~isempty(textind)
        postxt = scale*[ypos+yt xpos+xt];
        htmp = text('Position',postxt,'String',labelstr,'Visible','off', ...
          'Parent',imageaxeshandles,'HorizontalAlignment',halign);
        txtextent = htmp.Extent;
        xlim = imageaxeshandles.XLim;
        ylim = imageaxeshandles.YLim;

        adjustedy = postxt(2) - txtextent(4) / 2;

        textfits = (postxt(1) >= xlim(1)) && ...
           (postxt(1) + txtextent(3) <= xlim(2)) && ...
           (adjustedy >= ylim(1)) && ...
           (adjustedy + txtextent(4)  <= ylim(2));
      
        if textfits
          textloop = textind(1);
          DATA.Handles.roitext(panel,(textloop)).Position = scale*[ypos+yt xpos+xt];
          DATA.Handles.roitext(panel,(textloop)).HorizontalAlignment = halign;
          DATA.Handles.roitext(panel,(textloop)).String = labelstr;
          DATA.Handles.roitext(panel,(textloop)).BackgroundColor = bgcolor;
        end
        delete(htmp)
      end
    end
  end
end

%------------------------------
function drawplaneintersections 
%------------------------------
%draw plane intersections in the main GUI
global DATA

%Get global graphics handles
handles = DATA.Handles;

panels = find(DATA.ViewPanels);
for p1 = panels %this is the nos panel which we are going to plot in
  if ~any(strcmp(DATA.ViewPanelsType{p1},{'montage','montagesegmented','montagerow'})) %dont plot plane intersections within montage views
    scale = viewfunctions('getscale',p1);    
    for p2 = panels
      if p1~= p2
        [x,y] = calcfunctions('calcplaneintersections',...
          DATA.ViewPanels(p1),DATA.ViewPanels(p2),DATA.ViewPanelsType{p1},...
          DATA.ViewPanelsType{p2});%,'one','one',DATA.slices(p1),DATA.slices(p2));

        set([handles.planeintersection(p1,p2) handles.planeintersection(p1,p2)],...
          'YData',scale*x,'XData',scale*y);
        
        if p2 == DATA.CurrentPanel
          color = 'y';
        else
          color = 'w';
        end
        set(handles.planeintersection(p1,p2),'Color',color);
      end
    end
  else
    set(handles.planeintersection(p1,:),'YData', nan, 'XData', nan);
  end
end

%----------------------
function drawspeedimage 
%----------------------
%Ensure speedimage is updated.
global DATA

drawimages(find(strcmp(DATA.ViewPanelsType,'speedim')));

%-----------------------
function im = speedimage 
%-----------------------
%generates the 3dp speed image in rgb format ready for display in segment

global DATA SET NO

%--- generate 2D 'speed' image

tempimage = getimage(SET(NO).LevelSet.Pen.Color);

%---Remap the image

bw = segment3dp.tools('calcbwfromim',NO,tempimage);

tempremapped = single(bw)+1;

%convert colormap to uint8
cmap = DATA.LevelSet.colormap;
cmap = uint8(255*cmap);
temprgb = cat(2,...
  cmap(tempremapped,1),...
  cmap(tempremapped,2),...
  cmap(tempremapped,3));

im = reshape(temprgb,[size(tempimage,1), size(tempimage,2), 3]);

%-------------------------
function drawimages(panel) 
%-------------------------
%draw images in image display in main GUI

global DATA SET

if isempty(panel)
  return
end

no = DATA.ViewPanels(panel);

if DATA.issegment3dp

  switch DATA.ViewPanelsType{panel}
    case 'trans3DP'
      im = segment3dp.graphics('getoverlayimage','r');
      set(DATA.Handles.imagehandles(panel),'CData',im,'XData',[1 size(im,2)],'YData',[1 size(im,1)]);
      if isempty(SET(no).LevelSet.View.RZoomState)
        viewfunctions('getnewzoomstate',panel);
      end

    case 'sag3DP'
      im = segment3dp.graphics('getoverlayimage','g');
      set(DATA.Handles.imagehandles(panel),'CData',im,'XData',[1 size(im,2)],'YData',[1 size(im,1)]);
       if isempty(SET(no).LevelSet.View.GZoomState)
        viewfunctions('getnewzoomstate',panel);
      end
    case 'cor3DP'
      im = segment3dp.graphics('getoverlayimage','b');
      set(DATA.Handles.imagehandles(panel),'CData',im,'XData',[1 size(im,2)],'YData',[1 size(im,1)]);
       if isempty(SET(no).LevelSet.View.BZoomState)
        viewfunctions('getnewzoomstate',panel);
      end
    case 'speedim'
      im = drawfunctions('speedimage');
      set(DATA.Handles.imagehandles(panel),'CData',im,'XData',[1 size(im,2)],'YData',[1 size(im,1)]);
    case 'viewport'
      return
  end
    
  return

else
  if panel>length(DATA.ViewIM) || isempty(DATA.ViewIM{panel})
    createfunctions('createviewim',panel);
  end
end

sz = size(DATA.ViewIM{panel});
no = DATA.ViewPanels(panel);
if numel(sz) == 5
  im = squeeze(DATA.ViewIM{panel}(:,:,SET(no).CurrentTimeFrame,:,:));
  DATA.Handles.imagehandles(panel).CData = im;
elseif (numel(sz) == 3 && sz(3) == 3) ... %RGB images
    && ~contains(SET(no).DICOMImageType,'GT\AIF')... %Fix for Gadgeton diagrams (Aera 2018)
    && ~((isfield(SET(no),'T2preptime') && ~isempty(SET(no).T2preptime) && numel(unique(SET(no).T2preptime)) == 3) ||...
    (isfield(SET(no),'EchoTime') && ~isempty(SET(no).EchoTime) && numel(unique(SET(no).EchoTime)) == 3) ||...
    (isfield(SET(no),'InversionTime') && ~isempty(SET(no).InversionTime) && numel(unique(SET(no).InversionTime)) == 3)) %Fix for T1/T2 maps
  im = DATA.ViewIM{panel};
  DATA.Handles.imagehandles(panel).CData = im;
elseif ~isempty(SET(no).Colormap)
  im = DATA.ViewIM{panel}(:,:,SET(DATA.ViewPanels(panel)).CurrentTimeFrame);  
  DATA.Handles.imagehandles(panel).CData = ind2rgb(im,SET(DATA.ViewPanels(panel)).Colormap);
else
  im = DATA.ViewIM{panel}(:,:,SET(no).CurrentTimeFrame);  
  DATA.Handles.imagehandles(panel).CData = cat(3,im,im,im);
end

%Fix for interpolated. More work is required to verify this.
updatexydata('imagehandles',panel,size(im));
  
%--------------------------------
function drawselectedframe(panel) 
%--------------------------------
%Function that highlights the selected panel

global DATA

%turn all frames off
for p = 1:length(DATA.ViewPanels)
  DATA.Handles.selectedframe(p).Visible = 'off';
end

xl = DATA.Handles.imageaxes(panel).XLim;
yl = DATA.Handles.imageaxes(panel).YLim;
DATA.Handles.selectedframe(panel).XData =[xl,fliplr(xl),xl(1)];
DATA.Handles.selectedframe(panel).YData =[yl(1),yl(1),yl(2),yl(2),yl(1)];
DATA.Handles.selectedframe(panel).Visible = 'on';

%------------------------------------------------------------------------
function drawselectedslice(panel)
%------------------------------------------------------------------------
%Function that highlights the selected slice
global DATA SET

no = DATA.ViewPanels(panel);

%if panels are in mode 'one'/'orth' or
%if user unselected all slices (start and end slices are empty),
%then no frame is shown
if any(strcmp(DATA.ViewPanelsType{panel},{'one','orth'}))||((isempty(SET(no).StartSlice))&&isempty(SET(no).EndSlice))
  DATA.Handles.selectedslice(panel).XData = nan;
  DATA.Handles.selectedslice(panel).YData = nan;
  return
end

slicestoinclude = viewfunctions('slicesinpanel',panel);

ss = find(SET(no).StartSlice==slicestoinclude);
se = find(SET(no).EndSlice==slicestoinclude);

imdim = zoomfunctions.getxysize(no,panel);
numslices = se-ss+1;
x = nan(6*numslices,1);
y = x;
indstart = 1;
indend = 5;
for zloop = ss:se
  [x1,y1] = ind2sub(DATA.ViewPanelsMatrix{panel},zloop);
  x1 = (x1-1)*imdim.YSize + 1.5;
  y1 = (y1-1)*imdim.XSize + 1.5;
  x2 = -1.5+x1+imdim.YSize;
  y2 = -1.5+y1+imdim.XSize;
  x(indstart:indend) = [x1 x2 x2 x1 x1];
  y(indstart:indend) = [y1 y1 y2 y2 y1];
  indstart = indstart+6;
  indend = indend+6;
end

DATA.Handles.selectedslice(panel).XData = x;
DATA.Handles.selectedslice(panel).YData = y;

%fix so that linked images slices are plotted
nos = SET(no).Linked;
%nos(nos==no)=[]; %remove the current panel. EH: why?
panelstodo=find(ismember(DATA.ViewPanels,nos));
for loop=1:length(panelstodo)
  panel = panelstodo(loop);
  switch DATA.ViewPanelsType{panel}
    case {'montage','montagerow','montagefit','montagesegmented'}
      DATA.Handles.selectedslice(panel).XData = x;
      DATA.Handles.selectedslice(panel).YData = y;
      DATA.Handles.selectedslice(panel).LineStyle = '--';
      DATA.Handles.selectedslice(panel).Visible = 'on';
  end
end

%----------------------------
function drawviability(panel) 
%----------------------------
%Function to draw viability contour on screen used from drawimageslice,
%drawimagemontage
global DATA SET

hidecell={'hidescar'};
stateandicon = viewfunctions('iconson',hidecell);
state=[stateandicon{:,1}];
if state(1)
  return
end

type = DATA.ViewPanelsType{panel};
scale = viewfunctions('getscale',panel);

no = DATA.ViewPanels(panel);
%Get mask
switch type
  case {'montage','montagerow','montagefit','sax3'}
    result = viewfunctions('reshape2layout',SET(no).Scar.Result,panel);
  case {'one','orth'}
    result = SET(no).Scar.Result(:,:,SET(no).CurrentSlice);
  otherwise
    myfailed('Unknown viewtype in drawviabilityhelper.');
end

if sum(result(:))>0
  
  %Normal scar outline
  scarim = imresize(double(result),scale,'bilinear');
  DATA.Handles.scarcontour(panel).ZData = scarim;
  DATA.Handles.scarcontour(panel).Visible = 'on';
  %Fix for interpolated
  updatexydata('scarcontour',panel,size(scarim))
  %Weighted scar outline
  if SET(no).Scar.UseWeighting
    
    weighting = viability('viabilityweight',no); %This is a vector with weight and position as corresponding pixels in the whole volume.
    
    if ~isempty(weighting)
      %Get weighting
      temp = zeros(size(SET(no).Scar.Result));
      
      temp(SET(no).Scar.Result) = weighting;
      temp(SET(no).Scar.NoReflow) = 1; %Make sure MO is always included weighted graphically.
      
      %Sort the weighting
      sortedweighting = sort(weighting);
      weighted = sum(weighting(:)); %Calculate weighted percentage
      total = sum(SET(no).Scar.Result(:));
      f = weighted/total;
      
      %Find threshold so that it would include f% of the area
      f = 1-f; %Since we want the pixelse that are larger than the threshold
      thres = sortedweighting(min(max(round(length(sortedweighting)*f),1),length(sortedweighting)));
      
      %Convert to right size
      switch type
        case {'montage','montagerow','montagefit','sax3'}
          weighting = viewfunctions('reshape2layout',temp,panel);
        case {'one','orth'}
          weighting = temp(:,:,SET(no).CurrentSlice);
        otherwise
          myfailed('Unknown viewtype in drawviabilityhelper.');
      end
      if max(weighting(:)) < thres
        %no pink region to show
        DATA.Handles.weightedscarcontour(panel).ZData = [];
        DATA.Handles.weightedscarcontour(panel).XData = [];
        DATA.Handles.weightedscarcontour(panel).YData = [];
        DATA.Handles.weightedscarcontour(panel).Visible = 'off';
      else
        weightedscarim = imresize(weighting,scale,'bilinear');
        DATA.Handles.weightedscarcontour(panel).ZData = weightedscarim;
        DATA.Handles.weightedscarcontour(panel).Visible = 'on';
        DATA.Handles.weightedscarcontour(panel).LevelList = [thres thres];
        %Fix for interpolated
        updatexydata('weightedscarcontour',panel,size(weightedscarim));
      end      
    else
      DATA.Handles.weightedscarcontour(panel).ZData = [];
      DATA.Handles.weightedscarcontour(panel).XData = [];
      DATA.Handles.weightedscarcontour(panel).YData = [];
      DATA.Handles.weightedscarcontour(panel).Visible = 'off';
    end %nonempty weighting
  else
    DATA.Handles.weightedscarcontour(panel).ZData = [];
    DATA.Handles.weightedscarcontour(panel).XData = [];
    DATA.Handles.weightedscarcontour(panel).YData = [];
    DATA.Handles.weightedscarcontour(panel).Visible = 'off';
  end %use weighting
end

%Update taken mocontour and moextentcontour graphically
if sum(SET(no).Scar.NoReflow(:))>0
  
  if ~isfield(SET(no).Scar,'MOThreshold')
    SET(no).Scar.MOThreshold = 1.5;
  end
  
  %Extract and convert to right size
  %Take image and mask with no reflow.
  doit = true;
  
  switch type%DATA.ViewPanelsType{panel}
    case {'montage','montagerow','montagefit','sax3'}
      contourim = SET(no).Scar.IM;
      contourim(~SET(no).Scar.NoReflow) = 1; %Mask with no reflow. 1 is max.
      contourim = viewfunctions('reshape2layout',contourim,panel,NaN); %NaN is outsideelement.
    case {'one','orth'}
      if existfunctions('anyall',SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice))
        contourim = SET(no).Scar.IM(:,:,SET(no).CurrentSlice);
        contourim(~SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice)) = 1; %1 is max
      else
        doit = false;
      end
    otherwise
      myfailed('Unknown viewtype in drawviabilityhelper.');
  end
  
  if (not(isequal(length(SET(no).Scar.mthreshold),SET(no).ZSize)||isequal(length(SET(no).Scar.mthreshold),1)))&& doit
    %Length is not correct.
    viability('viabilitycalc'); %this will updat the size correctly.
  end
  
  if doit
    %If we need to do, then doit.
    
    %--- Find suitable threshold
    
    if isequal(SET(no).Scar.mthreshold,0)
      %Manual mode
      mthreshold = 1;
      logslice = SET(no).Scar.Result(:,:,SET(no).CurrentSlice);
      logslice = logslice & (~SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice));
      if sum(logslice(:))>0
        im = SET(no).Scar.IM(:,:,SET(no).CurrentSlice);
        mthreshold = min(im(logslice(:)));
        if mthreshold<0.01
          mthreshold=0.01;
        end
      end
    else
      %Some automated mode
      mthreshold = SET(no).Scar.mthreshold;
    end
    
    %Use remote intensity to calculate good threshold
    if length(mthreshold)<2
      thres = mthreshold;
    else
      thres = mthreshold(SET(no).CurrentSlice);
    end
    thres = thres*SET(no).Scar.MOThreshold; %OBS if you change here, you need ALSO to change in viabilitycalcvolume.
    mocoim = imresize(contourim,scale,'bilinear');
    DATA.Handles.mocontour(panel).ZData = mocoim;
    DATA.Handles.mocontour(panel).Visible = 'on';
    DATA.Handles.mocontour(panel).LevelList = [thres thres];
    DATA.Handles.moextentcontour(panel).ZData = imresize(contourim,scale,'bilinear');
    DATA.Handles.moextentcontour(panel).Visible = 'on';
    %Fix for interpolated
    szim = size(mocoim);
    updatexydata('mocontour',panel,szim);
    updatexydata('moextentcontour',panel,szim);    
  end
end

%----------------------
function drawmar(panel) 
%----------------------
%Function to draw MaR contour on screeen used from drawimageslice,
%drawimagemontage

global DATA SET
no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

%Get mask
switch DATA.ViewPanelsType{panel}
  case {'montage','montagerow','montagefit','sax3'}
    result = viewfunctions('reshape2layout',squeeze(SET(no).MaR.Result(:,:,SET(no).CurrentTimeFrame,:)),panel);
  case {'one','orth'}
    result = SET(no).MaR.Result(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
  otherwise
    myfailed('Unknown viewtype in mar-drawhelper.');
end
%Normal MaR outline
%Contour rendering not triggered for constant Zdata
if all(result(:) == 0)
  result(1) = 1;
end
marimg = imresize(double(result),scale,'bilinear');
DATA.Handles.marcontour(panel).ZData = marimg;
updatexydata('marcontour',panel,size(marimg));
stateandicon = viewfunctions('iconson','hidemar');
if ~stateandicon{1}
  DATA.Handles.marcontour(panel).Visible = 'on';
end
DATA.Handles.marcontour(panel).LevelList = [0.5 0.5];

%------------------------------
function drawcentercross(panel) 
%------------------------------
%draw the center cross in the image display

global SET DATA

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);
slices = viewfunctions('slicesinpanel',panel);

if isequal(DATA.ViewPanelsType{panel},'viewport')
  return
end

[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},slices-slices(1)+1);

imdim = zoomfunctions.getxysize(no,panel);

origxdata = SET(no).CenterY;
origxdata(origxdata > imdim.YEnd) = nan;
origxdata(origxdata < imdim.YStart) = nan;

origydata = SET(no).CenterX;
origydata(origydata > imdim.XEnd) = nan;
origydata(origydata < imdim.XStart) = nan;

xt = ((xl-1)*imdim.XSize- imdim.XStart +1);
yt = ((yl-1)*imdim.YSize- imdim.YStart +1);

DATA.Handles.centercross(panel).XData = scale*(origxdata+yt);
DATA.Handles.centercross(panel).YData = scale*(origydata+xt);

%------------------------------
function drawno(no, isonlyannotation) 
%------------------------------
%This function is intended to be used after calls that alter a no and may
%effect multiple panels and linked images. This will be used in all
%callbacks altering contours and rois for now. If more speed is needed then we will
%need to write specific calls to draw functions at each callback.

global DATA SET NO

if DATA.Silent
  return
end

if nargin < 1
  no = NO;
end

if nargin < 2
  isonlyannotation = false;
end

%get all panels which dont have no in it or are linked these are the ones
%that possible have intersections with the updated no and therefore needs
%updating.
panelstodo = find(~ismember(DATA.ViewPanels,SET(no).Linked) & DATA.ViewPanels>0);
panelslinked = find(ismember(DATA.ViewPanels,SET(no).Linked));

%Then we udpate intersections in intersecting panels
for p = panelstodo
  drawintersections(p)
end

%update altered panels fully
for p = panelslinked
  if strcmp(DATA.ViewPanelsType{p},'one') || strcmp(DATA.ViewPanelsType{p},'orth')
    DATA.ViewIM{p} = [];
  end
  viewfunctions('updatedrawlist',p,isonlyannotation)
  drawpanel(p)
  if isonlyannotation
    % add "drawimage" again, since it was skipped
    viewfunctions('addtodrawlist',p,sprintf('drawfunctions(''drawimages'',%d)',p));
  end
  if ~DATA.issegment3dp && (length(DATA.ViewIM)>=p) && isempty(DATA.ViewIM{p})
    createfunctions('createviewim',p);
  end
end

%if montage then we should update the selected slice.
if any(strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},{'montage','montagesegmented','montagerow'}))
  drawselectedslice(DATA.CurrentPanel)
end

%------------------------------
function drawpoint3D(panel) 
%------------------------------
%this draws and updates the panel image intersections.
global DATA NO
drawfunctions('drawplaneintersections')

if nargin == 1
  % add points to all panels
  scale = viewfunctions('getscale',panel);
  switch DATA.ViewPanelsType{panel}
    case 'hla'
      [y,z,x] = segment('getclickedcoords');
      DATA.Handles.point3D(panel).XData(end+1) = scale*y;
      DATA.Handles.point3D(panel).YData(end+1) = scale*z;
    case 'vla'
      [x,z,y] = segment('getclickedcoords');
      DATA.Handles.point3D(panel).XData(end+1) = scale*x;
      DATA.Handles.point3D(panel).YData(end+1) = scale*z;
    case 'gla'
      [ytemp,xtemp] = segment('getclickedcoords');
      [x,y,z] = calcfunctions('gla2sax',xtemp,ytemp,NO);
      [xo,yo] = calcfunctions('sax2gla',x,y,z,NO);
      DATA.Handles.point3D(panel).XData(end+1) = scale*xo;
      DATA.Handles.point3D(panel).YData(end+1) = scale*yo;
    case {'one','orth'}
      [y,x,z] = segment('getclickedcoords');
      DATA.Handles.point3D(panel).XData(end+1) = y;
      DATA.Handles.point3D(panel).YData(end+1) = x;
      x = x/scale;
      y = y/scale;
   case {'montage'}
      [y,x,z] = segment('getclickedcoords');
      [xofs,yofs] = calcfunctions('calcoffset',z,[],NO);
      DATA.Handles.point3D(panel).XData(end+1) = scale*(y+yofs);
      DATA.Handles.point3D(panel).YData(end+1) = scale*(x+xofs);
    case {'montagerow'}
      [y,x,z] = segment('getclickedcoords');
      [xofs,yofs] = calcfunctions('calcoffset',z,[],NO);
      DATA.Handles.point3D(panel).XData(end+1) = scale*(y+xofs);
      DATA.Handles.point3D(panel).YData(end+1) = scale*(x+yofs);
    otherwise
      [y,x,z] = segment('getclickedcoords');
      [xofs,yofs] = calcfunctions('calcoffset',z,[],NO);
      DATA.Handles.point3D(panel).XData(end+1) = scale*(y+yofs);
      DATA.Handles.point3D(panel).YData(end+1) = scale*(x+xofs);
  end 
  pos = calcfunctions('xyz2rlapfh',NO,x,y,z);
  for panel = find(DATA.ViewPanels ~= NO & DATA.ViewPanels > 0)
    no = DATA.ViewPanels(panel);
    scale = viewfunctions('getscale',panel);
    xyz = calcfunctions('rlapfh2xyz',no,pos(:,1),pos(:,2),pos(:,3));
    switch DATA.ViewPanelsType{panel}
      case {'one','orth'}
        newx = scale*xyz(2,:);
        newy = scale*xyz(1,:);
      case {'montage'}
        [xofs,yofs] = calcfunctions('calcoffset',xyz(3,:),[],no, panel);
        newx = scale*(xyz(2,:))+yofs;
        newy = scale*(xyz(1,:))+xofs;
      otherwise
        [xofs,yofs] = calcfunctions('calcoffset',xyz(3,:),[],no, panel);
        newx = scale*(xyz(2,:)+yofs);
        newy = scale*(xyz(1,:)+xofs);
    end
    
    DATA.Handles.point3D(panel).XData(end+1) = newx;%scale*xyz(2,:);
    DATA.Handles.point3D(panel).YData(end+1) = newy;%scale*xyz(1,:);
  end
else
  % delete all points be setting all values to nana
  for panel = find(DATA.ViewPanels > 0)
    DATA.Handles.point3D(panel).XData(:) = nan;
    DATA.Handles.point3D(panel).YData(:) = nan; 
  end
  
end

%------------------------
function drawpoint(panel) 
%------------------------
%draw annotation points in the main GUI

global DATA SET

if DATA.issegment3dp
  return
end

%Get global graphics handles
handles = DATA.Handles;
no = DATA.ViewPanels(panel);

%Clear all text
if ~DATA.Run %&& ~(length(SET(no).Point.Label) > 20) %text is displayed only for the first 20 points
  set(handles.pointtext(panel,:),'Position', [nan nan]);  
end
textoffset = handles.point(panel).MarkerSize/2;

if not(isempty(SET(no).Point)) 
  
  if DATA.Pref.BackgroundColor
    bgcolor = 'k';
  else
    bgcolor = 'none';
  end
  
  if ~any(strcmp(DATA.ViewPanelsType{panel},{'trans3DP','sag3DP','speedim','cor3DP'}))
    %Normal case (not 3D view)

    slices = viewfunctions('slicesinpanel',panel);
    scale = viewfunctions('getscale',panel);
    
    %this is the montageindex of the
    pointstodo = find(ismember(SET(no).Point.Z,slices)&...
      (ismember(SET(no).Point.T,SET(no).CurrentTimeFrame)|isnan(SET(no).Point.T))); %nan means constant over time
    
    [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},SET(no).Point.Z(pointstodo)-slices(1)+1);
    imdim = zoomfunctions.getxysize(no,panel);
    xt = (xl-1)*imdim.XSize - imdim.XStart +1;
    yt = (yl-1)*imdim.YSize - imdim.YStart +1;

    origxdata = SET(no).Point.Y(pointstodo);
    origxdata(origxdata > imdim.YEnd) = nan;
    origxdata(origxdata < imdim.YStart) = nan;

    origydata = SET(no).Point.X(pointstodo);
    origydata(origydata > imdim.XEnd) = nan;
    origydata(origydata < imdim.XStart) = nan;

    y = scale*(origydata+xt);
    x = scale*(origxdata+yt);
  else
    return
  end
  
  %Set points
  set(handles.point(panel), 'XData', x, 'YData', y);
  
  %Set point text
  numpointstodo = numel(pointstodo);
  if ~DATA.Run && ~(numpointstodo > 20) %display point text only for the first 20 points
    for i = 1:numpointstodo
      set(handles.pointtext(panel,i),...
        'Position', [x(i)+textoffset,y(i)],...
        'String', SET(no).Point.Label(pointstodo(i)),...
        'BackgroundColor', bgcolor);
    end
  end
  
end

%---------------------
function drawinterp(panel,type) 
%-----------------------
%Draws interpolation points in panel.
global DATA SET

if nargin == 1
  type = {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'};
else
  type = {type};
end

%Drawing parameters
no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);
slicestoinclude = viewfunctions('slicesinpanel',panel);
ctf = SET(no).CurrentTimeFrame;

%Check object ind in case of General Pen
currentobjectind = 1;
if strcmp(type,'GeneralPenInterp')
  currentobjectind = DATA.GeneralPenSettings.getcurrentobject;
  markercolor = SET(no).GeneralPenObjects(currentobjectind).getcolor;
else
  markercolor = 'none';
end

%when playing we dont want to see interpolation points
if DATA.Run == 1
  slicestoinclude = [];
end
imdim = zoomfunctions.getxysize(no,panel);
for loop = 1:length(type)
  x = nan;
  y = nan;
  typeInterpX = helperfunctions('parsesetfield',SET(no),type{loop},'X',currentobjectind); %get contour
  isinterpongoing = helperfunctions('isinterpongoing',SET(no),type{loop},currentobjectind);
  if isinterpongoing && ~strcmp(DATA.Handles.([lower(type{loop}(1:end-6)),'contour'])(panel).LineStyle, 'none')
    updateinterpolationsettings(panel,type{loop})
  end
  for slice = 1:length(slicestoinclude)
    if ~isempty(typeInterpX) && ~isempty(typeInterpX{ctf,slicestoinclude(slice)})
      typeInterpX = helperfunctions('parsesetfield',SET(no),type{loop},'X',currentobjectind);%,ctf,slicestoinclude(slice)); %get contour
      typeInterpY = helperfunctions('parsesetfield',SET(no),type{loop},'Y',currentobjectind);%,ctf,slicestoinclude(slice)); %get contour
      [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},slice);
      yl = (yl - 1) * imdim.XSize - imdim.XStart +1;
      xl = (xl - 1) * imdim.YSize - imdim.YStart +1;

      origx = typeInterpY{ctf,slicestoinclude(slice)};      
      origx(origx > imdim.YEnd) = nan;
      origx(origx < imdim.YStart) = nan;

      origy = typeInterpX{ctf,slicestoinclude(slice)};
      origy(origy > imdim.XEnd) = nan;
      origy(origy < imdim.XStart) = nan;

      x = cat(1,x,nan,xl + origx);
      y = cat(1,y,nan,yl + origy);
    end
  end
  DATA.Handles.(lower(type{loop}))(currentobjectind,panel).XData = scale*x;
  DATA.Handles.(lower(type{loop}))(currentobjectind,panel).YData = scale*y;
  DATA.Handles.(lower(type{loop}))(currentobjectind,panel).MarkerFaceColor = markercolor;

  if strcmp(type{loop}(1:end-6),'LA')
    %for LA in ED and ES: show AV plane and axis line
    drawlaparameters(no,ctf,panel);
  end
end

%---------------------
function drawcontours(panel,type) 
%-----------------------
%This function draws contours in the designated panel of the input type. It
%depends on that the naming of the coordinate fields of the new contour is homogenous to the
%classical EndoX EndoY representation.
global DATA SET

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

%To reduce code amount we use dynamic calling
if nargin == 1
  type = {'Endo','Epi','RVEndo','RVEpi'};
else
  type = {type};
end

slicestoinclude = viewfunctions('slicesinpanel',panel);
numslices = length(slicestoinclude);
tf = SET(no).CurrentTimeFrame;
imdim = zoomfunctions.getxysize(no,panel);

for i = 1:length(type)
  x = [];
  y = [];
  contourtypeX = [type{i},'X'];
  contourtypeY = [type{i},'Y']; 
  if ~isempty(SET(no).(contourtypeX))
    contourlength = size(SET(no).(contourtypeY), 1);
    % preallocate x and y arrays with NaNs
    x = nan((contourlength + 1) * numslices, 1);
    y = x;
    indstart = 1;
    indend = contourlength;
    
    for sl = 1:numslices
      [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},sl);
      yl = (yl-1)*imdim.XSize - imdim.XStart +1;
      xl = (xl-1)*imdim.YSize - imdim.YStart +1;
      margin = 0.5;
      origx = SET(no).(contourtypeY)(:,tf,slicestoinclude(sl));      
      origx(origx > imdim.YEnd+margin) = nan;
      origx(origx < imdim.YStart-margin) = nan;

      origy = SET(no).(contourtypeX)(:,tf,slicestoinclude(sl));
      origy(origy > imdim.XEnd+margin) = nan;
      origy(origy < imdim.XStart-margin) = nan;

      x(indstart:indend) = xl+origx;
      y(indstart:indend) = yl+origy;
      % update start and end indexes
      indstart = indstart + contourlength +1;
      indend = indend + contourlength +1;      
    end
  end

  % check if interpolation is ongoing in current NO/stack and adjust
  % contours appearance
  if ~SET(no).([type{i},'InterpOngoing'])
    tools('resetinterpolationcontours',panel,type{i})
  else
    updateinterpolationsettings(panel,type{i})
  end

  DATA.Handles.([lower(type{i}),'contour'])(panel).XData = scale*x;
  DATA.Handles.([lower(type{i}),'contour'])(panel).YData = scale*y;
end

%---------------------
function drawcontourslara(panel,type)
%-----------------------
%This function draws atrium pen contours in the designated panel
global DATA SET

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

if nargin == 1
  type = {'LA','RA'};
else
  type = {type};
end

slicestoinclude = viewfunctions('slicesinpanel',panel);
tf = SET(no).CurrentTimeFrame;
imdim = zoomfunctions.getxysize(no,panel);
margin = 0.5;
for loop = 1:length(type)
  if ~isempty(SET(no).(type{loop})) && ~isempty(SET(no).(type{loop}).X)
    x = [];
    y = [];
    for slice = 1:length(slicestoinclude)
      [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},slice);
      yl = (yl-1)*imdim.XSize - imdim.XStart +1;
      xl = (xl-1)*imdim.YSize - imdim.YStart +1;

      origx = SET(no).(type{loop}).Y(:,tf,slicestoinclude(slice));      
      origx(origx > imdim.YEnd+margin) = nan;
      origx(origx < imdim.YStart-margin) = nan;

      origy = SET(no).(type{loop}).X(:,tf,slicestoinclude(slice));
      origy(origy > imdim.XEnd+margin) = nan;
      origy(origy < imdim.XStart-margin) = nan;

      x = cat(1,x,nan,xl+origx);
      y = cat(1,y,nan,yl+origy);
    end

    % check if interpolation is ongoing in current NO/stack and adjust
    % contours appearance
    if ~SET(no).(type{loop}).InterpOngoing
      tools('resetinterpolationcontours',panel,type{loop})
    else
      updateinterpolationsettings(panel,type{loop})
    end

    linewidth = DATA.Pref.LineWidth;
    if isempty(linewidth)
      linewidth = 1;
    end

    DATA.Handles.([lower(type{loop}),'contour'])(panel).XData = scale*x;
    DATA.Handles.([lower(type{loop}),'contour'])(panel).YData = scale*y;
    DATA.Handles.([lower(type{loop}),'contour'])(panel).LineWidth = linewidth;
  end

  if strcmp(type{loop},'LA')
    %for LA in ED and ES: show AV plane and axis line
    drawlaparameters(no,tf,panel);
  end

   if strcmp(type{loop},'RA')
    %for RV in ED and ES: show AV plane
    drawraparameters(no,tf,panel);
  end
end

%-------------------------------------------------------------
function drawlaparameters(no,tf,panel)
%-------------------------------------------------------------
%Function to decide if LA parameters should be displayed or not
global DATA SET

%Check contour visibility status
stateandicon = viewfunctions('iconson','hidela');

if ismember(tf,[SET(no).EST, SET(no).EDT]) && ~DATA.Pref.RunFDAVersion && ismember(SET(no).ImageViewPlane,{'2CH','4CH'})
  ok = drawatrialparameters('LA',SET(no).ImageViewPlane,no);
  if ok && ~stateandicon{1}
    visiblestatus = 'on';
  else
    visiblestatus = 'off';
  end
else
  visiblestatus = 'off';
end
set([DATA.Handles.laavplane(panel), DATA.Handles.laaxis(panel)],'Visible',visiblestatus);

%-------------------------------------------------------------
function varargout = drawatrialparameters(heartpart,chamberview,no,tf)
%-------------------------------------------------------------
%Function to draw AV plane and left atrial axis
arguments
  heartpart {mustBeMember(heartpart,{'LA','RA'})}
  %chamberview {mustBeMember(chamberview,{'2CH','4CH'})}
  chamberview = []
  no = []
  tf {mustBeMember(tf,{'ed','es',''})} = ''
end
global DATA SET

ok = true;

%only available for research
if DATA.Pref.RunFDAVersion
  varargout{1} = false;
  return
end

if isempty(no)
  no = findfunctions('findlaxnowithheartpart',heartpart);
  nos = findfunctions('findnoXch',chamberview);
  no = no(ismember(no,nos));
  if isempty(no)
    ok = false;
    varargout{1} = ok;
    return
  end
end

%get time frame indice
timeframe = helperfunctions('gettimeframe',tf,no);

%find slice with largest delineated area

slicetouse = findfunctions('findslicewithlargestaera',heartpart,no,timeframe);

%other parameters
panel = find(DATA.ViewPanels == no,1);
slicestoinclude = viewfunctions('slicesinpanel',panel);

if ~ismember(slicetouse,slicestoinclude)
  ok = false;
  varargout{1} = ok;
  return
end

switch DATA.ViewPanelsType{panel}
  case 'one'
    %extract X and Y coordinates
    x = SET(no).(heartpart).X(:,timeframe,slicetouse);
    y = SET(no).(heartpart).Y(:,timeframe,slicetouse);
    indstart = 1;
    indend = length(x);

    if  strcmp(heartpart, 'LA')
      %get LA parameters
      [~,startcoord,endcoord] = calcfunctions('calculateatriallength',x,y,SET(no).ResolutionX,SET(no).ResolutionY);
    else
      % [area] = calcfunctions('calculateraarea',x,y,SET(no).ResolutionX,SET(no).ResolutionY); 
    end
  case {'montage','montagerow'}
    ok = false; %not implemented, see commented code below
    % imdim = zoomfunctions.getxysize(no,panel);
    % margin = 0.5;
    %  x = [];
    %  y = [];
    %  for slice = 1:length(slicestoinclude)
    %    if slicestoinclude(slice) == slicetouse
    %      indstart = length(y) + 2;
    %    end
    %    [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},slice);
    %    yl = (yl-1)*imdim.XSize - imdim.XStart +1;
    %    xl = (xl-1)*imdim.YSize - imdim.YStart +1;
    %
    %    origx = SET(no).(heartpart).Y(:,timeframe,slicestoinclude(slice));
    %    origx(origx > imdim.YEnd+margin) = nan;
    %    origx(origx < imdim.YStart-margin) = nan;
    %
    %    origy = SET(no).(heartpart).X(:,timeframe,slicestoinclude(slice));
    %    origy(origy > imdim.XEnd+margin) = nan;
    %    origy(origy < imdim.XStart-margin) = nan;
    %
    %    y = cat(1,x,nan,xl+origx);
    %    x = cat(1,y,nan,yl+origy);
    %
    %    if slicestoinclude(slice) == slicetouse
    %      indend = length(y);
    %    end
    %  end
    % %get LA parameters
    % [~,startcoord,endcoord] = calcfunctions('calculateatriallength',x(indstart:indend),y(indstart:indend),SET(no).ResolutionX,SET(no).ResolutionY);
end

%draw parameters
if ~isempty(panel) && ok
  scale = viewfunctions('getscale',panel);
  if strcmp(heartpart, 'LA')
  set(DATA.Handles.laavplane(panel),'XData',scale*[y(indstart) y(indend)],'YData',scale*[x(indstart) x(indend)]);
  set(DATA.Handles.laaxis(panel),'XData',scale*[endcoord.y startcoord.y],'YData',scale*[endcoord.x startcoord.x]);
  else
    set(DATA.Handles.raavplane(panel),'XData',scale*[y(indstart) y(indend)],'YData',scale*[x(indstart) x(indend)]);
  end
end

if nargout > 0
  varargout{1} = ok;
end


%-------------------------------------------------------------
function drawraparameters(no,tf,panel)
%-------------------------------------------------------------
%Function to decide if RA parameters should be displayed or not
global DATA SET

%Check contour visibility status
stateandicon = viewfunctions('iconson','hidela');

if ismember(tf,[SET(no).EST, SET(no).EDT]) && ~DATA.Pref.RunFDAVersion && ismember(SET(no).ImageViewPlane,{'2CH','4CH'})
  ok = drawatrialparameters('RA',SET(no).ImageViewPlane,no);
  if ok && ~stateandicon{1}
    visiblestatus = 'on';
  else
    visiblestatus = 'off';
  end
else
  visiblestatus = 'off';
end
set(DATA.Handles.raavplane(panel),'Visible',visiblestatus);

%---------------------
function drawcontoursgeneralpen(panel) 
%-----------------------
%This function draws general pen contours in the designated panel
global DATA SET

if DATA.GeneralPenSettings.getnumobjects == 0
  return
end

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);
slicestoinclude = viewfunctions('slicesinpanel',panel);
tf = SET(no).CurrentTimeFrame;
imdim = zoomfunctions.getxysize(no,panel);
for ind = 1:DATA.GeneralPenSettings.getnumobjects(no)
  currentobj = SET(no).GeneralPenObjects(ind);
  if ~isempty(currentobj.X)
    x = [];
    y = [];
    for slice = 1:length(slicestoinclude)
      [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},slice);
      yl = (yl-1)*imdim.XSize - imdim.XStart +1;
      xl = (xl-1)*imdim.YSize - imdim.YStart +1;

      origx = currentobj.Y(:,tf,slicestoinclude(slice));      
      origx(origx > imdim.YEnd+margin) = nan;
      origx(origx < imdim.YStart-margin) = nan;

      origy = currentobj.X(:,tf,slicestoinclude(slice));
      origy(origy > imdim.XEnd+margin) = nan;
      origy(origy < imdim.XStart-margin) = nan;

      x = cat(1,x,nan,xl+origx);
      y = cat(1,y,nan,yl+origy);
    end
  end

  % check if interpolation is ongoing in current NO/stack and adjust
  % contours appearance
  if ~currentobj.InterpOngoing
    tools('resetinterpolationcontours',panel,'GeneralPen');
  else
    updateinterpolationsettings(panel,'GeneralPen');
  end

  DATA.Handles.generalpencontour(ind,panel).XData = scale*x;
  DATA.Handles.generalpencontour(ind,panel).YData = scale*y;
  DATA.Handles.generalpencontour(ind,panel).Color = currentobj.getcolor;

  %Update selected object
  if ~currentobj.InterpOngoing
    if currentobj.Selected
      linewidth = DATA.Pref.LineWidth + 1;
    else
      linewidth = DATA.Pref.LineWidth;
    end
    DATA.Handles.generalpencontour(ind,panel).LineWidth = linewidth;
  end
  
  %Draw label
  if DATA.Pref.BackgroundColor
    bgcolor = 'k';
  else
    bgcolor = 'none';
  end
  [ymin,ix] = min(currentobj.Y(:,tf,slicestoinclude));
  ypos = ymin-1;
  xpos = currentobj.X(ix,tf,slicestoinclude);
  DATA.Handles.generalpentext(ind,panel).Position = scale*[ypos xpos];
  DATA.Handles.generalpentext(ind,panel).String = currentobj.getlabel;
  DATA.Handles.generalpentext(ind,panel).BackgroundColor = bgcolor;
end

%---------------------
function unselectcontoursgeneralpen(panel) 
%-----------------------
%This function unselect all general pen contours in the designated panel
global DATA

for ind = 1:DATA.GeneralPenSettings.getnumobjects
  DATA.Handles.generalpencontour(panel,ind).LineWidth = DATA.Pref.LineWidth;
end

%---------------------
function clearcontoursgeneralpen(panel,objind) 
%-----------------------
%This function clears general pen contours in the designated panel
global DATA

DATA.Handles.generalpeninterp(objind,panel).XData = nan;
DATA.Handles.generalpeninterp(objind,panel).YData = nan;
DATA.Handles.generalpencontour(objind,panel).XData = nan;
DATA.Handles.generalpencontour(objind,panel).YData = nan;
DATA.Handles.generalpencontour(objind,panel).Color = DATA.GeneralPenSettings.getdefaultcolor;
DATA.Handles.generalpencontour(objind,panel).LineWidth = DATA.Pref.LineWidth;

%---------------------
function drawintersections(panel)
%-----------------------
%draw intersection lines in main GUI

global DATA SET
persistent markersize

%Get global graphics handles
handles = DATA.Handles;

if nargin < 1
  panels = 1:length(DATA.ViewPanels);
else
  panels = panel;
end
  
for panel = panels
  no = DATA.ViewPanels(panel);
  viewtype = DATA.ViewPanelsType{panel};
  if not(no==0)
    scale = viewfunctions('getscale',panel);
    t = SET(no).CurrentTimeFrame;
    [endointersectionx,endointersectiony] = ...
      calcfunctions('calcsegmentationintersections',panel,'endo',t,viewtype);
    [epiintersectionx,epiintersectiony] = ...
      calcfunctions('calcsegmentationintersections',panel,'epi',t,viewtype);
    
    [rvendointersectionx,rvendointersectiony] = ...
      calcfunctions('calcsegmentationintersections',panel,'rvendo',t,viewtype);
    [rvepiintersectionx,rvepiintersectiony] = ...
      calcfunctions('calcsegmentationintersections',panel,'rvepi',t,viewtype);
  
    if isempty(markersize) ||(markersize ~= DATA.Pref.MarkerSize)
      markersize = DATA.Pref.MarkerSize;
   
      set([handles.endocontourintersection(panel),...
        handles.epicontourintersection(panel),...
        handles.rvendocontourintersection(panel),...
        handles.rvepicontourintersection(panel)],'markersize',markersize);
    end
    
    % need to be set separately otherwise they all have the same color
    set(handles.endocontourintersection(panel), ...
        'XData',scale*endointersectiony,'YData',scale*endointersectionx);
    set(handles.epicontourintersection(panel), ...
        'XData',scale*epiintersectiony,'YData',scale*epiintersectionx);
    set(handles.rvendocontourintersection(panel), ...
        'XData',scale*rvendointersectiony,'YData',scale*rvendointersectionx);
    set(handles.rvepicontourintersection(panel), ...
        'XData',scale*rvepiintersectiony,'YData',scale*rvepiintersectionx);
  end
end

%------------------------------
function showviabilityedits(panel) 
%------------------------------
%Show viability edits on screen as a temporary overlay.
global DATA SET

no = DATA.ViewPanels(panel);
hidecell={'hidescar'};
stateandicon = viewfunctions('iconson',hidecell);
state=[stateandicon{:,1}];

if isempty(SET(no).Scar)|| state(1)
  return;
  %viability('viabilityreset_Callback');
end

if isequal(get(DATA.Handles.viabilityshowinfarctaswhitemenu,'checked'),'on')
  showaswhite = true;
else
  showaswhite = false;
end

if ~isempty(SET(no).Scar.GreyZone.map) && ...
    isequal(get(DATA.Handles.viabilityshowgrayzonemenu,'checked'),'on')
  showgreyzone = true;
else
  showgreyzone = false;
end

% %If neither then we can just safely exit.
% if not(showaswhite || showgreyzone)
%     return;
% end

tempnos=no;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

panelstodo = find(DATA.ViewPanels==no);

for panel=panelstodo
  %Ok lets draw it
  switch DATA.ViewPanelsType{panel}
    case {'one','mmodespatial','orth'}
      
      temp = SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      if showgreyzone
        temp = imresize(temp,size(SET(no).Scar.GreyZone.map(:,:,1)),'bilinear');
      end
      colmap = SET(no).Colormap;
      if isempty(colmap)
        colmap = colormap('gray');
      end
      tempuint8 = min(uint8(255),uint8(1)+calcfunctions('remapuint8',temp,no,colmap));
      
      sz=size(tempuint8);
      scarimrgb = reshape(tempuint8,[prod(sz(1:2)) 3]);
      tmp = SET(no).Scar.Result(:,:,SET(no).CurrentSlice);
      tmp = logical(tmp(:));
      
      %Show infarct
      if DATA.Pref.LineWidth>0
        if showaswhite
          scarimrgb(tmp,:) = uint8(255);
          %Greyzone
        elseif showgreyzone
          greytmp = (SET(no).Scar.GreyZone.map(:,:,SET(no).CurrentSlice) == 1);
          coretmp = (SET(no).Scar.GreyZone.map(:,:,SET(no).CurrentSlice) == 2);
          scarimrgb(greytmp(:),:) = repmat(uint8([127 127 0]),sum(greytmp(:)),1);
          scarimrgb(coretmp(:),:) = repmat(uint8([127 0 0]),sum(coretmp(:)),1);
        end
      end
      
      if not(isfield(SET(no).Scar,'NoReflow'))
        SET(no).Scar.NoReflow = SET(no).Scar.Auto;
      end
      
      
      tmp=SET(no).Scar.Manual(:,:,SET(no).CurrentSlice);
      tmp=tmp(:);
      scarimrgb((tmp==int8( 1)),2) = uint8(255);
      scarimrgb((tmp==int8(-1)),3) = uint8(255);
      tmp=SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice);
      tmp=tmp(:);
      scarimrgb(logical(tmp),1) = uint8(255);
      
      scarimrgb=reshape(scarimrgb,[sz(1:2) 3]); %EH: added 3
      
    case {'montage','montagerow','montagefit','sax3'}
      % Convert to 2D and layout
      colmap = SET(no).Colormap;
      if isempty(colmap)
        colmap = colormap('gray');
      end
      tempuint8 = calcfunctions('remapuint8',...
        viewfunctions('reshape2layout',squeeze(SET(no).IM(:,:,SET(no).CurrentTimeFrame,:)),panel),...
        no,colmap);
      
      tempuint8 = min(uint8(255),uint8(1)+tempuint8);
      sz = size(tempuint8);
      scarimrgb = reshape(tempuint8,[prod(sz(1:2)) 3]);
      
      tmp = viewfunctions('reshape2layout',SET(no).Scar.Result,panel);%used to be NO but this may cause error
      tmp = logical(tmp(:));
      
      %Show infarct
      if DATA.Pref.LineWidth>0
        if showaswhite
          scarimrgb(tmp,:) = uint8(255);
          %Greyzone
        elseif showgreyzone
          temp = imresize(SET(no).Scar.GreyZone.map,[SET(no).XSize SET(no).YSize],'bilinear');
          greytmp = viewfunctions('reshape2layout',temp,panel);
          gztmp = greytmp(:) == 1;
          coretmp = greytmp(:) == 2;
          
          scarimrgb(coretmp,:) = repmat(uint8([127 0 0]),sum(coretmp),1);
          scarimrgb(gztmp(:),:) = repmat(uint8([127 127 0]),sum(gztmp(:)),1);
        end
      end
      
      if not(isfield(SET(no).Scar,'NoReflow'))
        SET(no).Scar.NoReflow = repmat(uint8(1),size(SET(no).Scar.Auto));
      end
      
      
      tmp=viewfunctions('reshape2layout',SET(no).Scar.Manual,panel);
      tmp=tmp(:);
      scarimrgb((tmp==int8( 1)),2) = uint8(255);
      scarimrgb((tmp==int8(-1)),3) = uint8(255);
      tmp=viewfunctions('reshape2layout',SET(no).Scar.NoReflow,panel);
      tmp=tmp(:);
      scarimrgb(logical(tmp),1) = uint8(255);
      
      scarimrgb=reshape(scarimrgb,sz);
  end
  if not(DATA.Silent)
    scale = viewfunctions('getscale',panel);
    imxsz = size(scarimrgb,1)*scale;
    imysz = size(scarimrgb,2)*scale;
    im = imresize(scarimrgb,[imxsz imysz],'bilinear');
    set(DATA.Handles.imagehandles(panel),'CData',im);
  end
end

%--------------------------------
function draw3dpline(panel,type) 
%--------------------------------
%Draws stack intersections in 3dp.

global DATA SET
no = DATA.ViewPanels(panel);

%Get global graphics handles
h = DATA.Handles;

if DATA.Handles.tabholder.isiconindented('showcrossicon')
  %Find max size
  [rmax,gmax,bmax] = segment3dp.tools('xyz2rgb',SET(no).XSize,SET(no).YSize,SET(no).ZSize);
  maxsize = max([rmax gmax bmax]);

  switch type
    case 'br'
      x = [0, maxsize];
      y = [SET(no).LevelSet.View.RSlice SET(no).LevelSet.View.RSlice];
    case 'gr'
      x = [0, maxsize];
      y = [SET(no).LevelSet.View.RSlice SET(no).LevelSet.View.RSlice];
    case 'gb'
      x = [SET(no).LevelSet.View.BSlice SET(no).LevelSet.View.BSlice];
      y = [0, maxsize];
    case 'rb'
      x = [0, maxsize];
      y = [SET(no).LevelSet.View.BSlice SET(no).LevelSet.View.BSlice];
    case 'rg'
      x = [SET(no).LevelSet.View.GSlice SET(no).LevelSet.View.GSlice];
      y = [0, maxsize];
    case 'bg'
      x = [SET(no).LevelSet.View.GSlice SET(no).LevelSet.View.GSlice];
      y = [0, maxsize];
  end
else
  x = NaN;
  y = NaN;
end

set(h.([type,'line'])(panel),'XData',x,'YData',y);

%--------------------------------
function displaypoints(view,ima) 
%--------------------------------
%Draws annotation points in 3dp imageaxes ima according to view.

global SET DATA NO


%add points
[r,g,b] = segment3dp.tools('xyz2rgb',SET(NO).Point.X,SET(NO).Point.Y, SET(NO).Point.Z);
switch view
  case 'r'
    pointstoplot = round(r)==SET(NO).LevelSet.View.RSlice;
    x = g;
    y = b;
  case 'g'
    pointstoplot = round(g)==SET(NO).LevelSet.View.GSlice;
    x = b;
    y = r;
  case 'b'
    pointstoplot = round(b)==SET(NO).LevelSet.View.BSlice;
    x = g;
    y = r;
end

delete(DATA.Handles.([view,'text3dp']))%This removes old points

DATA.Handles.([view,'text3dp']) = [];

if any(pointstoplot)
  for loop = find(pointstoplot)
    DATA.Handles.([view,'text3dp'])=[DATA.Handles.([view,'text3dp']),text('parent',ima,...
      'position',[x(loop) y(loop)],'HorizontalAlignment','center',...
      'string','+','color',[1,1,1]), text('parent',ima,...
      'position',[x(loop)+5 y(loop)],'HorizontalAlignment','left',...
      'string',SET(NO).Point.Label{loop},'color',[1,1,1])];
    set(DATA.Handles.([view,'text3dp']),'ButtonDownFcn','segment3dp.tools(''delete3dppoints'')');
  end
end

%------------------------------------------------------------------
function [xo,yo] = intersectlinehelper(n,type,rslice,gslice,bslice)
%------------------------------------------------------------------
%Helper function to compute set of intersection points with line n and
%plane given by type that is member of {'r','g','b'}. Also works for
%pointgroups.

global DATA SET NO

%Get object
O = SET(NO).LevelSet.Object;

objecttype = O.gettype(n);
switch lower(objecttype)
  case 'pointgroup'
    v = O.getmovedpointsxyz(n);
  case 'line'
    v = O.getmovedverticesxyz(n); %in XYZ system
  otherwise
    error(sprintf('Unknown type for intersecthelper %s',type)) %#ok<SPERR> 
end

%Get slices in r,g,b coordinates
if nargin < 3
  rslice = SET(NO).LevelSet.View.RSlice;
  gslice = SET(NO).LevelSet.View.GSlice;
  bslice = SET(NO).LevelSet.View.BSlice;
end

%Convert to xyz coordinates
[x,y,z] = segment3dp.tools('rgb2xyz',rslice,gslice,bslice);

delta = 0.5;
switch type
  case 'x'
    ind = find(abs(v(:,1)-x) < delta);
  case 'y'
    ind = find(abs(v(:,2)-y) < delta);
  case 'z'
    ind = find(abs(v(:,3)-z) < delta);
end

%Decide on what to send back depending on orientation
[xo,yo] = sendbackhelper(v(ind,:)',lower(DATA.LevelSet.imageorientation),type);

%-------------------------------------------------------------------------
function [xo,yo,markersize] = intersecthelper(n,type,rslice,gslice,bslice)
%-------------------------------------------------------------------------
%Helper function to compute set of intersection points with surface n and
%plane given by type that is member of {'r','g','b'}. The parameters
%rslice,gslice,bslice are optional, if not it is taken from SET struct.

global DATA SET NO

%Get object
O = SET(NO).LevelSet.Object;

xo = []; %default
yo = []; %default
markersize = 5;

objtype = lower(O.gettype(n));

%No need if it is a bw
if isequal(objtype,'bw')
  return
end

%Check if it is a line
if isequal(objtype,'line') || isequal(objtype,'pointgroup')
  [xo,yo] = intersectlinehelper(n,type);
  markersize = 20;
  return
end

v = O.getmovedverticesxyz(n); %get moved vertices in XYZ coordinates
f = O.getfaces(n); 

%Get slices in r,g,b coordinates
if nargin < 3
  rslice = SET(NO).LevelSet.View.RSlice;
  gslice = SET(NO).LevelSet.View.GSlice;
  bslice = SET(NO).LevelSet.View.BSlice;
end

%Convert to xyz coordinates
[x,y,z] = segment3dp.tools('rgb2xyz',rslice,gslice,bslice);

triangles = [v(f(:,1),:) v(f(:,2),:) v(f(:,3),:)]; %convert to triangles format, i.e 9 columns

switch type
  case 'x'
    %Find min/max x of each triangle
    minx = min(triangles(:,[1 4 7]),[],2);
    maxx = max(triangles(:,[1 4 7]),[],2);
    triangleinds = ((minx <= x) & (maxx > x)); 
  case 'y'
    %Find min/max y of each triangle
    miny = min(triangles(:,[2 5 8]),[],2);
    maxy = max(triangles(:,[2 5 8]),[],2);
    triangleinds = ((miny <= y) & (maxy > y)); 
  case 'z'
    %Find min/max z of each triangle
    minz = min(triangles(:,[3 6 9]),[],2);
    maxz = max(triangles(:,[3 6 9]),[],2);
    triangleinds = ((minz <= z) & (maxz > z));
end

%Get triangles
subsettriangles = triangles(triangleinds,:);
subsettriangles = subsettriangles';

%Get points
p1 = subsettriangles(1:3,:);
p2 = subsettriangles(4:6,:);
p3 = subsettriangles(7:9,:);

switch type
  case 'x'
    c = ones(1,size(p1,2))*x;
    P = [ones(1,size(p1,2));zeros(1, size(p1,2));zeros(1,size(p1,2))];
  case 'y'
    c = ones(1,size(p1,2))*y;
    P = [zeros(1,size(p1,2));ones(1,size(p1,2));zeros(1, size(p1,2))];
  case 'z'
    c = ones(1,size(p1,2))*z;
    P = [zeros(1,size(p1,2));zeros(1, size(p1,2));ones(1,size(p1,2))];
end

t1 = (c-sum(P.*p1))./sum(P.*(p1-p2));
t2 = (c-sum(P.*p2))./sum(P.*(p2-p3));
t3 = (c-sum(P.*p3))./sum(P.*(p3-p1));
intersect1 = p1 + times(p1-p2,t1);
intersect2 = p2 + times(p2-p3,t2);
intersect3 = p3 + times(p3-p1,t3);
i1 = intersect1(3,:) <= max(p1(3,:),p2(3,:)) & intersect1(3,:) >= min(p1(3,:),p2(3,:));
i2 = intersect2(3,:) <= max(p2(3,:),p3(3,:)) & intersect2(3,:) >= min(p2(3,:),p3(3,:));
i3 = intersect3(3,:) <= max(p3(3,:),p1(3,:)) & intersect3(3,:) >= min(p3(3,:),p1(3,:));

imain = i1+i2+i3 == 2;

ptsout = [[intersect1(:,i1&i2&imain);intersect2(:,i1&i2&imain)],[intersect2(:,i2&i3&imain);intersect3(:,i2&i3&imain)], [intersect3(:,i3&i1&imain);intersect1(:,i3&i1&imain)]];

%Decide on what to send back depending on orientation
[xo,yo] = sendbackhelper(ptsout,lower(DATA.LevelSet.imageorientation),type);

%----------------------------------------------------
function [xo,yo] = sendbackhelper(p,orientation,type)
%----------------------------------------------------
%Decide on what to send back depending on orientation

xo = [];
yo = [];

switch orientation
  case 'transversal'
    switch type
      case 'x'
        xo = p(3,:);
        yo = p(2,:);
      case 'y'
        xo = p(3,:);
        yo = p(1,:);
      case 'z'
        xo = p(1,:);
        yo = p(2,:);
    end
  case 'sagittal'
    switch type
      case 'x'
        xo = p(2,:);
        yo = p(3,:);
      case 'y'
        xo = p(1,:);
        yo = p(3,:);
      case 'z'
        xo = p(1,:);
        yo = p(2,:);
    end
  case 'coronal'
    switch type
      case 'x'
        xo = p(3,:);
        yo = p(2,:);
      case 'y'
        xo = p(1,:);
        yo = p(3,:);
      case 'z'
        xo = p(1,:);
        yo = p(2,:);
    end
end

%---------------------------------------------
function draw3dpsurfaceoutline(panel,varargin)
%---------------------------------------------
%Draws outline of visible surface objects. 

global DATA SET NO

if DATA.LevelSet.motionon  
  return
end

%Get object
O = SET(NO).LevelSet.Object;

maxsurfaces = length(DATA.Handles.surfaceintersectline{panel});

visibleobjects = O.getvisibleobjects({'surface','pin','line','cutplane','pointgroup'});

xmerge = [];
ymerge = [];

for loop = 1:length(visibleobjects)

  %get object
  n = visibleobjects(loop);
  paneltype = DATA.ViewPanelsType{panel};
  coord = '';
  switch lower(DATA.LevelSet.imageorientation)
    case 'transversal'
      switch paneltype
        case 'trans3DP'
          coord = 'z';
        case 'sag3DP'
          coord = 'y';
        case 'cor3DP'
          coord = 'x';
      end
    case 'sagittal'
      switch paneltype
        case 'trans3DP'
          coord = 'x';
        case 'sag3DP'
          coord = 'z';
        case 'cor3DP'
          coord = 'y';
      end
    case 'coronal'
      switch paneltype
        case 'trans3DP'
          coord = 'x';
        case 'sag3DP'
          coord = 'y';
        case 'cor3DP'
          coord = 'z';
      end
  end

  %Call intersect helper
  if ~isempty(coord)
    [x,y,markersize] = intersecthelper(n,coord);
  else
    x = [];
    y = [];
    markersize = 5;
  end

  if loop <= (maxsurfaces-1)
    %Normal when there are free handles
    set(DATA.Handles.surfaceintersectline{panel}(loop),'XData',y,'YData',x,'Color',O.getcolor(n)/255,'Marker','.','MarkerSize',markersize,'Visible','on');
  else
    %Merge them to one
    xmerge = [xmerge NaN x]; %#ok<AGROW> 
    ymerge = [ymerge NaN y]; %#ok<AGROW> 
  end

end %loop over visible objects

if ~isempty(xmerge)
  set(DATA.Handles.surfaceintersectline{panel}(maxsurfaces),'XData',y,'YData',x,'Color',[0.3 0.3 0.3],'Visible','on');  
end

%Loop over lines that are not used
for loop = (length(visibleobjects)+1):(length(DATA.Handles.surfaceintersectline{panel}))
  set(DATA.Handles.surfaceintersectline{panel}(loop),'Visible','off');
end

%--------------------------------------
function draw3dpoutline(panel,varargin) 
%--------------------------------------
%Shows contour of 3d segmentation if outline button is indented and motion function isnt running.

global DATA SET NO

if DATA.LevelSet.motionon
  return
end

%Get object
O = SET(NO).LevelSet.Object;
visibleobjects = O.getvisibleobjects('BW'); %Get all

im = [];

%Get image of first visible object
if ~isempty(visibleobjects)
  bw = O.getbw(visibleobjects(1));
  switch DATA.ViewPanelsType{panel}
    case 'trans3DP'
      im = double(getimagehelper(bw,'r'));
    case 'sag3DP'
      im = double(getimagehelper(bw,'g'));
    case 'cor3DP'
      im = double(getimagehelper(bw,'b'));
    case 'speedim'
      im = double(getimagehelper(bw,SET(NO).LevelSet.Pen.Color));
  end
end

%Loop over the rest of visible objects
for loop = 2:length(visibleobjects)
  bw = O.getbw(visibleobjects(loop));

  switch DATA.ViewPanelsType{panel}
    case 'trans3DP'
      im = max(im,double(getimagehelper(bw,'r')));
    case 'sag3DP'
      im = max(im,double(getimagehelper(bw,'g')));
    case 'cor3DP'
      im = max(im,double(getimagehelper(bw,'b')));
    case 'speedim'
      im = max(im,double(getimagehelper(bw,SET(NO).LevelSet.Pen.Color)));
  end
end

%no segmentation
if sum(im(:)>127)==0 %all(im==0)
  DATA.Handles.threedpcontour(panel).ZData = [127 0 ; 0 0];
  return
end

DATA.Handles.threedpcontour(panel).ZData = im;
DATA.Handles.threedpcontour(panel).Visible = 'on'; 

%-----------------------------------------
function im = getimagehelper(vol,view)
%----------------------------------------
%Similar as getimage, but do assme no TSize as stored locally in BW.

global DATA SET NO  

sz = size(vol);
sz1 = sz(1);
sz2 = sz(2);
sz3 = size(vol,3);

switch lower(DATA.LevelSet.imageorientation)
  case 'transversal'
    switch view
      case 'r'
        im = squeeze(vol(:,:,max(1,min(sz3,SET(NO).LevelSet.View.RSlice))));
      case 'g'
        im = squeeze(vol(:,max(1,min(sz2,SET(NO).LevelSet.View.GSlice)),:))';
      case 'b'
        im = squeeze(vol(max(1,min(sz1,SET(NO).LevelSet.View.BSlice)),:,:))';
    end
  case 'sagittal'
    switch view
      case 'r'
        im = squeeze(vol(max(1,min(sz1,SET(NO).LevelSet.View.RSlice)),:,:));
      case 'g'
        im = squeeze(vol(:,:,max(1,min(sz3,SET(NO).LevelSet.View.GSlice))));
      case 'b'
        im = squeeze(vol(:,max(1,min(sz2,SET(NO).LevelSet.View.BSlice)),:));
    end
  case 'coronal'
    switch view
      case 'r'
        im = squeeze(vol(max(1,min(sz1,SET(NO).LevelSet.View.RSlice)),:,:))';
      case 'g'
        im = squeeze(vol(:,max(1,min(sz2,SET(NO).LevelSet.View.GSlice)),:));
      case 'b'
        im = squeeze(vol(:,:,max(1,min(sz3,SET(NO).LevelSet.View.BSlice))));
    end
end

%---------------------------
function im = getimage(view)
%---------------------------
%Get image for view given imageorientation, r = transversal, g = sagittal,
%b=coronal

global DATA SET NO

tf = 1;
if SET(NO).TSize > 0
  tf = SET(NO).CurrentTimeFrame;
end

sz1 = SET(NO).XSize;
sz2 = SET(NO).YSize;
sz3 = SET(NO).ZSize;

orient = lower(DATA.LevelSet.imageorientation);

switch orient
  case 'transversal'
    switch view
      case 'r'
        im = SET(NO).IM(:,:,tf,max(1,min(sz3,SET(NO).LevelSet.View.RSlice)));
      case 'g'
        im = SET(NO).IM(:,max(1,min(sz2,SET(NO).LevelSet.View.GSlice)),tf,:);
      case 'b'
        im = SET(NO).IM(max(1,min(sz1,SET(NO).LevelSet.View.BSlice)),:,tf,:);
    end
  case 'sagittal'
    switch view
      case 'r'
        im = SET(NO).IM(max(1,min(sz1,SET(NO).LevelSet.View.RSlice)),:,tf,:);
      case 'g'
        im = SET(NO).IM(:,:,tf,max(1,min(sz3,SET(NO).LevelSet.View.GSlice)));
      case 'b'
        im = SET(NO).IM(:,max(1,min(sz2,SET(NO).LevelSet.View.BSlice)),tf,:);
    end
  case 'coronal'
    switch view
      case 'r'
        im = SET(NO).IM(max(1,min(sz1,SET(NO).LevelSet.View.RSlice)),:,tf,:);
      case 'g'
        im = SET(NO).IM(:,max(1,min(sz2,SET(NO).LevelSet.View.GSlice)),tf,:);
      case 'b'
        im = SET(NO).IM(:,:,tf,max(1,min(sz3,SET(NO).LevelSet.View.BSlice)));
    end
end

if ~ismatrix(im)
  im = squeeze(im);
  if (strcmp(orient,'transversal') && (strcmp(view,'g') || strcmp(view,'b'))) ...
      || (strcmp(orient,'coronal') && strcmp(view,'r'))
    im = im';
  end
end

%----------------------------------------
function drawintersectionpoints(no,panel)
%----------------------------------------
%Initiate handles and draw intersections with other contours
global DATA

if isempty(DATA.Handles.hideothercontouricon) || isequal(get(DATA.Handles.hideothercontouricon,'state'),'off')
  [endointersectline,maxintersect] = segment('getendointersection',no);
  viewtype = DATA.ViewPanelsType{panel};
  [endox,endoy] = calcfunctions('calcsegmentationintersections',no,'endo',viewtype);
  [epix,epiy] = calcfunctions('calcsegmentationintersections',no,'epi',viewtype);
  [rvendox,rvendoy] = calcfunctions('calcsegmentationintersections',no,'rvendo',viewtype);
  [rvepix,rvepiy] = calcfunctions('calcsegmentationintersections',no,'rvepi',viewtype);
  
  DATA.Handles.endointersectionline{panel} = zeros(1,maxintersect);
  for i=1:maxintersect
    if(i<=length(endointersectline))
      n = endointersectline(i).NPoints;
      DATA.Handles.endointersectionline{panel}(i) = plot(DATA.Handles.imageaxes(panel),...
        endointersectline(i).Y(1:n),endointersectline(i).X(1:n),'r');
    else
      % create empty handle for intersections in other time frames
      DATA.Handles.endointersectionline{panel}(i) = plot(DATA.Handles.imageaxes(panel),...
        [0 1],[0 1],'r');
      set(DATA.Handles.endointersectionline{panel}(i),'visible','off');
    end
  end
  DATA.Handles.endointersectionpoints{panel} = plot(DATA.Handles.imageaxes(panel),...
    endoy,endox,'r.');
  DATA.Handles.epiintersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    epiy,epix,'g.');
  DATA.Handles.rvendointersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    rvendoy,rvendox,'m.');
  DATA.Handles.rvepiintersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    rvepiy,rvepix,'c.');
else
  for i = 1:length(DATA.Handles.endointersectionline{panel})
    DATA.Handles.endointersectionline{panel}(i) = plot(DATA.Handles.imageaxes(panel),...
      NaN,NaN,'r');
  end
  DATA.Handles.endointersectionpoints{panel} = plot(DATA.Handles.imageaxes(panel),...
    NaN,NaN,'r.');
  DATA.Handles.epiintersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    NaN,NaN,'g.');
  DATA.Handles.rvendointersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    NaN,NaN,'m.');
  DATA.Handles.rvepiintersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    NaN,NaN,'c.');
end

if ~DATA.Pref.BlackWhite
  set(DATA.Handles.endointersectionpoints{panel},'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.epiintersection{panel},'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.rvendointersection{panel},'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.rvepiintersection{panel},'markersize',DATA.Pref.MarkerSize);
else
  set(DATA.Handles.endointersectionline{panel},'color',[1 1 1]);
  set(DATA.Handles.endointersectionpoints{panel},'markersize',DATA.Pref.MarkerSize,'color',[1 1 1]);
  set(DATA.Handles.epiintersection{panel},'markersize',DATA.Pref.MarkerSize,'color',[1 1 1]);
  set(DATA.Handles.rvendointersection{panel},'markersize',DATA.Pref.MarkerSize,'color',[1 1 1]);
  set(DATA.Handles.rvepiintersection{panel},'markersize',DATA.Pref.MarkerSize,'color',[1 1 1]);
end

%---------------------------------------
function updateintersectionpoints(panel)
%---------------------------------------
%Draw intersection with segmentation in other image stacks

global DATA
no = DATA.ViewPanels(panel(1));
viewtype = DATA.ViewPanelsType{panel};
if ~isequal(get(DATA.Handles.hideothercontouricon,'state'),'on')
  endointersect = segment('getendointersection',no);
  [endox,endoy] = calcfunctions('calcsegmentationintersections',no,'endo',viewtype);
  [epix,epiy] = calcfunctions('calcsegmentationintersections',no,'epi',viewtype);
  [rvendox,rvendoy] = calcfunctions('calcsegmentationintersections',no,'rvendo',viewtype);
  [rvepix,rvepiy] = calcfunctions('calcsegmentationintersections',no,'rvepi',viewtype);
  set(DATA.Handles.endointersectionline{panel},'visible','off');
  for i=1:length(endointersect)
    n = endointersect(i).NPoints;
    set(cellref(DATA.Handles.endointersectionline(panel),i),'xdata',endointersect(i).Y(1:n),...
      'ydata',endointersect(i).X(1:n),'visible','on');
  end
  if ~isempty([endox epix rvendox rvepix])
    set([DATA.Handles.endointersectionpoints{panel}], ...
      'XData',endoy,'YData',endox);
    set([DATA.Handles.epiintersection{panel}], ...
      'XData',epiy,'YData',epix);
    set([DATA.Handles.rvendointersection{panel}], ...
      'XData',rvendoy,'YData',rvendox);
    set([DATA.Handles.rvepiintersection{panel}], ...
      'XData',rvepiy,'YData',rvepix);
  end
end
%---------------------------------------
function updateinterpolationsettings(panel,type)
%---------------------------------------
%update the interpolation settings for contour display
global DATA

currentobjectind = 1;
if contains(type,'GeneralPen')
  currentobjectind = DATA.GeneralPenSettings.getcurrentobject;
end

DATA.Handles.([lower(type),'contour'])(currentobjectind,panel).LineStyle = 'none';
DATA.Handles.([lower(type),'contour'])(currentobjectind,panel).Marker = '.';
DATA.Handles.([lower(type),'contour'])(currentobjectind,panel).MarkerSize = 4;

%---------------------------------------
function updatexydata(handlename,panel,szim)
%---------------------------------------
% this function updates XData and YData depending on the scale and image
% size
global DATA

%---Get global graphics handles
handles = DATA.Handles;

if DATA.Pref.ViewInterpolated && strcmp(DATA.ViewPanelsType{panel},'one')
  %Assume scale = 2
  set(handles.(handlename)(panel),'XData',1.5:1:szim(2)+0.5,... %XData is corner coordinates for image
    'YData',1.5:1:szim(1)+0.5);
else
  set(handles.(handlename)(panel),'XData',1:1:szim(2),... %XData is corner coordinates for image
    'YData',1:1:szim(1));
end

%------------------------------------
function updateblackandwhite_Callback
%------------------------------------
%Sets colors of contours to match black and white settings

global DATA

%Extract handles
h = DATA.Handles;

if DATA.Pref.BlackWhite
  %I.e set them to white
   hvec = [...
    h.generalpencontour ...
    h.lacontour ...
    h.laavplane ...
    h.laaxis ...
    h.racontour ...
    h.raavplane ...
    h.endocontour ...
    h.epicontour ...
    h.rvendocontour ...
    h.rvepicontour ...
    h.scarcontour ...
    h.weightedscarcontour ...
    h.mocontour ...
    h.endointerp ...
    h.epiinterp ...
    h.rvendointerp ...
    h.rvepiinterp ...
    h.generalpeninterp ...
    h.lainterp ...
    h.rainterp ...
    h.moextentcontour ...
    h.endocontourintersection ...
    h.epicontourintersection ...
    h.rvendocontourintersection ...
    h.rvepicontourintersection];
   set(hvec,'Color',[1 1 1]);

    %Ignore the ones that are already white
    %h.marcontour ...
    %h.text ...
    %h.planeintersection];

else
  %Divide on colors

  %Red
  hvec = [...
    h.endocontour ...
    h.endointerp ...
    h.endocontourintersection ...
    h.mocontour ...
    h.moextentcontour];
  set(hvec,'Color',[1 0 0]);

  %Green
  hvec = [...
    h.epicontour ...
    h.epicontourintersection ...
    h.epiinterp];
  set(hvec,'Color',[0 1 0]);

  %Magenta
  hvec = [...
    h.rvendocontour ...
    h.rvendointerp ...
    h.rvendocontourintersection];
  set(hvec,'Color','m');

  %Cyan
  hvec = [...
      h.rvepicontour ...
      h.rvepiinterp ...
      h.rvepicontourintersection];
  set(hvec,'Color','c');

  %Yellow
  set(h.scarcontour,'Color','y');

  %Pink
  set(h.weightedscarcontour,'Color',[1 0.5 0.5]);

  %White => ignore
  %hvec = [...
  %  h.marcontour ...
  %  h.text ...
  %  h.planeintersection];

end