 function drawfunctions(varargin)
% DRAWFUNCTIONS
% Functions for drawing image panels, slices and montages
%
% Moved out from segment_main by Nisse Lundahl

%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
feval(varargin{:}); % FEVAL switchyard

%----------------------------------------
function drawimageview(nos,sz,panelstype)
%----------------------------------------
%Draws a full view specified by the arguments
global DATA NO

if nargin < 3
    panelstype = cell(sz);
    [panelstype{:}] = deal('one');
end
if DATA.CurrentPanel > numel(nos) || nos(DATA.CurrentPanel) ~= NO
    DATA.CurrentPanel = find(nos==NO,1,'last');
    if isempty(DATA.CurrentPanel)
        DATA.CurrentPanel = find(nos,1,'last');
        NO = nos(DATA.CurrentPanel);
    end
end

DATA.ViewPanelsType = panelstype;
for i = 1:numel(nos)
    if nos(i) > 0 && ismember(panelstype{i},{'montage'})
        [rows,cols] = calcfunctions('calcrowscols',nos(i));
        DATA.ViewPanelsMatrix{i} = [rows cols];
    else
        DATA.ViewPanelsMatrix{i} = [];
    end
end
DATA.ViewPanels = nos;
DATA.ViewIM = cell(1,numel(nos));

drawall(sz(1),sz(2));

%--------------------
function drawall(n,m)
%--------------------
%Draws all panels that should be visible or up to the number specified as n.
%  If called with n=[], then lookup from imagepanels
%  When n,m specified draw n*m panels.

global DATA SET NO

gotrowscols = false;

if DATA.Silent
    return;
end;

if nargin==0
    pos = find(DATA.ViewPanels==0);
    if isempty(pos)
        n = length(DATA.ViewPanels);
    else
        n = pos(1)-1;
    end;
end;

if isempty(n)
    n = length(DATA.ViewPanels);
end;

if length(n)>1
    temp = n;
    n = n(1);
    m = temp(2);
    n = n*m; %see below
    gotrowscols = true;
end;

%If called by two input arguments then it is row*cols
if nargin==2
    n = n*m;
    gotrowscols = true;
end;

% set([...
%     DATA.Handles.view1panelicon ...
%     DATA.Handles.view2panelicon ...
%     DATA.Handles.view2x1panelicon ...
%     DATA.Handles.view3panelicon ...
%     DATA.Handles.view1x3panelicon ...
%     DATA.Handles.view4panelicon ...
%     DATA.Handles.orthoviewicon ...
%     DATA.Handles.mipicon ...
%     DATA.Handles.view6panelicon ...
%     DATA.Handles.view9panelicon ...
%     DATA.Handles.view12panelicon ...
%     DATA.Handles.view16panelicon],'state','off');

%Delete the old ones first
try
    delete(DATA.Handles.imageaxes);
    delete(DATA.Handles.boxaxes);
catch %#ok<CTCH>
end;

%--- Find out configuration
if ~gotrowscols
    switch n
        case {0,1}
            rows = 1;
            cols = 1;
        case 2
            rows = 1;
            cols = 2;
        case 3
            rows = 3;
            cols = 1;
        case 4
            rows = 2;
            cols = 2;
        case {5,6}
            rows = 2;
            cols = 3;
        case {7,8}
            rows = 2;
            cols = 4;
        case 9
            rows = 3;
            cols = 3;
        case {10,11,12}
            rows = 3;
            cols = 4;
        case {13,14,15}
            rows = 3;
            cols = 5;
        case 16
            rows = 4;
            cols = 4;
        case {17,18}
            rows = 3;
            cols = 6;
        case {19,20}
            rows = 4;
            cols = 5;
        case 21
            rows = 3;
            cols = 7;
        case {22,23,24}
            rows = 4;
            cols = 6;
        otherwise
            cols = ceil(sqrt(n));
            rows = ceil(n/cols);
    end;
else
    %Two input arguments
    rows = n/m; %see above, nargin==2 clause before.
    cols = m;
end;

%Sets viewicon 
%DATA.setviewbuttons
% 
% %Fix icons
% if (rows==1)&&(cols==1)
%     set(DATA.Handles.view1panelicon,'state','on');
% elseif (rows==1)&&(cols==2)
%     set(DATA.Handles.view2panelicon,'state','on');
% elseif (rows==2)&&(cols==1)
%     set(DATA.Handles.view2x1panelicon,'state','on');
% elseif (rows==3)&&(cols==1)
%     set(DATA.Handles.view3panelicon,'state','on');
% elseif (rows==1)&&(cols==3)
%     set(DATA.Handles.view1x3panelicon,'state','on');
% elseif (rows==2)&&(cols==2)
%     set(DATA.Handles.view4panelicon,'state','on');
% elseif (rows==2)&&(cols==3)
%     set(DATA.Handles.view6panelicon,'state','on');
% elseif (rows==3)&&(cols==3)
%     set(DATA.Handles.view9panelicon,'state','on');
% elseif (rows==3)&&(cols==4)
%     set(DATA.Handles.view12panelicon,'state','on');
% elseif (rows==4)&&(cols==4)
%     set(DATA.Handles.view16panelicon,'state','on');
% end;
%Updates viewbuttons zero is crucial else stack overflow
%DATA.setviewbuttons(0)

%Layout has changed, reset zoom.
if ~isequal([rows cols],DATA.ViewMatrix)
    for loop=1:length(SET)
        SET(loop).NormalZoomState = [];
        SET(loop).MontageZoomState = [];
        SET(loop).MontageRowZoomState = [];
        SET(loop).MontageFitZoomState = [];
    end;
end
%Store for future use.
DATA.ViewMatrix = [rows cols];

%--- Add new axes

left = DATA.GUISettings.LeftGapWidth; %0.12;
right = 1-DATA.GUISettings.RightGapWidth-0.02;  %Based on that the report panel can never be more than 220.
bottom = DATA.GUISettings.BottomGapHeight; %0.013;
top = 1-DATA.GUISettings.TopGapHeight;
width = right-left;
height = top-bottom;

%Create image area with grid for boxes.
DATA.Handles.boxaxes = axes('position',...
    [left bottom width height],...
    'parent',DATA.imagefig);

%2 => 0.5
%3 => 0.333 0.666
h = [];
hold(DATA.Handles.boxaxes,'on');
for rloop=1:(rows-1)
    h = [h;plot(DATA.Handles.boxaxes,[0 1],[1/rows*rloop 1/rows*rloop],DATA.GUISettings.BoxAxesLineSpec)]; %#ok<AGROW>
end;
for cloop=1:(cols-1)
    h = [h;plot(DATA.Handles.boxaxes,[1/cols*cloop 1/cols*cloop],[0 1],DATA.GUISettings.BoxAxesLineSpec)]; %#ok<AGROW>
end;
h = [h;plot(DATA.Handles.boxaxes,[0 1],[0 0],DATA.GUISettings.BoxAxesLineSpec)];
h = [h;plot(DATA.Handles.boxaxes,[0 1],[1 1],DATA.GUISettings.BoxAxesLineSpec)];
h = [h;plot(DATA.Handles.boxaxes,[0 0],[0 1],DATA.GUISettings.BoxAxesLineSpec)];
h = [h;plot(DATA.Handles.boxaxes,[1 1],[0 1],DATA.GUISettings.BoxAxesLineSpec)];
hold(DATA.Handles.boxaxes,'off');
axis(DATA.Handles.boxaxes,'off');
set(h,'color',DATA.GUISettings.BoxAxesColor,'linewidth',2);

h = [];
for rloop=(rows-1):-1:0;
    for cloop=0:(cols-1);
        h = [h axes('position',...
            [left+cloop*width/cols bottom+rloop*height/rows width/cols height/rows],...
            'parent',DATA.imagefig)]; %#ok<AGROW>
    end;
end;

%Store handles
DATA.Handles.imageaxes = h;

totn = cols*rows;

%--- Hide all images
set(DATA.Handles.imageaxes,...
    'visible','off',...
    'xtick',[],'ytick',[]);

%---Add to list of what panels contains what
if (n==1)
    if ~isempty(DATA.ViewPanels)
        %If one the display current
        DATA.ViewPanels = DATA.ViewPanels(DATA.CurrentPanel);
        DATA.ViewPanelsType = DATA.ViewPanelsType(DATA.CurrentPanel);
        DATA.ViewPanelsMatrix = DATA.ViewPanelsMatrix(DATA.CurrentPanel);
        DATA.ViewIM = DATA.ViewIM(DATA.CurrentPanel);
        if ismember(DATA.ViewPanelsType{1},{'ortho','hla','vla','gla'});
            if ~strcmp(DATA.ViewPanelsType{1},'ortho')
                DATA.ViewIM{1} = [];
            end
            DATA.ViewPanelsType{1} = 'one';
        end
        if numel(DATA.Overlay) >= DATA.CurrentPanel
            DATA.Overlay = DATA.Overlay(DATA.CurrentPanel);
        end
        DATA.CurrentPanel = 1;
    else
        DATA.ViewPanels = [];
        DATA.ViewPanelsType = {};
        DATA.ViewPanelsMatrix = {};
        DATA.ViewIM = {};
        DATA.Overlay = [];
        segment('addtopanels',1);
    end;
else
    if length(DATA.ViewPanels)>totn
        DATA.ViewPanels = DATA.ViewPanels(1:totn); %remove
        DATA.ViewPanelsType = DATA.ViewPanelsType(1:totn); %remove
        DATA.ViewPanelsMatrix = DATA.ViewPanelsMatrix(1:totn); %remove
        DATA.ViewIM = DATA.ViewIM(1:totn); %remove
        if totn < 4
            [DATA.ViewIM{ismember(DATA.ViewPanelsType, ...
                {'hla','vla','gla'})}] = deal([]);
            [DATA.ViewPanelsType{ismember(DATA.ViewPanelsType, ...
                {'ortho','hla','vla','gla'})}] = deal('one');
        end
        if not(isempty(DATA.Overlay))
            DATA.Overlay = DATA.Overlay(1:totn);
        end
    else
        if length(DATA.ViewPanels)<totn
            %Case where images to show are fewer than view elements
            for loop=length(DATA.ViewPanels)+1:totn
                DATA.ViewPanels(loop) = 0; %assigning to new last element will create empty in between
                DATA.ViewPanelsType{loop} = '';
                DATA.ViewPanelsMatrix{loop} = [];
                DATA.ViewIM{loop} = [];
                if ~isempty(DATA.Overlay)
                    DATA.Overlay(loop) = struct('alphadata',[],'cdata',[]);
                end
            end
        end;
    end;
end;


%Reserve space for other handles.
%xxx
empty = NaN(1,totn);
DATA.Handles.imagehandle   = empty; %Handle to image
DATA.Handles.endocontour   = empty;
DATA.Handles.epicontour    = empty;
DATA.Handles.rvendocontour = empty; %Right ventricle contour
DATA.Handles.atrialscarcontour = empty; %Atrial scar contour with yellow dots.
DATA.Handles.rvepicontour = empty;
DATA.Handles.endopin = empty; %Pins at endocardium
DATA.Handles.epipin = empty;
DATA.Handles.rvendopin = empty;
DATA.Handles.rvepipin = empty;
DATA.Handles.cursor = empty; %White middle cursor
DATA.Handles.imagetypetext = empty;
DATA.Handles.seriesdescriptiontext = empty;
DATA.Handles.dicomimagetypetext = empty;
DATA.Handles.slicetimetext = empty;
DATA.Handles.endointerp = empty;
DATA.Handles.epiinterp = empty;
DATA.Handles.rvendointerp = empty;
DATA.Handles.rvepiinterp = empty;
DATA.Handles.overlayhandle = empty; %Handles to image overlays
DATA.Handles.endointersectionline = cell(1,totn);
DATA.Handles.endointersectionpoints = cell(1,totn);
DATA.Handles.epiintersection = cell(1,totn);
DATA.Handles.rvendointersection = cell(1,totn);
DATA.Handles.rvepiintersection = cell(1,totn);
DATA.Handles.center = cell(1,totn);
DATA.Handles.roicontour = cell(1,totn); %Roi's
DATA.Handles.roitext = cell(1,totn); %Roi's
DATA.Handles.scarcontour = cell(1,totn); %Scar contour
DATA.Handles.weightedscarcontour = cell(1,totn); %Weighted scar contour
DATA.Handles.mocontour = cell(1,totn); %microvascular obstruction contour
DATA.Handles.moextentcontour = cell(1,totn); %microvascular obstruction contour
DATA.Handles.marcontour = cell(1,totn); %myocardium at risk contour
DATA.Handles.planeintersectionline = cell(1,totn); %Plan intersection
DATA.Handles.pointo = cell(1,totn); %Points
DATA.Handles.pointp = cell(1,totn);
DATA.Handles.pointtext = cell(1,totn); %text of points
DATA.Handles.measureline = cell(1,totn); %measure line
[DATA.Handles.measureline{:}] = deal({});
DATA.Handles.measuretext = cell(1,totn); %measure text
DATA.Handles.selectslicehandle = cell(1,totn); %handle to yellow box sourrounding selected slices
DATA.Handles.sectorgrid = cell(1,totn);

DATA.Handles.phasecursor = []; %only one global sets by manualdrawbuttondown
%--- Call draw image

panelstodo = find(DATA.ViewPanels>0);

if (length(panelstodo)==1)&&(isequal(DATA.ViewPanelsType{1},'mmodespatial'))
    DATA.ViewPanelsType{1} = 'one';
end;

%Weird but necessary switcheroo.
temppanel=min(DATA.CurrentPanel,length(DATA.ViewPanels));
temppanel2=DATA.CurrentPanel;
DATA.CurrentPanel=temppanel;

for panel=panelstodo;
    no = DATA.ViewPanels(panel);
    drawimagepanel(panel);
    showedits(no);
end;

%DATA.CurrentPanel=temppanel2;
segment('switchtopanel',temppanel2,1);
%DATA.CurrentPanel=temppanel;

%Draw intersections
drawintersections;

%Enable/Disable play icon etc
DATA.updatetimethings;

% %assure DATA.LVNO is set when updating flow table
% findfunctions('setglobalstacks');

%Update volumeaxes
%segment('updatevolume');
%segment('updateflow');

%Update title on window
DATA.updatetitle;

%Update theme
%DATA.updateicons;
%Update callbacks. Needed to set buttondown on intersection line
updatetool;

%Render iconplaceholders aswell
DATA.Handles.toggleiconholder.render
DATA.Handles.permanenticonholder.render
if ~isempty(DATA.Handles.configiconholder.cdata)
  DATA.Handles.configiconholder.render
end
%Upate view icons
%segment('updateviewicons');

% %Updates viewbuttons zero is crucial else stack overflow
%DATA.setviewbuttons(0)

%Make now one orange
set(DATA.Handles.imageaxes,...
    'xcolor',[0 0 0],'ycolor',[0 0 0]);
set(DATA.Handles.imageaxes(DATA.CurrentPanel),...
    'xcolor',DATA.GUISettings.AxesColor,'ycolor',DATA.GUISettings.AxesColor,'linewidth',2.5);

%-----------------------------
function drawimagepanel(panel)
%-----------------------------
%Draws selected panels, used upon loading or when when major
%changes have occured such as change of slice, added measurement
%points etc. This fcn is the true workhorse in graphics etc.

global DATA SET

if DATA.Silent
    return;
end;

if nargin < 1
    panel = DATA.CurrentPanel;
end

if length(panel) == 1
    no = DATA.ViewPanels(panel);
else
    no = DATA.ViewPanels(panel(1));
    if ~isempty(SET(no).Parent)
        no = SET(no).Parent;
    end
end

%--- Find if to use viability
viashow = false;
if not(isempty(SET(no).Scar))
    if not(isequal(SET(no).Scar.Mode,'none'))
        viashow = true;
    end;
end;
stateandicon=segment('iconson',{'hidescar','hidemar','hideall'});
state=stateandicon{1,1};
if state%isequal(get(DATA.Handles.hidescaricon,'state'),'on')
    viashow = false;
end;

%--- Find out if to use MaR
marshow = false;
if not(isempty(SET(no).MaR))
    marshow = true;
end;
%stateandicon=segment('iconson','hidemar');
state=stateandicon{2,1};
if state%isequal(get(DATA.Handles.hidemaricon,'state'),'on')
    marshow = false;
end;

%--- Find out if to use overlays
olshow = false;
if not(isempty(SET(no).Overlay))
    olshow = true;
    %stateandicon=segment('iconson','hidoverlay');
    state=stateandicon{3,1};

    if state%isequal(get(DATA.Handles.hideoverlayicon,'state'),'on')
        olshow = false;
    end
end

%--- Find panels that match no
%for panelloop=panelstodo;
if isequal(DATA.ViewPanelsType{panel},'mmodetemporal')
    if (SET(no).TSize<2)
        DATA.ViewPanelsType{panel} = 'one';
    end;
    
    %Make sure that mmodespatial exists.
    foundit = false;
    for loop=1:panel
        if isequal(DATA.ViewPanelsType{loop},'mmodespatial')
            foundit = true;
        end;
    end;
    
    if not(foundit)
        myfailed('Can not show only mmode image.',DATA.GUI.Segment);
        DATA.ViewPanelsType{panel} = 'one';
    end;
end;

%Make the panel visible
set(DATA.Handles.imageaxes(panel),'visible','on');

%--- Check if normal/montage
switch DATA.ViewPanelsType{panel}
    case {'one','mmodespatial','ortho','orthomip'}
        drawimageone(panel,viashow,marshow,olshow)
        %     if strcmp(DATA.ViewPanelsType{panel},'hla')
        %      updatelongaxiscontours('HLA',no,panel);
        %       set(DATA.Handles.imageaxes(panel),'plotboxaspectratio',[...
        %         SET(no).ZSize*(SET(no).SliceThickness+SET(no).SliceGap) ...
        %         SET(no).YSize*SET(no).ResolutionY ...
        %         1]);
        %       set(DATA.Handles.imageaxes(panel),'plotboxaspectratio',[...
        %         size(SET(no).HLA.IM,1)*(SET(no).SliceThickness+SET(no).SliceGap) ...
        %         size(SET(no).HLA.IM,2)*SET(no).ResolutionX ...
        %         1]);
        %     elseif strcmp(DATA.ViewPanelsType{panel},'vla')
        %      updatelongaxiscontours('VLA',no,panel)
        %       set(DATA.Handles.imageaxes(panel),'plotboxaspectratio',[...
        %         SET(no).ZSize*(SET(no).SliceThickness+SET(no).SliceGap) ...
        %         SET(no).XSize*SET(no).ResolutionX ...
        %         1]);
        %       set(DATA.Handles.imageaxes(panel),'plotboxaspectratio',[...
        %         size(SET(no).VLA.IM,1)*SET(no).ResolutionY ...
        %         size(SET(no).VLA.IM,2)*(SET(no).SliceThickness+SET(no).SliceGap) ...
        %         1]);
        %     end
    case {'hla','vla','gla','hlamip','vlamip','glamip'}
        drawimagehlavla(panel,viashow,marshow,olshow)
    case {'montage','montagerow','montagefit','sax3','montagesegmented'}
        drawimagemontage(panel,viashow,marshow)
    case 'mmodetemporal'
        drawimagemmode(no,panel);
    otherwise
        myfailed(dprintf('Unknown image viewmode %s.',DATA.ViewPanelsType{panel}),DATA.GUI.Segment);
end;

%If is this panel is the current image panel then highlight
if isequal(DATA.CurrentPanel,panel)&&...
        (DATA.GUISettings.MontageBorder||~ismember(DATA.ViewPanelsType{panel},{'montage','montagerow','montagefit','sax3'}))
    
    set(DATA.Handles.imageaxes(panel), 'xcolor',DATA.GUISettings.AxesColor,...
        'ycolor',DATA.GUISettings.AxesColor, 'linewidth',2.5, 'visible','on');
else
    set(DATA.Handles.imageaxes(panel),...
        'xcolor',[0 0 0],'ycolor',[0 0 0],...
        'linewidth',0.5,...
        'visible','off');
end;

set([...
    DATA.Handles.endocontour(panel) ...
    DATA.Handles.rvendocontour(panel) ...
    DATA.Handles.epicontour(panel) ...
    DATA.Handles.rvepicontour(panel)],'linewidth',DATA.Pref.LineWidth);

%--- Reset popup menus
set(DATA.Handles.selectcontextmenu,'Visible','off');
set(DATA.Handles.endopinmenu,'Visible','off');
set(DATA.Handles.epipinmenu,'Visible','off');
set(DATA.Handles.rvendopinmenu,'Visible','off');
set(DATA.Handles.roicontextmenu,'Visible','off');
set(DATA.Handles.measurecontextmenu,'Visible','off');
set(DATA.Handles.pointcontextmenu,'Visible','off');
set(DATA.Handles.scarcontextmenu,'Visible','off');
set(DATA.Handles.rvcontextmenu,'Visible','off');
set(DATA.Handles.lvcontextmenu,'Visible','off');
set(DATA.Handles.contrastcontextmenu,'Visible','off');
%set(DATA.Handles.rvepipinmenu,'Visible','off');

viewupdateannotext(panel);
drawcolorbar(panel);

%set(DATA.Handles.colorbar,'visible','off')

if strcmp(DATA.CurrentTool,'orthoview') && ~DATA.isorthoview
    updatetool('select');
else
    updatetool(DATA.CurrentTool,panel);%sets buttondown functions correctly according to current tool
end

drawthumbnails(isempty(DATA.DATASETPREVIEW)); %EH

if DATA.Pref.AutoSave
    %Obsoleted
    if (now-DATA.LastSaved)*24*60>5
        %More than 5 min ago.
        h = msgbox('Autosaving');
        temp = DATA.Silent;
        filemenu('savesegmentation_Callback',pwd,'autosave.seg');
        DATA.Silent = temp;
        close(h);
        DATA.LastSaved = now;
    end;
end;

%-----------------------------
function drawcontrastimage(no) %#ok<DEFNU>
%-----------------------------
%draws contrast image

global DATA NO

if nargin < 1
    no = NO;
end

panels = DATA.ViewPanels == no;
[DATA.ViewIM{panels}] = deal([]);
[DATA.Overlay(panels)] = deal(struct('alphadata',[],'cdata',[]));
%segment('makeviewim',DATA.CurrentPanel,no);
drawimageno(no);
segment('update_thumbnail',no);
if DATA.Pref.UseLight
    DATA.BalloonLevel = -1; %Force update of ballonimage
end;


%-----------------------
function drawthumbnails(calculatepreview,sliderupdated)
%-----------------------
%Draw all thumbnails. Calculatepreview is a boolean
%indicating if thumbnails needs to be redrawn.

global DATA SET
persistent setlength

if DATA.Silent || (~DATA.DataLoaded) || isempty(DATA.Handles.datasetaxes)
    return;
end;

if isempty(DATA.VisibleThumbnails)
    DATA.VisibleThumbnails=1:min(DATA.Pref.NumberVisibleThumbnails,length(SET));
    setlength=length(SET);
end

thumbsize=DATA.GUISettings.ThumbnailSize;
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
    calcfunctions('calcdatasetpreview');
end

DATA.Handles.datasetpreviewimage = image(DATA.DATASETPREVIEW,'parent',DATA.Handles.datasetaxes);
set(DATA.Handles.datasetpreviewimage,'ButtondownFcn',...
    sprintf('%s(''thumbnail_Buttondown'')','segment'));
%colormap(DATA.Colormap); %no longer needed, now truecolor /JU
axis(DATA.Handles.datasetaxes,'image','off');
ylim=[(DATA.VisibleThumbnails(1)-1) DATA.VisibleThumbnails(end)]*thumbsize+1;
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
%DATA.Handles.datasetselectline = [];
DATA.Handles.datasetflowline = [];
for loop=1:length(SET)
    ypos = (loop-1)*thumbsize+1;
    DATA.Handles.datasetflowline =  ...
        [DATA.Handles.datasetflowline ...
        plot(DATA.Handles.datasetaxes,...
        [1    1                thumbsize        thumbsize 1   ],...
        [ypos ypos+thumbsize-1 ypos+thumbsize-1 ypos      ypos],...
        'color',DATA.GUISettings.ThumbFlowLineColor)...
        ];
end;
set(DATA.Handles.datasetflowline,'visible','off');

%Show frames around linked images
if length(SET(NO).Linked) > 1
    linkies = SET(NO).Linked(SET(NO).Linked ~= NO);
    set(DATA.Handles.datasetflowline(linkies),'visible','on');
end

%draw frame around current image
ypos = (NO-1)*thumbsize+1;
DATA.Handles.datasetpreviewline =  plot(DATA.Handles.datasetaxes,...
    [1    1                thumbsize        thumbsize 1   ],...
    [ypos ypos+thumbsize-1 ypos+thumbsize-1 ypos      ypos],...
    'color',DATA.GUISettings.ThumbLineColor);
%  'linewidth',2);

% %Check current imageviewpanel
%     montagestateandicon=segment('iconson','viewall');
%      onestateandicon=segment('iconson','viewone');
%      currenttype=DATA.ViewPanelsType{DATA.CurrentPanel};
% switch currenttype
%   case 'one'
%     if ~onestateandicon{1}
%     onestateandicon{2}.isindented=1;
%     onestateandicon{2}.cdataDisplay=onestateandicon{2}.cdataIndent;
%     montagestateandicon{2}.isindented=0;
%     montagestateandicon{2}.cdataDisplay=montagestateandicon{2}.cdata;
%     DATA.Handles.permanenticonholder.render
%     end
%   case 'montage'
%     if ~montagestateandicon{1}
%     onestateandicon{2}.isindented=0;
%     onestateandicon{2}.cdataDisplay=onestateandicon{2}.cdata;
%     montagestateandicon{2}.isindented=1;
%     montagestateandicon{2}.cdataDisplay=montagestateandicon{2}.cdataIndent;
%     DATA.Handles.permanenticonholder.render
%     end
%end
%Updates viewbuttons zero is crucial else stack overflow
DATA.setviewbuttons(0)



%------------------------
function drawimageno(no)
%------------------------
%Call drawimagepanel for the panels that contain no including flow
%panels. If no is not specified then the current image stack NO is
%used. This function is typically called when new objects have been
%created or modified.

global DATA SET NO

if DATA.Silent
    return;
end;

if nargin==0
    no = NO;
end;

nos = SET(no).Linked;

ind = [];
for loop=1:length(nos)
    ind = [ind find(DATA.ViewPanels==nos(loop))]; %#ok<AGROW>
end;

ind = unique(ind); %remove duplicates.
for loop=1:length(ind)
    drawimagepanel(ind(loop));
    pause(0.05)
end;
%
% if nargin==0
%   %To fix intersection points at different image stacks.
%   for loop=1:length(DATA.ViewPanels)
%     % If not an index we've sent to imagepanel above, then slicepanel
%     if ~ismember(DATA.ViewPanels(loop),DATA.ViewPanels(ind))&&(DATA.ViewPanels(loop)~=0)
%       drawslicepanel(loop);
%     end;
%   end;
% end;

% Contains enabler/disablers that need to be updated whenever T/ZSize
% has been changed. Rather than make the check at the instances (which
% could easily be forgotten when adding a new instances), they are called
% here. /JU
%Updates viewbuttons zero is crucial else stack overflow
%DATA.setviewbuttons(0)

%segment('updateviewicons');
%DATA.updateaxestables('flowclearall');
%DATA.updateaxestables('volumeclearall');
    
%DATA.updateaxestables('measure');
%DATA.updateaxestables('volume');
%DATA.updateaxestables('flow');

showedits(no);
drawintersections;

%--------------------------
function drawintersections
%--------------------------
%Draw intersections between image visible stacks. Uses
%calcplaneintersection to calculate intersection lines.
%This function respects view settings.

global DATA

if DATA.Silent
    return;
end;

for panel=1:length(DATA.ViewPanels)
    
    no = DATA.ViewPanels(panel);
    
    %Check if valid panel
    if (DATA.ViewPanels(panel)>0)&&(ismember(DATA.ViewPanelsType{panel}, ...
            {'one','ortho','hla','vla','gla','orthomip','hlamip','vlamip'}))
        
        %--- Delete old ones if existing.
        try
            delete(DATA.Handles.planeintersectionline{panel})
            delete(DATA.Handles.glarotatehandle(panel))
        catch %#ok<CTCH>
            %do nothing
        end;
        
        %Draw new one
        stateandicon=segment('iconson','hideintersections');
        state=stateandicon{1};
        if not(state)%isequal(get(DATA.Handles.hideintersectionsicon,'state'),'off')
            h = [];
            hold(DATA.Handles.imageaxes(panel),'on');
            %current panel last to make sure orange line is not disabled
            panelseq = [1:length(DATA.ViewPanels) DATA.CurrentPanel];
            panelseq(DATA.CurrentPanel) = [];
            for loop = panelseq
                %if (DATA.ViewPanels(loop)~=no)&&...
                %    (DATA.ViewPanels(loop)>0)&&...
                %    isequal(DATA.ViewPanelsType{loop},'one')
                [x,y] = calcfunctions('calcplaneintersections',no, ...
                    DATA.ViewPanels(loop),DATA.ViewPanelsType{panel},DATA.ViewPanelsType{loop});
                h = [h;plot(DATA.Handles.imageaxes(panel),y,x,'w-');]; %#ok<AGROW>
                if isequal(loop,DATA.CurrentPanel)&&~isempty(h)&&~isempty(x)
                    set(h(end),'color',DATA.GUISettings.IntersectColor);
                end;
                if ismember(DATA.ViewPanelsType{panel},{'one','ortho'}) && ...
                        strcmp(DATA.ViewPanelsType{loop},'gla')
                    [~,i] = sort(y);
                    DATA.Handles.glarotatehandle(panel) = plot(DATA.Handles.imageaxes(panel), ...
                        [nan y(i([1 end]))*[0.3;0.7]],[nan x(i([1 end]))*[0.3;0.7]],'bd');
                    set(DATA.Handles.glarotatehandle(panel),'ButtonDownFcn', ...
                        'segment(''glarotatehandle_Buttondown'')');
                end
            end;
            hold(DATA.Handles.imageaxes(panel),'off');
            DATA.Handles.planeintersectionline{panel} = h';
            if DATA.Pref.LineWidth>0
                set(h,'linewidth',DATA.Pref.LineWidth);
            end;
            %sets the same buttondownfunction for planeintersection as for
            %imagehandle must be set here since drawintersections is called after drawintersections (which calls updatetool) in which buttondownfcn is updated.
            if ishandle(DATA.Handles.imagehandle(panel))
                set(DATA.Handles.planeintersectionline{panel},'ButtonDownFcn',get(DATA.Handles.imagehandle(panel),'ButtonDownFcn'));
            end
        else
            DATA.Handles.planeintersectionline{panel} = [];
        end;
        
    end; %Valid panel
    
end; %Loop over panel

%----------------------
function drawallslices %#ok<DEFNU>
%----------------------
%This fcn updates graphics in all visible image panels.

global DATA NO

noncurrentnos = setdiff(DATA.ViewPanels(DATA.ViewPanels > 0),NO);
for noloop = [noncurrentnos NO]
    updatenopanels(noloop);
end

%-----------------------
function drawsliceno(no)
%-----------------------
%Call updatenopanels
global NO

if nargin < 1
    no = NO;
end
updatenopanels(no);

%--------------------------
function updatenopanels(no)
%--------------------------
%Update panels containing image stack no
global DATA SET

panel = find(DATA.ViewPanels == no);
if isempty(panel)
    %No panels to update
    return
end
panelstype = DATA.ViewPanelsType(panel);
onepanels = panel(ismember(panelstype,{'one','mmodespatial','ortho'}));
mmodepanels = panel(strcmp(panelstype,'mmodespatial'));
montagepanels = panel(strcmp(panelstype,'montage'));
montagerowpanels = panel(strcmp(panelstype,'montagerow'));
montagefitpanels = panel(strcmp(panelstype,'montagefit'));
montagesegpanels = panel(strcmp(panelstype,'montagesegmented'));
montageallpanels = {montagepanels,montagerowpanels,montagefitpanels,montagesegpanels};
sax3panels = panel(strcmp(panelstype,'sax3'));
orthomippanels = panel(strcmp(panelstype,'orthomip'));
hlapanels = panel(strcmp(panelstype,'hla'));
hlamippanels = panel(strcmp(panelstype,'hlamip'));
vlapanels = panel(strcmp(panelstype,'vla'));
vlamippanels = panel(strcmp(panelstype,'vlamip'));
glapanels = panel(strcmp(panelstype,'gla'));

stateandicon=segment('iconson',{'hidescar','hidemar','hideall','play'});
%Find if viability, for drawing contour only.
viashow = false;
if not(isempty(SET(no).Scar))
    if not(isequal(SET(no).Scar.Mode,'none'))
        viashow = true;
    end;
end;
if stateandicon{1,1}% isequal(get(DATA.Handles.hidescaricon,'state'),'on')
    viashow = false;
end;

%--- Find out if to use MaR
marshow = false;
if not(isempty(SET(no).MaR))
    marshow = true;
end;
if stateandicon{2,1}%isequal(get(DATA.Handles.hidemaricon,'state'),'on')
    marshow = false;
end;

%--- Find out if to use overlays
olshow = false;
if not(isempty(SET(no).Overlay))
    olshow = true;
end;
if stateandicon{3,1}%isequal(get(DATA.Handles.hideoverlayicon,'state'),'on')
    olshow = false;
end;

if DATA.Silent
    return;
end;

%stateandicon=segment('iconson','play');
state=stateandicon{4,1};

if DATA.Run&&~state%not(isequal(get(DATA.Handles.playallicon,'state'),'on'))
    SET(no).CurrentTimeFrame = segment('getframenumber');
    if length(SET(no).Linked) > 1
        nos = SET(no).Linked;
        for nloop=1:length(nos)
            SET(nos(nloop)).CurrentTimeFrame = SET(no).CurrentTimeFrame;
        end;
    end
end;


%--- Change definition of no so that ROI and contours points to
% magnitude data set.
cno = no;
if ~isempty(SET(no).Parent)
    no = SET(no).Parent;
end;

%Do things that require iteration
for panelloop = panel
    %Make view im
    if isempty(DATA.ViewIM{panelloop})
        segment('makeviewim',panelloop,no);
    end
    
    %--- Viability contour
    if (viashow)
        delete(DATA.Handles.scarcontour{panelloop});
        try
            delete(DATA.Handles.weightedscarcontour{panelloop});
            delete(DATA.Handles.moextentcontour{panelloop});
            delete(DATA.Handles.mocontour{panelloop});
        catch %#ok<CTCH>
        end;
        drawviabilityhelper(no,panelloop);
    end
    
    %--- Mar contour
    if (marshow)
        try
            delete(DATA.Handles.marcontour(panelloop));
            DATA.Handles.marcontour{panelloop} = [];
        catch %#ok<CTCH>
        end;
        drawmarhelper(no,panelloop);
    end;
    
    
    if ismember(panelloop,onepanels)
        %Atrial contour
        if ~isempty(SET(no).AtrialScar)
            atrialscar('drawsliceone',no,panel);
        end;
    end
    
    drawcolorbar;
end;

%--- Do different update depending on viewMode
if ~isempty(onepanels)
    if size(DATA.ViewIM{onepanels(1)},3) == 1
        if size(DATA.ViewIM{onepanels(1)},4) == 1
            set(DATA.Handles.imagehandle(onepanels),'cdata',...
                DATA.ViewIM{onepanels(1)}(:,:,:));
        else
            set(DATA.Handles.imagehandle(onepanels),'cdata',...
                squeeze(DATA.ViewIM{onepanels(1)}(:,:,:,SET(no).CurrentSlice,:)));
        end
    else
        if size(DATA.ViewIM{onepanels(1)},4) == 1
            set(DATA.Handles.imagehandle(onepanels),'cdata',...
                squeeze(DATA.ViewIM{onepanels(1)}(:,:,SET(no).CurrentTimeFrame,:)));
        else
            set(DATA.Handles.imagehandle(onepanels),'cdata',...
                squeeze(DATA.ViewIM{onepanels(1)}(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice,:)));
        end
    end
    
    if olshow
        %Need to redraw and update panel frame if overlay
        drawimageone(panel,viashow,marshow,olshow);
        set(DATA.Handles.imageaxes,...
            'xcolor',[0 0 0],'ycolor',[0 0 0]);
        set(DATA.Handles.imageaxes(DATA.CurrentPanel),...
            'xcolor',DATA.GUISettings.AxesColor,'ycolor',DATA.GUISettings.AxesColor,'linewidth',2.5);
        return
    end
    
    %--- Endo and epi contours
    updatecontours(no,onepanels);
    
    %--- Pins
    %updatepins(no,onepanels);
    
    %--- Interpolation points
    updateinterp(no,onepanels);
    
    %Draw intersection with segmentation
    updateintersectionpoints(onepanels);
    
    %--- Draw rois
    updaterois(no,onepanels);
    
    %--- Display points
    updateannotationpoints(no,onepanels);
    
    %--- Plot measures, time specific
    updatemeasures(no,onepanels);
    
    %Update slicetimetext
    stri = dprintf('Slice:%02d Time:%03d ms',SET(no).CurrentSlice,round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
    set(DATA.Handles.slicetimetext(onepanels),'string',stri);
    viewupdatetextposition(onepanels);
    viewupdateannotext(onepanels);
end
for i = 1:4
    %loop over montage, montagerow, montagefit, montagesegmented panels
    montagepanels = montageallpanels{i};
    if ~isempty(montagepanels)
        %%%%%%%%%%%%%%%
        %%% Montage %%%
        %%%%%%%%%%%%%%%
        
        updatemodeldisplay(no,montagepanels(1));
        set(DATA.Handles.imagehandle(montagepanels),'cdata',squeeze(DATA.ViewIM{montagepanels(1)}(:,:,SET(no).CurrentTimeFrame,:)));
        
        % --- Endo and epi contours
        updatemontagecontours(no,montagepanels);
        
        %--- Pins
        updatemontagepins(no,montagepanels);
        
        %--- Interp Pts
        updatemontageinterp(no,montagepanels);
        
        %--- Display ROI
        updatemontagerois(no,montagepanels);
        
        %--- Display points
        updatemontagepoints(no,montagepanels);
        
        %--- Plot measures, time specific, JU
        updatemontagemeasures(no,montagepanels);
        
        %Update slicetimetext
        stri = dprintf('Time:%03d ms',round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
        set(DATA.Handles.slicetimetext(montagepanels),'string',stri);
        viewupdatetextposition(montagepanels);
        viewupdateannotext(montagepanels);
    end
end

if ~isempty(sax3panels)
    %%%%%%%%%%%%%%
    %%%% SAX3 %%%%
    %%%%%%%%%%%%%%
    
    DATA.updatesax3display(no,sax3panels(1));
    set(DATA.Handles.imagehandle(sax3panels),'cdata',squeeze(DATA.ViewIM{sax3panels(1)}(:,:,SET(no).CurrentTimeFrame,:)));
    
    % --- Endo and epi contours
    updatesax3contours(no,sax3panels);
    
    %Update slicetimetext
    stri = dprintf('Slice:%02d Time:%03d ms',SET(no).CurrentSlice,round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
    set(DATA.Handles.slicetimetext(sax3panels),'string',stri);
    viewupdatetextposition(sax3panels);
    viewupdateannotext(sax3panels);
    
end; %Swith type of panel

if ~isempty(hlapanels)
    %Image data
    set(DATA.Handles.imagehandle(hlapanels),'cdata',squeeze(DATA.ViewIM{hlapanels(1)}(:,:,SET(no).CurrentTimeFrame,:)));
    
    %Draw intersection with segmentation
    updateintersectionpoints(hlapanels);
    
    %Update measures
    updatemeasures(no,hlapanels);
    
    %Update slicetimetext
    stri = dprintf('Slice:%02d Time:%03d ms',SET(no).HLA.slice,round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
    set(DATA.Handles.slicetimetext(hlapanels),'string',stri);
    viewupdatetextposition(hlapanels);
    viewupdateannotext(hlapanels);
end
if ~isempty(vlapanels)
    %Image data
    set(DATA.Handles.imagehandle(vlapanels),'cdata',squeeze(DATA.ViewIM{vlapanels(1)}(:,:,SET(no).CurrentTimeFrame,:)));
    
    %Draw intersection with segmentation
    updateintersectionpoints(vlapanels);
    
    %Update measures
    updatemeasures(no,vlapanels);
    
    %Update slicetimetext
    stri = dprintf('Slice:%02d Time:%03d ms',SET(no).VLA.slice,round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
    set(DATA.Handles.slicetimetext(vlapanels),'string',stri);
    viewupdatetextposition(vlapanels);
    viewupdateannotext(vlapanels);
end

mippanels = {orthomippanels,hlamippanels,vlamippanels};
for i = 1:3
    panels = mippanels{i};
    if ~isempty(panels)
        %Image data
        set(DATA.Handles.imagehandle(panels),'cdata',squeeze(DATA.ViewIM{panels(1)}(:,:,SET(no).CurrentTimeFrame)));
        
        %Update slicetimetext
        stri = dprintf('Time:%03d ms',round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
        set(DATA.Handles.slicetimetext(panels),'string',stri);
        viewupdatetextposition(panels);
    end
end
if ~isempty(glapanels)
    %Image data
    imsz = size(DATA.ViewIM{glapanels(1)});
    xsz = SET(no).ZSize;
    ysz = size(DATA.ViewIM{glapanels(1)},2); %floor(abs(SET(pno).YSize*cos(glaangle))+SET(pno).XSize*sin(glaangle));
    glaangle = SET(no).GLA.angle;
    xres = SET(no).SliceThickness + SET(no).SliceGap;
    yres = SET(no).ResolutionY*cos(glaangle)+SET(no).ResolutionX*abs(sin(glaangle));
    updateglazoomstate(no,imsz(2));
    set(DATA.Handles.imagehandle(glapanels), ...
        'cdata',squeeze(DATA.ViewIM{glapanels(1)}(:,:,SET(no).CurrentTimeFrame,:)), ...
        'XData',[0.5 imsz(2)+0.5],'YData',[0.5 imsz(1)+0.5]);
    set(DATA.Handles.imageaxes(glapanels),...
        'xlim',SET(no).GLA.ZoomState(1:2),...
        'ylim',SET(no).GLA.ZoomState(3:4),...
        'plotboxaspectratio',[ysz*yres xsz*xres 1], ...
        'dataaspectratio',[1/yres 1/xres 1]);
    
    %Draw intersection with segmentation
    updateintersectionpoints(glapanels);
    
    %Update measures
    updatemeasures(no,glapanels);
    
    %Update slicetimetext
    stri = dprintf('Time:%03d ms',round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
    set(DATA.Handles.slicetimetext(glapanels),'string',stri);
    viewupdatetextposition(glapanels);
    viewupdateannotext(glapanels);
end

%--- Update timebar

%Only update time slider and time text on current image stack
if SET(no).TSize>1
  %DATA.updatetimebaraxes;
  t = SET(no).TimeVector*1000;
  %profile on
  if no == DATA.LVNO
    set(DATA.Handles.timebarlv,'xdata',...
      t(SET(no).CurrentTimeFrame)*[1 1])
    
  end
  if no == DATA.FlowNO
    set(DATA.Handles.timebarflow,'xdata',...
      t(SET(no).CurrentTimeFrame)*[1 1])
    
  end
  %profile report
  set(DATA.Handles.timebar,'xdata',...
    t(SET(no).CurrentTimeFrame)*[1 1]);
end;

showedits(no);
%updatetool(DATA.CurrentTool,panel)
updatevisibility;

%-----------------------------------------------
function drawimagemontage(panel,viashow,marshow)
%-----------------------------------------------
%Main workhorse for creating montage view.
% - viashow is if viability data should be updated.
% - marshow is if MaR data should be updated.

global DATA SET

pno = DATA.ViewPanels(panel);
type = DATA.ViewPanelsType{panel};

if (~isempty(SET(pno).Scar))&&(SET(pno).ZSize>100)
    if yesno('Viewing scar in all slices might take very long time. Revert to one slice?',[],DATA.GUI.Segment);
        DATA.ViewPanelsType{panel} = 'one';
        drawimagepanel(panel);
        return;
    end;
end;

no = pno;
if ~isempty(SET(pno).Parent)
    no = SET(pno).Parent;
end;

if isempty(DATA.ViewIM{panel})
    segment('makeviewim',panel,no);
    if ismember(type,{'montagefit','sax3'})
        segment('updatemodeldisplay');
    end
elseif ismember(type,{'montagefit','sax3'}) && ...
        (size(DATA.ViewIM{panel},2)-size(DATA.ViewIM{panel},1)) * ...
        (DATA.ViewMatrix(1)-DATA.ViewMatrix(2)) < 0
    %Wrong orientation of view, make new viewim
    segment('makeviewim',panel,no);
    segment('updatemodeldisplay');
end;
if strcmp(type,'sax3')
    DATA.updatesax3display(no,panel);
else
    updatemodeldisplay(no,panel);
end

DATA.Handles.imagehandle(panel) = image(squeeze(DATA.ViewIM{panel}(:,:,SET(no).CurrentTimeFrame,:)),...
    'parent',DATA.Handles.imageaxes(panel));
DATA.Handles.scarcontour{panel} = [];
DATA.Handles.weightedscarcontour{panel} = [];
DATA.Handles.mocontour{panel} = [];
DATA.Handles.moextentcontour{panel} = [];
DATA.Handles.marcontour{panel} = [];
%colormap(DATA.Handles.imageaxes(panel),DATA.Colormap);

set(DATA.Handles.imageaxes(panel),...
    'color',[0 0 0],...
    'xtick',[],...
    'ytick',[],...
    'plotboxaspectratio',[...
    size(DATA.ViewIM{panel},2)*SET(no).ResolutionY ...
    size(DATA.ViewIM{panel},1)*SET(no).ResolutionX 1]);

%Update zoom state
switch DATA.ViewPanelsType{panel}
    case 'montage'
        if isempty(SET(no).MontageZoomState)
            SET(no).MontageZoomState = [0.5;size(DATA.ViewIM{panel},2)+0.5;0.5;size(DATA.ViewIM{panel},1)+0.5];
        end;
        set(DATA.Handles.imageaxes(panel),...
            'xlim',SET(no).MontageZoomState(1:2),...
            'ylim',SET(no).MontageZoomState(3:4));
    case 'montagerow'
        %Montagerow
        if isempty(SET(no).MontageRowZoomState)
            SET(no).MontageRowZoomState = [0.5;size(DATA.ViewIM{panel},2)+0.5;0.5;size(DATA.ViewIM{panel},1)+0.5];
        end;
        set(DATA.Handles.imageaxes(panel),...
            'xlim',SET(no).MontageRowZoomState(1:2),...
            'ylim',SET(no).MontageRowZoomState(3:4));
    case 'montagesegmented'
        %Montagesegmented
        if isempty(SET(no).MontageFitZoomState)
            SET(no).MontageFitZoomState = [0.5;size(DATA.ViewIM{panel},2)+0.5;0.5;size(DATA.ViewIM{panel},1)+0.5];
        end;
        set(DATA.Handles.imageaxes(panel),...
            'xlim',SET(no).MontageFitZoomState(1:2),...
            'ylim',SET(no).MontageFitZoomState(3:4));
    case {'montagefit','sax3'}
        %Montagefit
        viewimsz = size(DATA.ViewIM{panel});
        if isempty(SET(no).MontageFitZoomState)
            %ugly hack to get size of panel
            set(DATA.Handles.imageaxes(panel),'Units','points');
            boxpos = get(DATA.Handles.imageaxes(panel),'Position');
            set(DATA.Handles.imageaxes(panel),'Units','normalized');
            
            if viewimsz(1) < viewimsz(2)
                ratio = (viewimsz(1)/boxpos(4))/(viewimsz(2)/boxpos(3));
                SET(no).MontageFitZoomState(1:2) = [0.5;viewimsz(2)*ratio+0.5];
                SET(no).MontageFitZoomState(3:4) = [0.5;viewimsz(1)+0.5];
            else
                ratio = (viewimsz(2)/boxpos(3))/(viewimsz(1)/boxpos(4));
                SET(no).MontageFitZoomState(1:2) = [0.5;viewimsz(2)+0.5];
                SET(no).MontageFitZoomState(3:4) = [0.5;viewimsz(1)*ratio+0.5];
            end
        end
        
        %     if isempty(SET(no).MontageFitZoomState)
        %       SET(no).MontageFitZoomState = [0.5;viewimsz(2)*ratio2+0.5;0.5;viewimsz(1)*ratio1+0.5];
        %     end;
        
        mfzs = SET(no).MontageFitZoomState;
        xdif = mfzs(2)-mfzs(1);
        ydif = mfzs(4)-mfzs(3);
        %Display has changed, turn zoom state 90 degrees
        if (viewimsz(2)-viewimsz(1))*(xdif-ydif) < 0
            SET(no).MontageFitZoomState = mfzs([3 4 1 2]);
            mfzs = SET(no).MontageFitZoomState;
            xdif = mfzs(2)-mfzs(1);
            ydif = mfzs(4)-mfzs(3);
        end
        xlim = mean(mfzs(1:2)) + xdif*[-0.6 0.6];
        ylim = mean(mfzs(3:4)) + ydif*[-0.6 0.6];
        set(DATA.Handles.imageaxes(panel),...
            'xlim',xlim,...
            'ylim',ylim,...
            'plotboxaspectratio',[xdif * SET(no).ResolutionY ...
            ydif * SET(no).ResolutionX 1]);
end;

set(DATA.Handles.imageaxes(panel),'dataaspectratio',...
    [1/SET(no).ResolutionY ...
    1/SET(no).ResolutionX 1]);

%Draw box around slices
DATA.Handles.selectslicehandle{panel} = zeros(1,SET(no).ZSize);
hold(DATA.Handles.imageaxes(panel),'on');

%On a detailed scale, this box is in fact wrongly placed, but otherwise the
%sliceline ends up below other lines. /JU
if strcmp(DATA.ViewPanelsType{panel},'montagesegmented')
    slicestoinclude = find(findfunctions('findslicewithendo',no))';
    if min(slicestoinclude) > 1
        slicestoinclude = [min(slicestoinclude)-1 slicestoinclude];
    end
    if max(slicestoinclude) < SET(no).ZSize;
        slicestoinclude = [slicestoinclude max(slicestoinclude)+1];
    end
else
    slicestoinclude = 1:SET(no).ZSize;
end
for zloop=1:numel(slicestoinclude)
    loop = slicestoinclude(zloop);
    x1 = 2.25+mod((zloop-1),DATA.ViewPanelsMatrix{panel}(2))*SET(no).YSize;
    x2 = -2.5+x1+SET(no).YSize;
    y1 = 2.25+floor((zloop-1)/DATA.ViewPanelsMatrix{panel}(2))*SET(no).XSize;
    y2 = -2.5+y1+SET(no).XSize;
    DATA.Handles.selectslicehandle{panel}(loop) = plot(DATA.Handles.imageaxes(panel),...
        [x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],DATA.GUISettings.SliceLineSpec);
    set(DATA.Handles.selectslicehandle{panel}(loop),'visible','off','linewidth',3);
end;

%Activate currentslice
if isempty(SET(no).StartSlice)
    SET(no).StartSlice = NaN;
end;
if isnan(SET(no).StartSlice)||(SET(no).StartSlice<0)||(SET(no).StartSlice>SET(no).ZSize)
    SET(no).StartSlice=1;
end;
if isempty(SET(no).EndSlice)
    SET(no).EndSlice = NaN;
end;
if isnan(SET(no).EndSlice)||(SET(no).EndSlice<0)||(SET(no).EndSlice>SET(no).ZSize)
    SET(no).EndSlice=SET(no).ZSize;
end;

set(DATA.Handles.selectslicehandle{panel}(SET(no).StartSlice:SET(no).EndSlice),'visible','on');

%--- Write text describing image type
drawmontageimagetypetext(no,pno,panel);

%Add contours and pins

drawmontagecontours(no,panel);
if strcmp(type,'sax3')
    updatesax3contours(no,panel);
end

%Intersections are not visible now.
DATA.Handles.planeintersectionline{panel} = [];
DATA.Handles.endointersectionline{panel} = [];
DATA.Handles.endointersectionpoints{panel} = [];
DATA.Handles.epiintersection{panel} = [];
DATA.Handles.rvendointersection{panel} = [];
DATA.Handles.rvepiintersection{panel} = [];

% --- Viability contour
if (viashow)
    drawviabilityhelper(pno,panel);
end

% --- MaR contour
if (marshow)
    drawmarhelper(pno,panel);
end

%Endo & Epi pins
drawmontagepins(no,panel);

%Endo & Epi interp
drawmontageinterp(no,panel);

%Atrial contour
if ~isempty(SET(no).AtrialScar)
    atrialscar('drawimagemontage',no,panel);
end;

%--- Add measures
drawmontagemeasures(no,panel);

%--- Add points
drawmontagepoints(no,panel);

%--- Plot roi's if existing
DATA.plotrois(panel,no);

%--- Plot center
hold(DATA.Handles.imageaxes(panel),'on');
[tempx,tempy] = ndgrid(0:(DATA.ViewPanelsMatrix{panel}(1)-1),0:(DATA.ViewPanelsMatrix{panel}(2)-1));
DATA.Handles.center{panel} = plot(...
    DATA.Handles.imageaxes(panel),...
    SET(no).CenterY+tempy(:)*SET(no).YSize,...
    SET(no).CenterX+tempx(:)*SET(no).XSize,'+',...
    'color',DATA.centercrossdef);

clear tempx tempy


%--- Plot cursor
DATA.Handles.cursor(panel) = plot(...
    DATA.Handles.imageaxes(panel),...
    SET(no).CenterY,SET(no).CenterX,'r-');
set(DATA.Handles.cursor(panel),'visible','off');
hold(DATA.Handles.imageaxes(panel),'off');

updatevisibility;

%--------------------------------
function drawimagemmode(no,panel)
%--------------------------------
%Draw temporal mmode view
global DATA SET

% Show the temporal part of the mmode view

DATA.Handles.imagehandle(panel) = image([-1 2]*SET(no).TSize,[0 SET(no).XSize], ...
    zeros(SET(no).XSize,3*SET(no).TSize,3),'parent',...
    DATA.Handles.imageaxes(panel));

hold(DATA.Handles.imageaxes(panel),'on');
temp = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
hold(DATA.Handles.imageaxes(panel),'off');

%Assign to an "empty handle"
DATA.Handles.endocontour(panel) = temp;
DATA.Handles.epicontour(panel) = temp;
DATA.Handles.rvendocontour(panel) = temp;
DATA.Handles.rvepicontour(panel) = temp;
DATA.Handles.cursor(panel) = temp;
DATA.Handles.roicontour{panel} = temp;
DATA.Handles.roitext{panel} = temp;
DATA.Handles.scarcontour{panel} = temp;
DATA.Handles.weightedscarcontour{panel} = temp;
DATA.Handles.mocontour{panel} = temp;
DATA.Handles.moextentcontour{panel} = temp;
DATA.Handles.marcontour{panel} = temp;
DATA.Handles.pointp{panel} = temp;
DATA.Handles.pointo{panel} = temp;
DATA.Handles.pointtext{panel} = temp;
DATA.Handles.measureline{panel} = {temp};
DATA.Handles.measuretext{panel} = temp;
DATA.Handles.planeintersectionline{panel} = [];

%--- Add lines in mmode display
hold(DATA.Handles.imageaxes(panel),'on');
DATA.Handles.mmode1line = plot(DATA.Handles.imageaxes(panel),...
    [0 SET(no).TSize+1],[SET(no).XSize/2-10 SET(no).XSize/2-10],'w:');
DATA.Handles.mmode2line = plot(DATA.Handles.imageaxes(panel),...
    [0 SET(no).TSize+1],[SET(no).XSize/2+10 SET(no).XSize/2+10],'w:');
DATA.Handles.mmodetimebar1 = plot(DATA.Handles.imageaxes(panel),[1 1],[0 SET(no).XSize],'w:');
DATA.Handles.mmodetimebar2 = plot(DATA.Handles.imageaxes(panel),[SET(no).TSize SET(no).TSize],[0 SET(no).XSize],'w:');
hold(DATA.Handles.imageaxes(panel),'off');

%--- Button downs
set(DATA.Handles.mmode1line,'ButtonDownFcn',...
    sprintf('%s(''mmode1line_Buttondown'',%d)','segment',panel));
set(DATA.Handles.mmode2line,'ButtonDownFcn',...
    sprintf('%s(''mmode2line_Buttondown'',%d)','segment',panel));
set(DATA.Handles.mmodetimebar1,'ButtonDownFcn',...
    sprintf('%s(''mmodetimebar1_Buttondown'',%d)','segment',panel));
set(DATA.Handles.mmodetimebar2,'ButtonDownFcn',...
    sprintf('%s(''mmodetimebar2_Buttondown'',%d)','segment',panel));

set(DATA.Handles.imageaxes(panel),...
    'color',[0 0 0],...
    'xtick',[],...
    'ytick',calcfunctions('calcticks',SET(no).XSize,SET(no).ResolutionX),...
    'xticklabel',[],'yticklabel',[]);

segment('updatemmode'); %Displays mmode image
segment('updatemmodevisibility'); %Hide parts
segment('updatemmodeline'); %Draw correct mmode line

%--Main image buttondown
pause(0.05);
updatetool;
updatevisibility;

%--------------------------------------------------
function drawimageone(panel,viashow,marshow,olshow)
%--------------------------------------------------
%Main workhorse for creating view of one image slice.
%- viashow is whether to show viability
%- marshow is whether to show MaR
%- olshow is whether to show overlays or not.

global DATA SET

%%%%%%%%%%%%%%%
%%% ONE %%%%%%%
%%%%%%%%%%%%%%%

pno = DATA.ViewPanels(panel);
no = pno;

if isempty(DATA.ViewIM{panel})
    segment('makeviewim',panel,no);
end;

%
if ~isempty(SET(pno).Parent)
    no = SET(pno).Parent;
end;

if length(SET(pno).Linked)<2
    linkaxes(DATA.Handles.imageaxes(panel),'off');
end;

%--- Draw image
axes(DATA.Handles.imageaxes(panel)); %Make sure draw in correct figure.
hold off

xsz = SET(pno).XSize;
ysz = SET(pno).YSize;
xres = SET(pno).ResolutionX;
yres = SET(pno).ResolutionY;

if size(DATA.ViewIM{panel},3) == 1
    if size(DATA.ViewIM{panel},4) == 1
        DATA.Handles.imagehandle(panel) = image([0.5 ysz+0.5],[0.5 xsz+0.5],DATA.ViewIM{panel}(:,:,:),'parent',DATA.Handles.imageaxes(panel));
    else
        DATA.Handles.imagehandle(panel) = image([0.5 ysz+0.5],[0.5 xsz+0.5],squeeze(DATA.ViewIM{panel}(:,:,:,SET(no).CurrentSlice,:)),'parent',DATA.Handles.imageaxes(panel));
    end
else
    if size(DATA.ViewIM{panel},4) == 1
        DATA.Handles.imagehandle(panel) = image([0.5 ysz+0.5],[0.5 xsz+0.5],squeeze(DATA.ViewIM{panel}(:,:,SET(no).CurrentTimeFrame,:)),'parent',DATA.Handles.imageaxes(panel));
    else
        DATA.Handles.imagehandle(panel) = image([0.5 ysz+0.5],[0.5 xsz+0.5],squeeze(DATA.ViewIM{panel}(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice,:)),'parent',DATA.Handles.imageaxes(panel));
    end
end
set(DATA.Handles.imageaxes(panel),...
    'color',DATA.GUISettings.BackgroundColor,...
    'xtick',calcfunctions('calcticks',ysz,yres),...
    'ytick',calcfunctions('calcticks',xsz,xres),...
    'xticklabel',[],'yticklabel',[]);
hold(DATA.Handles.imageaxes(panel),'on');
DATA.Handles.scarcontour{panel} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
DATA.Handles.weightedscarcontour{panel} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
DATA.Handles.mocontour{panel} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
DATA.Handles.moextentcontour{panel} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
DATA.Handles.marcontour{panel} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
hold(DATA.Handles.imageaxes(panel),'off');

set(DATA.Handles.imageaxes(panel),'plotboxaspectratio',[...
    ysz*yres ...
    xsz*xres 1]);


%--- Update zoom state
if isempty(SET(no).NormalZoomState)
    SET(no).NormalZoomState = [0.5;ysz+0.5;0.5;xsz+0.5];
end;
set(DATA.Handles.imageaxes(panel),...
    'xlim',SET(no).NormalZoomState(1:2),...
    'ylim',SET(no).NormalZoomState(3:4));

set(DATA.Handles.imageaxes(panel),'dataaspectratio',...
    [1/yres ...
    1/xres 1]);

%--- Turn on hold
hold(DATA.Handles.imageaxes(panel),'on');

%--- Overlay
if (olshow)
    overlay('draw',pno,panel,DATA.CurrentTheme);
end

%Image type text
drawimagetypetext(no,pno,panel);

%Draw contours of Endo, Epi, RVEndo and RVEpi.
drawcontours(no,panel);

%Atrial contour
if ~isempty(SET(no).AtrialScar)
    atrialscar('drawimageone',no,panel);
end;

%--- Viability contour
if (viashow)
    drawviabilityhelper(pno,panel);
end

%--- MaR contour
if (marshow)
    drawmarhelper(pno,panel);
end;

%--- plot Interp points
drawinterp(no,panel,olshow);

%--- plot pins
drawpins(no,panel);

if ~ismember(DATA.ViewPanelsType{panel},{'orthomip'})
    %Plot intersection points
    drawintersectionpoints(no,panel);
end

%Used for ROI's etc.
DATA.Handles.cursor(panel) = plot(DATA.Handles.imageaxes(panel), ...
    SET(no).CenterY,...
    SET(no).CenterX,'y-');
set(DATA.Handles.cursor(panel),'visible','off');

%Plot center cross
DATA.Handles.center{panel} = plot(DATA.Handles.imageaxes(panel),...
    SET(no).CenterY, ...
    SET(no).CenterX,'+', ...
    'color',DATA.centercrossdef);

%--- Plot intersecting planes
DATA.Handles.planeintersectionline{panel} = [];

%--- Plot measures
drawmeasures(no,panel);

%--- Add points
drawannotationpoints(no,panel)

%Add selected frame
DATA.Handles.selectslicehandle{panel} = [];

%--- Fix axis and color of image
%colormap(DATA.Colormap);

nop = DATA.ViewPanels(panel);
if length(SET(nop).Linked)>1
    %Find what panels are the different data sets
    temppanels = SET(nop).Linked;
    panels = [];
    for loop=1:length(temppanels)
        temp = find(DATA.ViewPanels==temppanels(loop));
        panels = [panels;temp(:)]; %#ok<AGROW>
    end;
    linkaxes(DATA.Handles.imageaxes(panels),'xy');
else
    linkaxes(DATA.Handles.imageaxes(panel),'off');
end;

%--- Plot roi's if existing
DATA.drawroiinpanel(panel);

if isequal(DATA.ViewPanelsType{panel},'mmodespatial')
    %--- Add mmode line
    hold(DATA.Handles.imageaxes(panel),'on');
    DATA.Handles.mmodeline = plot(DATA.Handles.imageaxes(panel),...
        [SET(no).Mmode.X SET(no).Mmode.X+1],...
        [SET(no).Mmode.Y SET(no).Mmode.Y+1],'w-');
    DATA.Handles.mmodepoint1 = plot(DATA.Handles.imageaxes(panel),SET(no).Mmode.X,SET(no).Mmode.Y,'bo');
    DATA.Handles.mmodepoint2 = plot(DATA.Handles.imageaxes(panel),SET(no).Mmode.X,SET(no).Mmode.Y,'bo');
    DATA.Handles.mmodepointcenter = plot(DATA.Handles.imageaxes(panel),SET(no).Mmode.X,SET(no).Mmode.Y,'bd');
    DATA.Handles.mmodempoint1 = plot(DATA.Handles.imageaxes(panel),SET(no).Mmode.X,SET(no).Mmode.Y,'w+');
    DATA.Handles.mmodempoint2 = plot(DATA.Handles.imageaxes(panel),SET(no).Mmode.X,SET(no).Mmode.Y,'w+');
    hold(DATA.Handles.imageaxes(panel),'off');
    
    %Set up button down things in mmode
    set(DATA.Handles.mmodepointcenter,'ButtonDownFcn',...
        sprintf('%s(''mmodecenter_Buttondown'',%d)','segment',panel));
    set(DATA.Handles.mmodepoint1,'ButtonDownFcn',...
        sprintf('%s(''mmode1_Buttondown'',%d)','segment',panel));
    set(DATA.Handles.mmodepoint2,'ButtonDownFcn',...
        sprintf('%s(''mmode2_Buttondown'',%d)','segment',panel));
    set(DATA.Handles.mmodempoint1,'ButtonDownFcn',...
        sprintf('%s(''mmodepoint1_Buttondown'',%d)','segment',panel));
    set(DATA.Handles.mmodempoint2,'ButtonDownFcn',...
        sprintf('%s(''mmodepoint2_Buttondown'',%d)','segment',panel));
    
end;

%--- Add centerline if rotated data set
hold(DATA.Handles.imageaxes(panel),'on');
if SET(no).Rotated
    DATA.Handles.centerline = plot(DATA.Handles.imageaxes(panel),...
        [SET(no).RotationCenter SET(no).RotationCenter],...
        [0 SET(no).XSize],'c:');
end;
hold(DATA.Handles.imageaxes(panel),'off');

set(DATA.imagefig,'WindowButtonUpFcn',...
    sprintf('%s(''buttonup_Callback'')','segment'));

updatevisibility;

%-----------------------------------------------------
function drawimagehlavla(panel,viashow,marshow,olshow)
%-----------------------------------------------------
%Main workhorse for creating view of one image slice.
%- viashow is whether to show viability
%- marshow is whether to show MaR
%- olshow is whether to show overlays or not.

global DATA SET

%%%%%%%%%%%%%%%
%%% HLA/VLA %%%
%%%%%%%%%%%%%%%

pno = DATA.ViewPanels(panel);
no = pno;
type = DATA.ViewPanelsType{panel};

if isempty(DATA.ViewIM{panel})
    segment('makeviewim',panel,no);
end;

if ~isempty(SET(pno).Parent)
    no = SET(pno).Parent;
end;

if length(SET(pno).Linked)<2
    linkaxes(DATA.Handles.imageaxes(panel),'off');
end;

%--- Draw image
axes(DATA.Handles.imageaxes(panel)); %Make sure draw in correct figure.
hold off

switch type(1:3)
    case 'hla'
        xsz = SET(pno).ZSize;
        ysz = SET(pno).YSize;
        xres = SET(pno).SliceThickness + SET(pno).SliceGap;
        yres = SET(pno).ResolutionY;
    case 'vla'
        xsz = SET(pno).ZSize;
        ysz = SET(pno).XSize;
        xres = SET(pno).SliceThickness + SET(pno).SliceGap;
        yres = SET(pno).ResolutionX;
    case 'gla'
        %General longaxis. glaangle = 0 <=> HLA; glaangle = pi/2 <=> VLA
        glaangle = SET(pno).GLA.angle;
        xsz = SET(pno).ZSize;
        ysz = size(DATA.ViewIM{panel},2); %floor(abs(SET(pno).YSize*cos(glaangle))+SET(pno).XSize*sin(glaangle));
        xres = SET(pno).SliceThickness + SET(pno).SliceGap;
        %yres = SET(pno).ResolutionY*cos(glaangle)^2+SET(pno).ResolutionX*sin(glaangle)^2;
        yres = SET(pno).ResolutionY*cos(glaangle)+SET(pno).ResolutionX*abs(sin(glaangle));
end
if size(DATA.ViewIM{panel},3) == 1
    if size(DATA.ViewIM{panel},4) == 1
        DATA.Handles.imagehandle(panel) = image([0.5 ysz+0.5],[0.5 xsz+0.5],DATA.ViewIM{panel}(:,:,:),'parent',DATA.Handles.imageaxes(panel));
    else
        DATA.Handles.imagehandle(panel) = image([0.5 ysz+0.5],[0.5 xsz+0.5],squeeze(DATA.ViewIM{panel}(:,:,:,SET(no).CurrentSlice,:)),'parent',DATA.Handles.imageaxes(panel));
    end
else
    if size(DATA.ViewIM{panel},4) == 1
        DATA.Handles.imagehandle(panel) = image([0.5 ysz+0.5],[0.5 xsz+0.5],squeeze(DATA.ViewIM{panel}(:,:,SET(no).CurrentTimeFrame,:)),'parent',DATA.Handles.imageaxes(panel));
    else
        DATA.Handles.imagehandle(panel) = image([0.5 ysz+0.5],[0.5 xsz+0.5],squeeze(DATA.ViewIM{panel}(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice,:)),'parent',DATA.Handles.imageaxes(panel));
    end
end
set(DATA.Handles.imageaxes(panel),...
    'color',DATA.GUISettings.BackgroundColor,...
    'xtick',calcfunctions('calcticks',ysz,yres),...
    'ytick',calcfunctions('calcticks',xsz,xres),...
    'xticklabel',[],'yticklabel',[]);
hold(DATA.Handles.imageaxes(panel),'on');
DATA.Handles.scarcontour{panel} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
DATA.Handles.weightedscarcontour{panel} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
DATA.Handles.mocontour{panel} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
DATA.Handles.moextentcontour{panel} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
DATA.Handles.marcontour{panel} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
hold(DATA.Handles.imageaxes(panel),'off');

set(DATA.Handles.imageaxes(panel),'plotboxaspectratio',[...
    ysz*yres ...
    xsz*xres 1]);

%--- Turn on hold
hold(DATA.Handles.imageaxes(panel),'on');

%Image type text
drawimagetypetext(no,pno,panel);

if ismember(type,{'hla','vla','gla'})
    %Plot intersection points
    drawintersectionpoints(no,panel);
    
    %Plot measures
    drawmeasures(no,panel);
end

%Used for ROI's etc.
DATA.Handles.cursor(panel) = plot(DATA.Handles.imageaxes(panel),...
    SET(no).CenterY,...
    SET(no).CenterX,'y-');
set(DATA.Handles.cursor(panel),'visible','off');

DATA.Handles.center{panel} = plot(DATA.Handles.imageaxes(panel),...
    SET(no).CenterY, ...
    SET(no).CenterX,'+', ...
    'color',DATA.centercrossdef);

%--- Set up contour handles
DATA.Handles.endocontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,'r');
DATA.Handles.epicontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,'g');
DATA.Handles.rvendocontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,'m');
DATA.Handles.rvepicontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,'c');

%--- Plot intersecting planes
DATA.Handles.planeintersectionline{panel} = [];

%--- Set up points handles
DATA.Handles.pointp{panel} = [];
DATA.Handles.pointo{panel} = [];
DATA.Handles.pointtext{panel} = [];

%Add selected frame
DATA.Handles.selectslicehandle{panel} = [];

nop = DATA.ViewPanels(panel);
if length(SET(nop).Linked)>1
    %Find what panels are the different data sets
    temppanels = SET(nop).Linked;
    panels = [];
    for loop=1:length(temppanels)
        temp = find(DATA.ViewPanels==temppanels(loop));
        panels = [panels;temp(:)]; %#ok<AGROW>
    end;
    linkaxes(DATA.Handles.imageaxes(panels),'xy');
else
    linkaxes(DATA.Handles.imageaxes(panel),'off');
end;

%--- Update zoom state
viewtype = upper(type(1:3));
if isempty(SET(no).(viewtype).ZoomState)
    if ~strcmp(viewtype,'GLA')
        SET(no).(viewtype).ZoomState = [0.5;ysz+0.5;0.5;xsz+0.5];
    else
        updateglazoomstate(no,ysz);
    end
end;
set(DATA.Handles.imageaxes(panel),...
    'xlim',SET(no).(viewtype).ZoomState(1:2),...
    'ylim',SET(no).(viewtype).ZoomState(3:4));
set(DATA.Handles.imageaxes(panel),'dataaspectratio',...
    [1/yres ...
    1/xres 1]);

%Set up ROI handles do not plot
DATA.Handles.roitext{panel} = zeros(1,SET(no).RoiN);
DATA.Handles.roicontour{panel} = zeros(1,SET(no).RoiN);
hold(DATA.Handles.imageaxes(panel),'on');
DATA.Handles.roicontour{panel}(:) = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
DATA.Handles.roitext{panel}(:) = text(DATA.Handles.imageaxes(panel),0,0,'');
hold(DATA.Handles.imageaxes(panel),'off');

set(DATA.imagefig,'WindowButtonUpFcn',...
    sprintf('%s(''buttonup_Callback'')','segment'));
updatevisibility;

%---------------------------------------
function drawimagetypetext(no,pno,panel)
%---------------------------------------
global DATA SET
%--- Write text describing image type
x = SET(no).YSize*0.9;
y = SET(no).XSize*0.9;

%Image type
DATA.Handles.imagetypetext(panel) = text(...
    x,...
    y,...
    strcat(SET(pno).ImageType,',',SET(pno).ImageViewPlane),...
    'interpreter','none',...
    'parent',DATA.Handles.imageaxes(panel));

extent = get(DATA.Handles.imagetypetext(panel),'extent');
set(DATA.Handles.imagetypetext(panel),...
    'color','white','position',...
    [size(DATA.ViewIM{panel},2)-extent(3) ...
    size(DATA.ViewIM{panel},1)-extent(4) 0]);

%seriesdescription
DATA.Handles.seriesdescriptiontext(panel) = text(...
    x,...
    y,...
    SET(pno).SeriesDescription,...
    'interpreter','none',...
    'parent',DATA.Handles.imageaxes(panel));
extent = get(DATA.Handles.seriesdescriptiontext(panel),'extent');
set(DATA.Handles.seriesdescriptiontext(panel),...
    'color','white','position',...
    [size(DATA.ViewIM{panel},2)-extent(3) ...
    size(DATA.ViewIM{panel},1)-3*extent(4) 0]);

%dicomimagetype
DATA.Handles.dicomimagetypetext(panel) = text(...
    x,...
    y,...
    SET(pno).DICOMImageType,'interpreter','none');
extent = get(DATA.Handles.dicomimagetypetext(panel),'extent');
set(DATA.Handles.dicomimagetypetext(panel),...
    'interpreter','none',...
    'color','white','position',...
    [size(DATA.ViewIM{panel},2)-extent(3) ...
    size(DATA.ViewIM{panel},1)-2*extent(4) 0]);

%slicetimeframetext
if strcmp(DATA.ViewPanelsType{panel},'hla')
    stri = dprintf('Slice:%02d Time:%03d ms',SET(no).HLA.slice,round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
elseif strcmp(DATA.ViewPanelsType{panel},'vla')
    stri = dprintf('Slice:%02d Time:%03d ms',SET(no).VLA.slice,round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
elseif ismember(DATA.ViewPanelsType{panel},{'orthomip','hlamip','vlamip'})
    stri = dprintf('Time:%03d ms',round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
else
    stri = dprintf('Slice:%02d Time:%03d ms',SET(no).CurrentSlice,round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
end
DATA.Handles.slicetimetext(panel) = text(...
    x,...
    y,...
    stri,'interpreter','none');

viewupdatetextposition(panel);

%------------------------------
function drawcontours(no,panel)
%------------------------------
%Initiate handles and draw endo- and epicardial contours of the LV and RV.
global DATA SET
%persistent endohandle
if DATA.Pref.BlackWhite
    endocolor = 'w-';
    epicolor = 'w-';
    rvendocolor = 'w-';
    rvepicolor = 'w-';
else
    endocolor = 'r-';
    epicolor = 'g-';
    rvendocolor = 'm-';
    rvepicolor = 'c-';
end
notlax = ~ismember(DATA.ViewPanelsType{panel}, ...
    {'hla','vla','orthomip','hlamip','vlamip'});

%These are not always properly disposed of. But when?
%Endo
if ~isempty(SET(no).EndoX) && notlax
    DATA.Handles.endocontour(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),...
        SET(no).EndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),endocolor);
    %--- Sector grid
    DATA.drawsectorgrid(no,panel);
else
    DATA.Handles.endocontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,endocolor);
end;%endohandle = DATA.Handles.endocontour;

%Epi
if ~isempty(SET(no).EpiX) && notlax
    DATA.Handles.epicontour(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),...
        SET(no).EpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),epicolor);
else
    DATA.Handles.epicontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,epicolor);%,'linewidth',2);
end;

%RVEndo
if ~isempty(SET(no).RVEndoX) && notlax
    DATA.Handles.rvendocontour(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),...
        SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),rvendocolor);
else
    DATA.Handles.rvendocontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,rvendocolor);
end;

%RVEpi
if ~isempty(SET(no).RVEpiX) && notlax
    DATA.Handles.rvepicontour(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),...
        SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),rvepicolor);
else
    DATA.Handles.rvepicontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,rvepicolor);
end;

%--------------------------------
function updatecontours(no,panel)
%--------------------------------
%Update coordinates of contour handles of one view panels
global DATA SET

if ~isempty(SET(no).EndoX)
    set(DATA.Handles.endocontour(panel),...
        'xdata',SET(no).EndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),...
        'ydata',SET(no).EndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice));
    
    %Sector grid
    DATA.drawsectorgrid(no,panel);
else
    set(DATA.Handles.endocontour(panel),...
        'xdata',[],'ydata',[]);
end;
if ~isempty(SET(no).RVEndoX)
    set(DATA.Handles.rvendocontour(panel),...
        'xdata',SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),...
        'ydata',SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice));
else
    set(DATA.Handles.rvendocontour(panel),...
        'xdata',[],'ydata',[]);
end;
if ~isempty(SET(no).RVEpiX)
    set(DATA.Handles.rvepicontour(panel),...
        'xdata',SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),...
        'ydata',SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice));
else
    set(DATA.Handles.rvepicontour(panel),...
        'xdata',[],'ydata',[]);
end;
if ~isempty(SET(no).EpiX)
    set(DATA.Handles.epicontour(panel),...
        'xdata',SET(no).EpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),...
        'ydata',SET(no).EpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice));
else
    set(DATA.Handles.epicontour(panel),...
        'xdata',[],'ydata',[]);
end;


%-------------------------------------
function drawmontagecontours(no,panel)
%-------------------------------------
%Initiate handles and draw contours for montage view.
global DATA SET

if DATA.Pref.BlackWhite
    endocolor = 'w-';
    epicolor = 'w-';
    rvendocolor = 'w-';
    rvepicolor = 'w-';
else
    endocolor = 'r-';
    epicolor = 'g-';
    rvendocolor = 'm-';
    rvepicolor = 'c-';
end

%Endo
if ~isempty(SET(no).EndoX)
    DATA.Handles.endocontour(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EndoYView(:,SET(no).CurrentTimeFrame),...
        SET(no).EndoXView(:,SET(no).CurrentTimeFrame),endocolor);
else
    DATA.Handles.endocontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,endocolor);
end;
%Epi
if ~isempty(SET(no).EpiX)
    DATA.Handles.epicontour(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EpiYView(:,SET(no).CurrentTimeFrame),...
        SET(no).EpiXView(:,SET(no).CurrentTimeFrame),epicolor);
else
    DATA.Handles.epicontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,epicolor);
end;
%RVEndo
if ~isempty(SET(no).RVEndoX)
    DATA.Handles.rvendocontour(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEndoYView(:,SET(no).CurrentTimeFrame),...
        SET(no).RVEndoXView(:,SET(no).CurrentTimeFrame),rvendocolor);
else
    DATA.Handles.rvendocontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,rvendocolor);
end;
%RVEpi
if ~isempty(SET(no).RVEpiX)
    DATA.Handles.rvepicontour(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEpiYView(:,SET(no).CurrentTimeFrame),...
        SET(no).RVEpiXView(:,SET(no).CurrentTimeFrame),rvepicolor);
else
    DATA.Handles.rvepicontour(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,rvepicolor);
end;

%---------------------------------------
function updatemontagecontours(no,panel)
%---------------------------------------
%Update coordinates of contour handles for montage image view
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

if ~isempty(SET(no).EndoX)
    set(DATA.Handles.endocontour(panel),...
        'xdata',SET(no).EndoYView(:,SET(no).CurrentTimeFrame),...
        'ydata',SET(no).EndoXView(:,SET(no).CurrentTimeFrame));
else
    set(DATA.Handles.endocontour(panel),...
        'xdata',[],'ydata',[]);
end;
if ~isempty(SET(no).RVEndoX)
    set(DATA.Handles.rvendocontour(panel),...
        'xdata',SET(no).RVEndoYView(:,SET(no).CurrentTimeFrame),...
        'ydata',SET(no).RVEndoXView(:,SET(no).CurrentTimeFrame));
else
    set(DATA.Handles.rvendocontour(panel),...
        'xdata',[],'ydata',[]);
end;
if ~isempty(SET(no).RVEpiX)
    set(DATA.Handles.rvepicontour(panel),...
        'xdata',SET(no).RVEpiYView(:,SET(no).CurrentTimeFrame),...
        'ydata',SET(no).RVEpiXView(:,SET(no).CurrentTimeFrame));
else
    set(DATA.Handles.rvepicontour(panel),...
        'xdata',[],'ydata',[]);
end;
if ~isempty(SET(no).EpiX)
    set(DATA.Handles.epicontour(panel),...
        'xdata',SET(no).EpiYView(:,SET(no).CurrentTimeFrame),...
        'ydata',SET(no).EpiXView(:,SET(no).CurrentTimeFrame));
else
    set(DATA.Handles.epicontour(panel),...
        'xdata',[],'ydata',[]);
end;

%------------------------------------
function updatesax3contours(no,panel)
%------------------------------------
%Update coordinates of contour handles for SAX3 image view
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

if ~isempty(SET(no).EndoX)
    set(DATA.Handles.endocontour(panel),...
        'xdata',SET(no).SAX3.EndoY(:,SET(no).CurrentTimeFrame),...
        'ydata',SET(no).SAX3.EndoX(:,SET(no).CurrentTimeFrame));
else
    set(DATA.Handles.endocontour(panel),...
        'xdata',[],'ydata',[]);
end;
if ~isempty(SET(no).EpiX)
    set(DATA.Handles.epicontour(panel),...
        'xdata',SET(no).SAX3.EpiY(:,SET(no).CurrentTimeFrame),...
        'ydata',SET(no).SAX3.EpiX(:,SET(no).CurrentTimeFrame));
else
    set(DATA.Handles.epicontour(panel),...
        'xdata',[],'ydata',[]);
end;

%--------------------------------------------
function updatelongaxiscontours(arg,no,panel)
%--------------------------------------------
%Update coordinates of contour handles for HLA or VLA image view
global DATA SET

if nargin < 3
    panel = find(DATA.ViewPanels == no);
end

if ~isempty(SET(no).EndoX)
    set(DATA.Handles.endocontour(panel),...
        'xdata',SET(no).(arg).EndoY{SET(no).CurrentTimeFrame},...
        'ydata',SET(no).(arg).EndoX{SET(no).CurrentTimeFrame});
else
    set(DATA.Handles.endocontour(panel),...
        'xdata',[],'ydata',[]);
end;
set(DATA.Handles.epicontour(panel),...
    'xdata',[],'ydata',[]);

%-----------------------------------
function drawinterp(no,panel,olshow)
%-----------------------------------
%Initiate handles and draw interpolation points.
%Do not work well with overlays
global DATA SET

if olshow || (isempty(SET(no).EndoInterpX)||isempty(SET(no).EndoInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice}))
    DATA.Handles.endointerp(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,...
        'rs','MarkerFaceColor','r','MarkerSize',4);
else
    DATA.Handles.endointerp(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EndoInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        SET(no).EndoInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        'rs','MarkerFaceColor','r','MarkerSize',4);
end;
if olshow || (isempty(SET(no).EpiInterpX)||isempty(SET(no).EpiInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice}))
    DATA.Handles.epiinterp(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,...
        'gs','MarkerFaceColor','g','MarkerSize',4);
else
    DATA.Handles.epiinterp(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EpiInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        SET(no).EpiInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        'gs','MarkerFaceColor','g','MarkerSize',4);
end;
if olshow || (isempty(SET(no).RVEndoInterpX)||isempty(SET(no).RVEndoInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice}))
    DATA.Handles.rvendointerp(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,...
        'ms','MarkerFaceColor','m','MarkerSize',4);
else
    DATA.Handles.rvendointerp(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEndoInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        SET(no).RVEndoInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        'ms','MarkerFaceColor','m','MarkerSize',4);
end;
if olshow || (isempty(SET(no).RVEpiInterpX)||isempty(SET(no).RVEpiInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice}))
    DATA.Handles.rvepiinterp(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,...
        'cs','MarkerFaceColor','c','MarkerSize',4);
else
    DATA.Handles.rvepiinterp(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEpiInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        SET(no).RVEpiInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        'cs','MarkerFaceColor','c','MarkerSize',4);
end;

%------------------------------
function updateinterp(no,panel)
%------------------------------
%Update coordinates of interpolation point handles
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

if ~isempty(SET(no).EndoInterpX)
  if size(SET(no).EndoInterpX,1)~=SET(no).TSize
    tmpX = cell(SET(no).TSize,SET(no).ZSize);
    tmpY = tmpX;
    for i=1:size(SET(no).EndoInterpX,1)
      for j=1:size(SET(no).EndoInterpX,2)
        tmpX{i,j}=SET(no).EndoInterpX{i,j};
        tmpY{i,j}=SET(no).EndoInterpY{i,j};
      end
    end
    SET(no).EndoInterpX=tmpX;
    SET(no).EndoInterpY=tmpY;
  end
  
  set(DATA.Handles.endointerp(panel),...
    'xdata',SET(no).EndoInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
    'ydata',SET(no).EndoInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice});
else
  set(DATA.Handles.endointerp(panel),...
    'xdata',[],...
    'ydata',[]);
end;
if ~isempty(SET(no).EpiInterpX)
  if size(SET(no).EpiInterpX,1)~=SET(no).TSize
    tmpX = cell(SET(no).TSize,SET(no).ZSize);
    tmpY = tmpX;
    for i=1:size(SET(no).EpiInterpX,1)
      for j=1:size(SET(no).EpiInterpX,2)
        tmpX{i,j}=SET(no).EpiInterpX{i,j};
        tmpY{i,j}=SET(no).EpiInterpY{i,j};
      end
    end
    SET(no).EpiInterpX=tmpX;
    SET(no).EpiInterpY=tmpY;
  end
  
  set(DATA.Handles.epiinterp(panel),...
    'xdata',SET(no).EpiInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
    'ydata',SET(no).EpiInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice});
else
  set(DATA.Handles.epiinterp(panel),...
    'xdata',[],...
    'ydata',[]);
end;
if ~isempty(SET(no).RVEndoInterpX)
  if size(SET(no).RVEndoInterpX,1)~=SET(no).TSize
    tmpX = cell(SET(no).TSize,SET(no).ZSize);
    tmpY = tmpX;
    for i=1:size(SET(no).RVEndoInterpX,1)
      for j=1:size(SET(no).RVEndoInterpX,2)
        tmpX{i,j}=SET(no).RVEndoInterpX{i,j};
        tmpY{i,j}=SET(no).RVEndoInterpY{i,j};
      end
    end
    SET(no).RVEndoInterpX=tmpX;
    SET(no).RVEndoInterpY=tmpY;
  end
  
  set(DATA.Handles.rvendointerp(panel),...
    'xdata',SET(no).RVEndoInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
    'ydata',SET(no).RVEndoInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice});
else
  set(DATA.Handles.rvendointerp(panel),...
    'xdata',[],...
    'ydata',[]);
end;
if ~isempty(SET(no).RVEpiInterpX)
  if size(SET(no).RVEpiInterpX,1)~=SET(no).TSize
    tmpX = cell(SET(no).TSize,SET(no).ZSize);
    tmpY = tmpX;
    for i=1:size(SET(no).RVEpiInterpX,1)
      for j=1:size(SET(no).RVEpiInterpX,2)
        tmpX{i,j}=SET(no).RVEpiInterpX{i,j};
        tmpY{i,j}=SET(no).RVEpiInterpY{i,j};
      end
    end
    SET(no).RVEpiInterpX=tmpX;
    SET(no).RVEpiInterpY=tmpY;
  end
  
  set(DATA.Handles.rvepiinterp(panel),...
    'xdata',SET(no).RVEpiInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
    'ydata',SET(no).RVEpiInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice});
else
  set(DATA.Handles.rvepiinterp(panel),...
    'xdata',[],...
    'ydata',[]);
end;

%-----------------------------------
function drawmontageinterp(no,panel)
%-----------------------------------
%Initiate handles and draw interpolation points for montage view
global DATA SET

if isempty(SET(no).EndoInterpX)||isempty(SET(no).EndoInterpXView{SET(no).CurrentTimeFrame})
    DATA.Handles.endointerp(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,...
        'rs','MarkerFaceColor','r','MarkerSize',3);
else
    DATA.Handles.endointerp(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EndoInterpYView{SET(no).CurrentTimeFrame},...
        SET(no).EndoInterpXView{SET(no).CurrentTimeFrame},...
        'rs','MarkerFaceColor','r','MarkerSize',3);
end;
if isempty(SET(no).EpiInterpX)||isempty(SET(no).EpiInterpXView{SET(no).CurrentTimeFrame})
    DATA.Handles.epiinterp(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,...
        'gs','MarkerFaceColor','g','MarkerSize',3);
else
    DATA.Handles.epiinterp(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EpiInterpYView{SET(no).CurrentTimeFrame},...
        SET(no).EpiInterpXView{SET(no).CurrentTimeFrame},...
        'gs','MarkerFaceColor','g','MarkerSize',3);
end;
if isempty(SET(no).RVEndoInterpX)||isempty(SET(no).RVEndoInterpXView{SET(no).CurrentTimeFrame})
    DATA.Handles.rvendointerp(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,...
        'ms','MarkerFaceColor','m','MarkerSize',3);
else
    DATA.Handles.rvendointerp(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEndoInterpYView{SET(no).CurrentTimeFrame},...
        SET(no).RVEndoInterpXView{SET(no).CurrentTimeFrame},...
        'ms','MarkerFaceColor','m','MarkerSize',3);
end;
if isempty(SET(no).RVEpiInterpX)||isempty(SET(no).RVEpiInterpXView{SET(no).CurrentTimeFrame})
    DATA.Handles.rvepiinterp(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,...
        'cs','MarkerFaceColor','c','MarkerSize',3);
else
    DATA.Handles.rvepiinterp(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEpiInterpYView{SET(no).CurrentTimeFrame},...
        SET(no).RVEpiInterpXView{SET(no).CurrentTimeFrame},...
        'cs','MarkerFaceColor','c','MarkerSize',3);
end;



%-------------------------------------
function updatemontageinterp(no,panel)
%-------------------------------------
%Update coordinates of interpolation point handles for montage view.
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

if ~isempty(SET(no).EndoInterpX)
    set(DATA.Handles.endointerp(panel),...
        'xdata',SET(no).EndoInterpYView{SET(no).CurrentTimeFrame},...
        'ydata',SET(no).EndoInterpXView{SET(no).CurrentTimeFrame});
else
    set(DATA.Handles.endointerp(panel),...
        'xdata',[],...
        'ydata',[]);
end;
if ~isempty(SET(no).EpiInterpX)
    set(DATA.Handles.epiinterp(panel),...
        'xdata',SET(no).EpiInterpYView{SET(no).CurrentTimeFrame},...
        'ydata',SET(no).EpiInterpXView{SET(no).CurrentTimeFrame});
else
    set(DATA.Handles.epiinterp(panel),...
        'xdata',[],...
        'ydata',[]);
end;
if ~isempty(SET(no).RVEndoInterpX)
    set(DATA.Handles.rvendointerp(panel),...
        'xdata',SET(no).RVEndoInterpYView{SET(no).CurrentTimeFrame},...
        'ydata',SET(no).RVEndoInterpXView{SET(no).CurrentTimeFrame});
else
    set(DATA.Handles.rvendointerp(panel),...
        'xdata',[],...
        'ydata',[]);
end;
if ~isempty(SET(no).RVEpiInterpX)
    set(DATA.Handles.rvepiinterp(panel),...
        'xdata',SET(no).RVEpiInterpYView{SET(no).CurrentTimeFrame},...
        'ydata',SET(no).RVEpiInterpXView{SET(no).CurrentTimeFrame});
else
    set(DATA.Handles.rvepiinterp(panel),...
        'xdata',[],...
        'ydata',[]);
end;

%--------------------------
function drawpins(no,panel)
%--------------------------
%Initiate handles and draw pins.
global DATA SET

rpinstri = 'r+';
gpinstri = 'g+';
mpinstri = 'm+';
cpinstri = 'c+';

if isempty(SET(no).EndoPinX)||isempty(SET(no).EndoPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice})
    DATA.Handles.endopin(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,rpinstri);
else
    DATA.Handles.endopin(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EndoPinY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        SET(no).EndoPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},rpinstri);
end;
if isempty(SET(no).EpiPinX)||isempty(SET(no).EpiPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice})
    DATA.Handles.epipin(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,gpinstri);
else
    DATA.Handles.epipin(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EpiPinY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        SET(no).EpiPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},gpinstri);
end;
if isempty(SET(no).RVEndoPinX)||isempty(SET(no).RVEndoPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice})
    DATA.Handles.rvendopin(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,mpinstri);
else
    DATA.Handles.rvendopin(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEndoPinY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        SET(no).RVEndoPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},mpinstri);
end;
if isempty(SET(no).RVEpiPinX)||isempty(SET(no).RVEpiPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice})
    DATA.Handles.rvepipin(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,cpinstri);
else
    DATA.Handles.rvepipin(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEpiPinY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        SET(no).RVEpiPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},cpinstri);
end;

set(DATA.Handles.endopin(panel),'ButtonDownFcn',...
    sprintf('%s(''clickedpin_Callback'',''endo'')','segment'));
set(DATA.Handles.epipin(panel),'ButtonDownFcn',...
    sprintf('%s(''clickedpin_Callback'',''epi'')','segment'));
set(DATA.Handles.rvendopin(panel),'ButtonDownFcn',...
    sprintf('%s(''clickedpin_Callback'',''rvendo'')','segment'));
set(DATA.Handles.rvepipin(panel),'ButtonDownFcn',...
    sprintf('%s(''clickedpin_Callback'',''rvepi'')','segment'));

%----------------------------
function updatepins(no,panel)
%----------------------------
%Update coordinates of pin handles.
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

if ~isempty(SET(no).EndoPinX)
    set(DATA.Handles.endopin(panel),...
        'xdata',SET(no).EndoPinY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        'ydata',SET(no).EndoPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice});
else
    set(DATA.Handles.endopin(panel),...
        'xdata',[],...
        'ydata',[]);
end;
if ~isempty(SET(no).EpiPinX)
    set(DATA.Handles.epipin(panel),...
        'xdata',SET(no).EpiPinY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        'ydata',SET(no).EpiPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice});
else
    set(DATA.Handles.epipin(panel),...
        'xdata',[],...
        'ydata',[]);
end;
if ~isempty(SET(no).RVEndoPinX)
    set(DATA.Handles.rvendopin(panel),...
        'xdata',SET(no).RVEndoPinY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        'ydata',SET(no).RVEndoPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice});
else
    set(DATA.Handles.rvendopin(panel),...
        'xdata',[],...
        'ydata',[]);
end;
if ~isempty(SET(no).RVEpiPinX)
    set(DATA.Handles.rvepipin(panel),...
        'xdata',SET(no).RVEpiPinY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice},...
        'ydata',SET(no).RVEpiPinX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice});
else
    set(DATA.Handles.rvepipin(panel),...
        'xdata',[],...
        'ydata',[]);
end;

%---------------------------------
function drawmontagepins(no,panel)
%---------------------------------
%Initiate handles and draw pins for montage view
global DATA SET

%Color of pins
rpinstri = 'r+';
gpinstri = 'g+';
mpinstri = 'm+';
cpinstri = 'c+';

if isempty(SET(no).EndoPinX)||isempty(SET(no).EndoPinXView{SET(no).CurrentTimeFrame})
    DATA.Handles.endopin(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,rpinstri);
else
    DATA.Handles.endopin(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EndoPinYView{SET(no).CurrentTimeFrame},...
        SET(no).EndoPinXView{SET(no).CurrentTimeFrame},rpinstri);
end;
if isempty(SET(no).EpiPinX)||isempty(SET(no).EpiPinXView{SET(no).CurrentTimeFrame})
    DATA.Handles.epipin(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,gpinstri);
else
    DATA.Handles.epipin(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).EpiPinYView{SET(no).CurrentTimeFrame},...
        SET(no).EpiPinXView{SET(no).CurrentTimeFrame},gpinstri);
end;
if isempty(SET(no).RVEndoPinX)||isempty(SET(no).RVEndoPinXView{SET(no).CurrentTimeFrame})
    DATA.Handles.rvendopin(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,mpinstri);
else
    DATA.Handles.rvendopin(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEndoPinYView{SET(no).CurrentTimeFrame},...
        SET(no).RVEndoPinXView{SET(no).CurrentTimeFrame},mpinstri);
end;
if isempty(SET(no).RVEpiPinX)||isempty(SET(no).RVEpiPinXView{SET(no).CurrentTimeFrame})
    DATA.Handles.rvepipin(panel) = plot(DATA.Handles.imageaxes(panel),NaN,NaN,cpinstri);
else
    DATA.Handles.rvepipin(panel) = plot(DATA.Handles.imageaxes(panel),...
        SET(no).RVEpiPinYView{SET(no).CurrentTimeFrame},...
        SET(no).RVEpiPinXView{SET(no).CurrentTimeFrame},cpinstri);
end;
set(DATA.Handles.endopin(panel),'ButtonDownFcn',...
    sprintf('%s(''clickedpin_Callback'',''endo'')','segment'));
set(DATA.Handles.epipin(panel),'ButtonDownFcn',...
    sprintf('%s(''clickedpin_Callback'',''epi'')','segment'));
set(DATA.Handles.endopin(panel),'ButtonDownFcn',...
    sprintf('%s(''clickedpin_Callback'',''rvendo'')','segment'));
set(DATA.Handles.epipin(panel),'ButtonDownFcn',...
    sprintf('%s(''clickedpin_Callback'',''rvepi'')','segment'));



%-----------------------------------
function updatemontagepins(no,panel)
%-----------------------------------
%Update coordinates of pin handles for montage view
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

if ~isempty(SET(no).EndoPinX)
    set(DATA.Handles.endopin(panel),...
        'xdata',SET(no).EndoPinYView{SET(no).CurrentTimeFrame},...
        'ydata',SET(no).EndoPinXView{SET(no).CurrentTimeFrame});
else
    set(DATA.Handles.endopin(panel),...
        'xdata',[],...
        'ydata',[]);
end;
if ~isempty(SET(no).EpiPinX)
    set(DATA.Handles.epipin(panel),...
        'xdata',SET(no).EpiPinYView{SET(no).CurrentTimeFrame},...
        'ydata',SET(no).EpiPinXView{SET(no).CurrentTimeFrame});
else
    set(DATA.Handles.epipin(panel),...
        'xdata',[],...
        'ydata',[]);
end;
if ~isempty(SET(no).RVEndoPinX)
    set(DATA.Handles.rvendopin(panel),...
        'xdata',SET(no).RVEndoPinYView{SET(no).CurrentTimeFrame},...
        'ydata',SET(no).RVEndoPinXView{SET(no).CurrentTimeFrame});
else
    set(DATA.Handles.rvendopin(panel),...
        'xdata',[],...
        'ydata',[]);
end;
if ~isempty(SET(no).RVEpiPinX)
    set(DATA.Handles.rvepipin(panel),...
        'xdata',SET(no).RVEpiPinYView{SET(no).CurrentTimeFrame},...
        'ydata',SET(no).RVEpiPinXView{SET(no).CurrentTimeFrame});
else
    set(DATA.Handles.rvepipin(panel),...
        'xdata',[],...
        'ydata',[]);
end;

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
        end;
    end;
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
    end;
    DATA.Handles.endointersectionpoints{panel} = plot(DATA.Handles.imageaxes(panel),...
        NaN,NaN,'r.');
    DATA.Handles.epiintersection{panel} = plot(DATA.Handles.imageaxes(panel),...
        NaN,NaN,'g.');
    DATA.Handles.rvendointersection{panel} = plot(DATA.Handles.imageaxes(panel),...
        NaN,NaN,'m.');
    DATA.Handles.rvepiintersection{panel} = plot(DATA.Handles.imageaxes(panel),...
        NaN,NaN,'c.');
end;

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
end;


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
    end;
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
end;

%----------------------------
function updaterois(no,panel)
%----------------------------
%Update coordinates of ROI handles.
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

for loop=1:SET(no).RoiN
    if ismember(SET(no).CurrentSlice,SET(no).Roi(loop).Z) && ...
            ismember(SET(no).CurrentTimeFrame,SET(no).Roi(loop).T)
        if SET(no).Roi(loop).Sign > 0
            roisign = '';
        else
            roisign = ' (-)';
        end
        set(cellref(DATA.Handles.roicontour(panel),loop),...
            'xdata',SET(no).Roi(loop).Y(:,SET(no).CurrentTimeFrame),...
            'ydata',SET(no).Roi(loop).X(:,SET(no).CurrentTimeFrame));
        [ymin,ix] = min(SET(no).Roi(loop).Y(:,SET(no).CurrentTimeFrame));
        set(cellref(DATA.Handles.roitext(panel),loop),...
            'position',...
            [ymin-1 SET(no).Roi(loop).X(ix,SET(no).CurrentTimeFrame)], ...
            'String',{sprintf('%s%s',SET(no).Roi(loop).Name,roisign), ...
            sprintf('%3.1f [cm^2]', ...
            SET(no).Roi(loop).Area(SET(no).CurrentTimeFrame)), ...
            sprintf('%3.1f ± %3.1f', ...
            SET(no).Roi(loop).Mean(SET(no).CurrentTimeFrame), ...
            SET(no).Roi(loop).StD(SET(no).CurrentTimeFrame))},...
            'HorizontalAlignment','right','VerticalAlignment','middle');
        if ismember(loop,SET(no).RoiCurrent)
            set(cellref(DATA.Handles.roicontour(panel),loop), ...
                'LineWidth',DATA.Pref.LineWidth+1);
        else
            set(cellref(DATA.Handles.roicontour(panel),loop), ...
                'LineWidth',DATA.Pref.LineWidth);
        end
    else
        set(cellref(DATA.Handles.roicontour(panel),loop), ...
            'XData',nan,'YData',nan);
        set(cellref(DATA.Handles.roitext(panel),loop), ...
            'Position',[nan nan]);
    end
end;


%-----------------------------------
function updatemontagerois(no,panel)
%-----------------------------------
%Update coordnates of ROI handles for montage view.
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

roih = cat(1,DATA.Handles.roicontour{panel});
for loop=1:SET(no).RoiN
    [xofs,yofs] = calcfunctions('calcoffset', ...
        SET(no).Roi(loop).Z,[],no,panel(1));
    set(roih(:,loop),...
        'xdata',SET(no).Roi(loop).Y(:,SET(no).CurrentTimeFrame)+yofs,...
        'ydata',SET(no).Roi(loop).X(:,SET(no).CurrentTimeFrame)+xofs);
    if ismember(loop,SET(no).RoiCurrent)
        set(DATA.Handles.roicontour{panel}(loop), ...
            'LineWidth',DATA.Pref.LineWidth+1);
    else
        set(DATA.Handles.roicontour{panel}(loop), ...
            'LineWidth',DATA.Pref.LineWidth);
    end
end;


%------------------------------
function drawmeasures(no,panel)
%------------------------------
%Draw measurements, if available
global DATA SET

if not(isempty(SET(no).Measure))
    DATA.Handles.measureline{panel} = cell(1,length(SET(no).Measure));
    DATA.Handles.measuretext{panel} = zeros(1,length(SET(no).Measure));
    
    [measure,slice] = segment('getmeasurecoords',no,panel);
    for loop=1:length(measure)
        % JU: Measures now time specific
        % NL: Measures can now cross multiple slices
        ziv = round(measure(loop).Z);
        ziv = min(ziv):max(ziv);
        if ismember(slice,ziv) && ...
                ((SET(no).Measure(loop).T==SET(no).CurrentTimeFrame)||(isnan(SET(no).Measure(loop).T)))
            DATA.Handles.measureline{panel}{loop} = plot(DATA.Handles.imageaxes(panel),...
                measure(loop).Y,measure(loop).X,DATA.GUISettings.MeasureLineSpec);
            if ~all(ziv == slice)
                set(DATA.Handles.measureline{panel}{loop},'LineStyle','--');
            end
            set(DATA.Handles.measureline{panel}{loop},'markersize',DATA.GUISettings.MeasureLineMarkerSize);
            [ymax,ix] = max(measure(loop).Y);
            DATA.Handles.measuretext{panel}(loop) = text(...
                'position',[ymax+1 measure(loop).X(ix)],...
                'string',sprintf('%s\n%0.1f [mm]',SET(no).Measure(loop).Name,SET(no).Measure(loop).Length),...
                'parent',DATA.Handles.imageaxes(panel),...
                'color',[1 1 1],...
                'interpreter','none','horizontalalignment','left');
            DATA.measurefontsize(panel,loop);
        else
            DATA.Handles.measureline{panel}{loop} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
            DATA.Handles.measuretext{panel}(loop) = text(...
                'position',[NaN NaN],...
                'string','',...
                'parent',DATA.Handles.imageaxes(panel));
        end;
    end;
    set([DATA.Handles.measureline{panel}{:}; DATA.Handles.measuretext{panel}],'ButtondownFcn',...
        sprintf('segment(''measurepoint_Buttondown'',%d)',panel));
else
    DATA.Handles.measureline{panel} = {};
    DATA.Handles.measuretext{panel} = [];
end;

%--------------------------------
function updatemeasures(no,panel)
%--------------------------------
%Update coordinates of measurement handles
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

if not(isempty(SET(no).Measure))
    [measure,slice] = segment('getmeasurecoords',no,panel);
    for loop=1:length(SET(no).Measure)
        ziv = round(measure(loop).Z);
        ziv = min(ziv):max(ziv);
        if ismember(slice,ziv) && ...
                (SET(no).Measure(loop).T==SET(no).CurrentTimeFrame)||(isnan(SET(no).Measure(loop).T))
            mlh = cellref(DATA.Handles.measureline(panel),loop);
            set([mlh{:}],...
                'xdata', measure(loop).Y,...
                'ydata', measure(loop).X,...
                'Color',DATA.GUISettings.MeasureLineSpec(1),...
                'Marker', DATA.GUISettings.MeasureLineSpec(2),...
                'LineStyle',DATA.GUISettings.MeasureLineSpec(3),...
                'markersize',DATA.GUISettings.MeasureLineMarkerSize);
            if ~all(ziv == slice)
                set([mlh{:}], ...
                    'LineStyle','--');
            end
            [ymax,ix] = max(measure(loop).Y);
            set(cellref(DATA.Handles.measuretext(panel),loop),...
                'position',[ymax+1 measure(loop).X(ix)],...
                'string',sprintf('%s\n%0.1f [mm]',SET(no).Measure(loop).Name,SET(no).Measure(loop).Length),...
                'color',[1 1 1],...
                'interpreter','none',...
                'HorizontalAlignment','left','VerticalAlignment','middle');
            %'parent',DATA.Handles.imageaxes(panel),...
            DATA.measurefontsize(panel,loop);
        else
            mlh = cellref(DATA.Handles.measureline(panel),loop);
            set([mlh{:}],...
                'XData',nan,'YData',nan);
            set(cellref(DATA.Handles.measuretext(panel),loop),...
                'Position',[nan nan]);
        end;
    end;
    
    
end;

%-------------------------------------
function drawmontagemeasures(no,panel)
%-------------------------------------
%Initiate handles and draw measurements for montage view.
global DATA SET

if not(isempty(SET(no).Measure))
    DATA.Handles.measureline{panel} = cell(1,length(SET(no).Measure));
    DATA.Handles.measuretext{panel} = zeros(1,length(SET(no).Measure));
    
    for loop=1:length(SET(no).Measure)
        if (SET(no).Measure(loop).T==SET(no).CurrentTimeFrame)||(isnan(SET(no).Measure(loop).T))
            %offsetx = SET(no).XSize*floor((SET(no).Measure(loop).Z-1)/DATA.ViewPanelsMatrix{panel}(2));
            %offsety = SET(no).YSize*mod(SET(no).Measure(loop).Z-1,DATA.ViewPanelsMatrix{panel}(2));
            ziv = min(SET(no).Measure(loop).Z):max(SET(no).Measure(loop).Z);
            [offsetx,offsety] = calcfunctions('calcoffset', ...
                ziv,[],no,panel);
            [OFFSETY,MEASUREY] = ndgrid(SET(no).Measure(loop).Y,offsety);
            [OFFSETX,MEASUREX] = ndgrid(SET(no).Measure(loop).X,offsetx);
            DATA.Handles.measureline{panel}{loop} = plot(DATA.Handles.imageaxes(panel),...
                OFFSETY+MEASUREY,...
                OFFSETX+MEASUREX,...
                DATA.GUISettings.MeasureLineSpec)';
            if size(OFFSETX,2) > 1
                set(DATA.Handles.measureline{panel}{loop},'LineStyle','--')
            end
            set(DATA.Handles.measureline{panel}{loop},'markersize',DATA.GUISettings.MeasureLineMarkerSize);
            [ymax,ix] = max(SET(no).Measure(loop).Y);
            DATA.Handles.measuretext{panel}(loop) = text(...
                'position',[offsety(end)+ymax+1 offsetx(end)+SET(no).Measure(loop).X(ix)],...
                'string',sprintf('%s\n%0.1f [mm]',SET(no).Measure(loop).Name,SET(no).Measure(loop).Length),...
                'color',[1 1 1],...
                'parent',DATA.Handles.imageaxes(panel),...
                'interpreter','none',...
                'horizontalalignment','center');
            DATA.measurefontsize(panel,loop);
        else
            DATA.Handles.measureline{panel}{loop} = plot(DATA.Handles.imageaxes(panel),NaN,NaN);
            DATA.Handles.measuretext{panel}(loop) = text(...
                'position',[NaN NaN],...
                'string','',...
                'parent',DATA.Handles.imageaxes(panel));
        end;
    end;
    set([DATA.Handles.measureline{panel}{:} DATA.Handles.measuretext{panel}],'ButtondownFcn',...
        sprintf('segment(''measurepoint_Buttondown'',%d)',panel));
    
    
else
    DATA.Handles.measureline{panel} = {};
    DATA.Handles.measuretext{panel} = [];
end;

%---------------------------------------
function updatemontagemeasures(no,panel)
%---------------------------------------
%Update coordinates of measurement handles, for montage view.
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

lineh = cat(1,DATA.Handles.measureline{panel});
texth = cat(1,DATA.Handles.measuretext{panel});
if not(isempty(SET(no).Measure))
    
    for loop=1:length(SET(no).Measure)
        if ((SET(no).Measure(loop).T==SET(no).CurrentTimeFrame)||(isnan(SET(no).Measure(loop).T)))
            ziv = min(SET(no).Measure(loop).Z):max(SET(no).Measure(loop).Z);
            [offsetx,offsety] = calcfunctions('calcoffset', ...
                ziv,[],no,panel(1));
            [OFFSETY,MEASUREY] = ndgrid(SET(no).Measure(loop).Y,offsety);
            [OFFSETX,MEASUREX] = ndgrid(SET(no).Measure(loop).X,offsetx);
            sz = size(OFFSETY);
            set(lineh{:,loop},...
                {'xdata'}, mat2cell(OFFSETY+MEASUREY,sz(1),ones(1,sz(2)))',...
                {'ydata'}, mat2cell(OFFSETX+MEASUREX,sz(1),ones(1,sz(2)))',...
                'Color',DATA.GUISettings.MeasureLineSpec(1),...
                'Marker', DATA.GUISettings.MeasureLineSpec(2),...
                'LineStyle',DATA.GUISettings.MeasureLineSpec(3),...
                'markersize',DATA.GUISettings.MeasureLineMarkerSize);
            if sz(2) > 1
                set(lineh{:,loop},'LineStyle','--')
            end
            [ymax,ix] = max(SET(no).Measure(loop).Y);
            set(texth(:,loop),...
                'position',[offsety(end)+ymax+1 offsetx(end)+SET(no).Measure(loop).X(ix)],...
                'string',sprintf('%s\n%0.1f [mm]',SET(no).Measure(loop).Name,SET(no).Measure(loop).Length),...
                'color',[1 1 1],...
                'interpreter','none','horizontalalignment','center');
            %'parent',DATA.Handles.imageaxes(panel),...
            DATA.measurefontsize(panel,loop);
        else
            set(lineh(:,loop),'XData',nan,'YData',nan);
            set(texth(:,loop),'Position',[nan nan]);
        end;
    end;
    
end;

%--------------------------------------
function drawannotationpoints(no,panel)
%--------------------------------------
%Draw annotation points, if available
global DATA SET

if not(isempty(SET(no).Point))
    
    markerSize=7;
    lineSize=1;
    
    DATA.Handles.pointp{panel} = zeros(1,length(SET(no).Point.X));
    DATA.Handles.pointo{panel} = DATA.Handles.pointp{panel};
    DATA.Handles.pointtext{panel} = DATA.Handles.pointp{panel};
    for loop=1:length(SET(no).Point.X)
        if strcmp(SET(no).Point.Label(loop),'Ant track') || strcmp(SET(no).Point.Label(loop),'Inf track') || strcmp(SET(no).Point.Label(loop),'RV track') || strcmp(SET(no).Point.Label(loop),'Tracked Fwd')
            DATA.Handles.pointp{panel}(loop) = plot(DATA.Handles.imageaxes(panel),...
                SET(no).Point.Y(loop),SET(no).Point.X(loop),'b+','linewidth', lineSize,'markersize', markerSize);
            DATA.Handles.pointo{panel}(loop) = plot(DATA.Handles.imageaxes(panel),...
                SET(no).Point.Y(loop),SET(no).Point.X(loop),'bo','linewidth', lineSize,'markersize', markerSize);
            DATA.Handles.pointtext{panel}(loop) = text(...
                'position',[SET(no).Point.Y(loop)+2 SET(no).Point.X(loop)],...
                'string',SET(no).Point.Label{loop},...
                'parent',DATA.Handles.imageaxes(panel));
        else
            DATA.Handles.pointp{panel}(loop) = plot(DATA.Handles.imageaxes(panel),...
                SET(no).Point.Y(loop),SET(no).Point.X(loop),'w+');
            DATA.Handles.pointo{panel}(loop) = plot(DATA.Handles.imageaxes(panel),...
                SET(no).Point.Y(loop),SET(no).Point.X(loop),'wo');
            DATA.Handles.pointtext{panel}(loop) = text(...
                'position',[SET(no).Point.Y(loop)+2 SET(no).Point.X(loop)],...
                'string',SET(no).Point.Label{loop},...
                'parent',DATA.Handles.imageaxes(panel));
        end
        
        if isnan(SET(no).Point.T(loop))
            set(DATA.Handles.pointtext{panel}(loop),'fontweight','bold')
        end;
        hide = false;
        if not(round(SET(no).Point.Z(loop))==SET(no).CurrentSlice)
            hide = true;
        end;
        if (SET(no).Point.T(loop)~=SET(no).CurrentTimeFrame)&&not(isnan(SET(no).Point.T(loop)))
            hide = true;
        end;
        if hide
            set([DATA.Handles.pointp{panel}(loop) DATA.Handles.pointo{panel}(loop)], ...
                'XData',nan,'YData',nan);
            set(DATA.Handles.pointtext{panel}(loop),'Position',[nan nan]);
        end;
    end;
    set(DATA.Handles.pointtext{panel},'Color',[1 1 1]);
    set([...
        DATA.Handles.pointp{panel} ...
        DATA.Handles.pointo{panel} ...
        DATA.Handles.pointtext{panel}],'ButtonDownFcn',...
        sprintf('%s(''pointat_Buttondown'', %d)','annotationpoint',panel));
    
    
else
    DATA.Handles.pointp{panel} = [];
    DATA.Handles.pointo{panel} = [];
    DATA.Handles.pointtext{panel} = [];
end;

%----------------------------------------
function updateannotationpoints(no,panel)
%----------------------------------------
%Update coordinates of annotation point handles
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

if not(isempty(SET(no).Point))
    if ~isempty(SET(no).Point.Z)
        ind = (round(SET(no).Point.Z)==SET(no).CurrentSlice)&...
            (isnan(SET(no).Point.T)|(SET(no).Point.T==SET(no).CurrentTimeFrame));
        pos = find(ind);
        for loop=1:length(pos)
            set([...
                cellref(DATA.Handles.pointo(panel),pos(loop)) ...
                cellref(DATA.Handles.pointp(panel),pos(loop))],...
                'xdata',SET(no).Point.Y(pos(loop)),...
                'ydata',SET(no).Point.X(pos(loop)));
            set(cellref(DATA.Handles.pointtext(panel),pos(loop)),...
                'position',....
                [SET(no).Point.Y(pos(loop))+2 SET(no).Point.X(pos(loop)) 0]);
        end;
        if any(~ind)
            set([...
                cellref(DATA.Handles.pointo(panel),~ind) ...
                cellref(DATA.Handles.pointp(panel),~ind)], ...
                'XData',nan,'YData',nan);
            set(cellref(DATA.Handles.pointtext(panel),~ind), ...
                'Position',[nan nan]);
        end
    end;
    
    
end;

%-----------------------------------
function drawmontagepoints(no,panel)
%-----------------------------------
%Initiate handles and draw annotation points for montage view.
global DATA SET

if not(isempty(SET(no).Point))
    DATA.Handles.pointp{panel} = zeros(1,length(SET(no).Point.X));
    DATA.Handles.pointo{panel} = DATA.Handles.pointp{panel};
    DATA.Handles.pointtext{panel} = DATA.Handles.pointp{panel};
    for loop=1:length(SET(no).Point.X)
        [offsety,offsetx] = calcfunctions('calcoffset',...
            SET(no).Point.Z(loop),[],no,panel(1));
        DATA.Handles.pointp{panel}(loop) = plot(DATA.Handles.imageaxes(panel),...
            offsetx+SET(no).Point.Y(loop),offsety+SET(no).Point.X(loop),'w+');
        DATA.Handles.pointo{panel}(loop) = plot(DATA.Handles.imageaxes(panel),...
            offsetx+SET(no).Point.Y(loop),offsety+SET(no).Point.X(loop),'wo');
        DATA.Handles.pointtext{panel}(loop) = text(...
            'parent',DATA.Handles.imageaxes(panel),...
            'position',[offsetx+SET(no).Point.Y(loop)+2 offsety+SET(no).Point.X(loop)],...
            'string',SET(no).Point.Label{loop});
        if isnan(SET(no).Point.T(loop))
            set(DATA.Handles.pointtext{panel}(loop),'fontweight','bold')
        end;
        hide = false;
        if (SET(no).Point.T(loop)~=SET(no).CurrentTimeFrame)&&not(isnan(SET(no).Point.T(loop)))
            hide = true;
        end;
        if hide
            set([DATA.Handles.pointp{panel}(loop) DATA.Handles.pointo{panel}(loop)], ...
                'XData',nan,'YData',nan);
            set(DATA.Handles.pointtext{panel}(loop),'Position',[nan nan]);
        end;
    end;
    set(DATA.Handles.pointtext{panel},'Color',[1 1 1]);
    set([...
        DATA.Handles.pointp{panel} ...
        DATA.Handles.pointo{panel} ...
        DATA.Handles.pointtext{panel}],'ButtonDownFcn',...
        sprintf('%s(''pointat_Buttondown'',%d)','annotationpoint',panel));
    
else
    DATA.Handles.pointp{panel} = [];
    DATA.Handles.pointo{panel} = [];
    DATA.Handles.pointtext{panel} = [];
end;

%-------------------------------------
function updatemontagepoints(no,panel)
%-------------------------------------
%Update coordinates of annotation point handles for montage view.
global DATA SET

if nargin < 2
    panel = find(DATA.ViewPanels == no);
end

oh = cat(1,DATA.Handles.pointo{panel});
ph = cat(1,DATA.Handles.pointp{panel});
if not(isempty(SET(no).Point))
    if ~isempty(SET(no).Point.Z)
        ind = (SET(no).Point.T==SET(no).CurrentTimeFrame)|isnan(SET(no).Point.T);
        pos = find(ind);
        for loop=1:length(pos)
            [xofs,yofs] = calcfunctions('calcoffset',...
                SET(no).Point.Z(pos(loop)),[],no,panel(1));
            set([...
                oh(:,pos(loop)) ...
                ph(:,pos(loop))],...
                'xdata',yofs+SET(no).Point.Y(pos(loop)),...
                'ydata',xofs+SET(no).Point.X(pos(loop)));
            set(cellref(DATA.Handles.pointtext(panel),pos(loop)),...
                'position',....
                [SET(no).Point.Y(pos(loop))+2+yofs SET(no).Point.X(pos(loop))+xofs 0]);
        end;
        if any(~ind)
            set([oh(:,~ind) ph(:,~ind)], 'XData',nan,'YData',nan);
            set(cellref(DATA.Handles.pointtext(panel),~ind),'Position',[nan nan]);
        end
    end;
    
end;

%----------------------------------------------
function drawmontageimagetypetext(no,pno,panel)
%----------------------------------------------
%Initiate handles and draw image type text
global DATA SET

DATA.Handles.imagetypetext(panel) = text(...
    size(DATA.ViewIM{panel},2)*0.8,...
    size(DATA.ViewIM{panel},1)*0.95,...
    strcat(SET(pno).ImageType,',',SET(pno).ImageViewPlane),...
    'interpreter','none',...
    'parent',DATA.Handles.imageaxes(panel));

%seriesdesciption
DATA.Handles.seriesdescriptiontext(panel) = text(...
    size(DATA.ViewIM{panel},2)*0.8,...
    size(DATA.ViewIM{panel},1)*0.95,...
    SET(pno).SeriesDescription,...
    'interpreter','none',...
    'parent',DATA.Handles.imageaxes(panel));

%dicomimagetypetext
DATA.Handles.dicomimagetypetext(panel) = text(...
    size(DATA.ViewIM{panel},2)*0.8,...
    size(DATA.ViewIM{panel},1)*0.95,...
    SET(pno).DICOMImageType,...
    'interpreter','none',...
    'parent',DATA.Handles.imageaxes(panel));

%slicetimetext
stri = dprintf('Time:%03d ms',round(1000*SET(no).TimeVector(SET(no).CurrentTimeFrame)));
DATA.Handles.slicetimetext(panel) = text(...
    size(DATA.ViewIM{panel},2)*0.8,...
    size(DATA.ViewIM{panel},1)*0.95,...
    stri,...
    'interpreter','none',...
    'parent',DATA.Handles.imageaxes(panel));
viewupdatetextposition(panel);


%---------------------------
function drawcolorbar(panel)
%---------------------------
%Draw color bar in an image panel
global DATA SET
if nargin < 1
    panel = DATA.CurrentPanel;
end

%Clean up
try
  delete(DATA.Handles.colorbar(panel));
catch %#ok<CTCH>
end

if not(DATA.GUISettings.ShowColorbar) % strcmp(get(DATA.Handles.colorbaricon,'state'),'off')
  DATA.Handles.colorbar(panel) = nan;
  return
end

no = DATA.ViewPanels(panel);
parenth = DATA.Handles.imageaxes(panel);
hold(parenth,'on');
slice = SET(no).CurrentSlice;

% cmap = flipud(SET(no).Colormap);
% if isempty(cmap)
%   cmap = flipud(colormap('gray'));
% end
% im = cat(3,cmap(:,1),cmap(:,2),cmap(:,3));

%Use remapuint8
if isa(SET(no).IM,'single')
  imcur = SET(no).IM(:,:,SET(no).CurrentTimeFrame,slice);
  barminmax = quantile(imcur(:),[0.1 0.9]);
  im = calcfunctions('remapuint8',linspace( ...
    barminmax(2),barminmax(1),DATA.GUISettings.ColorMapSize)');
else
  imcur = single(SET(no).IM(:,:,SET(no).CurrentTimeFrame,slice));
  barminmax = [0 quantile(imcur(:),0.9)];
  im = linspace(barminmax(2),barminmax(1),DATA.GUISettings.ColorMapSize)';
end

ax = axis(parenth);
DATA.Handles.colorbar(panel) = image(ax(2)-[(ax(2)-ax(1))/50 0],ax(3:4),im,...
  'parent',parenth);
minmax = calcfunctions('calctruedata',barminmax,no);

rvec = [10 50 100 500 1000];
for r = rvec;
  scalevec = r*round(minmax(2)/r):-r:minmax(1);
  if length(scalevec) <= 12
    break
  end
end
scalevec = num2cell(scalevec);
strifun = @(x)sprintf('%0.0f',x);
scalevec = cellfun(strifun,scalevec,'UniformOutput',false);
tickvec = linspace(ax(3),ax(4),numel(scalevec));
set(parenth,...
  'ytick',tickvec,...
  ... %'yticklabel',scalevec, ...
  'ycolor','w','yaxislocation','right','visible','on');
texth = zeros(numel(tickvec),1);
for i = 1:numel(tickvec)
  texth(i) = text(ax(2)-3,tickvec(i),scalevec{i});
end
set(texth,'Color','w','HorizontalAlignment','right');

%-------------------------------------
function drawviabilityhelper(no,panel)
%-------------------------------------
%Function to draw viability contour on screeen used from drawimageslice,
%drawimagemontage

global DATA SET

%Get mask
switch DATA.ViewPanelsType{panel}
    case {'montage','montagerow','montagefit','sax3'}
        result = segment('reshape2layout',SET(no).Scar.Result,no,panel);
    case {'one','ortho'}
        result = SET(no).Scar.Result(:,:,SET(no).CurrentSlice);
    otherwise
        myfailed('Unknown viewtype in drawviabilityhelper.');
end;

if sum(result(:))>0
    
    %Ensure that we delete stuff
    set(DATA.Handles.imageaxes(panel),'NextPlot','add');
    
    %Normal scar outline
    [~,DATA.Handles.scarcontour{panel}] = contour(DATA.Handles.imageaxes(panel),double(result),[0.5 0.5]);
    if DATA.Pref.BlackWhite
        set(DATA.Handles.scarcontour{panel},'linecolor',[1 1 1]);
    else
        set(DATA.Handles.scarcontour{panel},'linecolor',[1 1 0]);
    end;
    if DATA.Pref.LineWidth>0
        set(DATA.Handles.scarcontour{panel},'linewidth',DATA.Pref.LineWidth);
    else
        set(DATA.Handles.scarcontour{panel},'visible','off');
    end;
    
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
            %disp(sprintf('noweight:%0.5g weight:%0.5g',sum(result(:)),sum(weighting(:))));
            
            %Find threshold so that it would include f% of the area
            f = 1-f; %Since we want the pixelse that are larger than the threshold
            thres = sortedweighting(min(max(round(length(sortedweighting)*f),1),length(sortedweighting)));
            
            %Convert to right size
            switch DATA.ViewPanelsType{panel}
                case {'montage','montagerow','montagefit','sax3'}
                    weighting = segment('reshape2layout',temp,no,panel);
                case {'one','ortho'}
                    weighting = temp(:,:,SET(no).CurrentSlice);
                otherwise
                    myfailed('Unknown viewtype in drawviabilityhelper.');
            end;
            
            [~,DATA.Handles.weightedscarcontour{panel}] = contour(DATA.Handles.imageaxes(panel),weighting,[thres thres]);
            if DATA.Pref.BlackWhite
                set(DATA.Handles.weightedscarcontour{panel},'linecolor',[1 1 1]);
            else
                set(DATA.Handles.weightedscarcontour{panel},'linecolor',[1 0.5 0.5]);
            end;
            if DATA.Pref.LineWidth>0
                set(DATA.Handles.weightedscarcontour{panel},'linewidth',DATA.Pref.LineWidth);
            else
                set(DATA.Handles.weightedscarcontour{panel},'visible','off');
            end;
        else
            DATA.Handles.weightedscarcontour{panel} = [];
        end; %nonempty weighting
    else
        DATA.Handles.weightedscarcontour{panel} = [];
    end %use weighting
else
    %No scar to display
    DATA.Handles.scarcontour{panel} = [];
    DATA.Handles.weightedscarcontour{panel} = [];
end;

%Update taken mocontour and moextentcontour graphically
if sum(SET(no).Scar.NoReflow(:))>0
    
    if ~isfield(SET(no).Scar,'MOThreshold')
        SET(no).Scar.MOThreshold = 1.5;
    end;
    
    %Extract and convert to right size
    %Take image and mask with no reflow.
    doit = true;
    
    switch DATA.ViewPanelsType{panel}
        case {'montage','montagerow','montagefit','sax3'}
            contourim = SET(no).Scar.IM;
            contourim(~SET(no).Scar.NoReflow) = 1; %Mask with no reflow. 1 is max.
            contourim = segment('reshape2layout',contourim,no,panel,NaN); %NaN is outsideelement.
        case {'one','ortho'}
            if existfunctions('anyall',SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice))
                contourim = SET(no).Scar.IM(:,:,SET(no).CurrentSlice);
                contourim(~SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice)) = 1; %1 is max
            else
                doit = false;
            end;
        otherwise
            myfailed('Unknown viewtype in drawviabilityhelper.');
    end;
    
    if (not(isequal(length(SET(no).Scar.mthreshold),SET(no).ZSize)||isequal(length(SET(no).Scar.mthreshold),1)))&& doit
        %Length is not correct.
        viability('viabilitycalc'); %this will updat the size correctly.
    end;
    
    if doit
        %If we need to do, then doit.
        
        %Use remote intensity to calculate good threshold
        if length(SET(no).Scar.mthreshold)<2
            thres=SET(no).Scar.mthreshold;
        else
            thres = SET(no).Scar.mthreshold(SET(no).CurrentSlice);
        end
        thres = thres*SET(no).Scar.MOThreshold; %OBS if you change here, you need ALSO to change in viabilitycalcvolume.
        
        hold(DATA.Handles.imageaxes(panel),'on');
        [~,DATA.Handles.mocontour{panel}] = contour(DATA.Handles.imageaxes(panel),...
            contourim,...
            [thres thres]);
        
        [~,DATA.Handles.moextentcontour{panel}] = contour(DATA.Handles.imageaxes(panel),...
            contourim,...
            [0.99 0.99]);
        
        moh = DATA.Handles.mocontour{panel};
        moextenth = DATA.Handles.moextentcontour{panel};
        
        if DATA.Pref.BlackWhite
            set(moh,'linecolor',[1 1 1],'linestyle',':');
            set(moextenth,'linecolor',[1 1 1],'linestyle','-');
        else
            set(moh,'linecolor',[1 0 0],'linestyle',':');
            set(moextenth,'linecolor',[1 0 0],'linestyle','-');
        end;
        if DATA.Pref.LineWidth>0
            set(moh,'linewidth',DATA.Pref.LineWidth);
            set(moextenth,'linewidth',DATA.Pref.LineWidth);
        else
            set(moh,'visible','off');
            set(moextenth,'visible','off');
        end;
    else
        %No mo in current slice and one view, do nothing
        DATA.Handles.mocontour{panel} = [];
        DATA.Handles.moextentcontour{panel} = [];
    end;
    
else
    DATA.Handles.mocontour{panel} = [];
    DATA.Handles.moextentcontour{panel} = [];
end;

%------------------------------
function showviabilityedits(no)
%------------------------------
%Show viability edits on screen as a temporary overlay.
global DATA SET NO

if nargin==0
    no = NO;
end;

if isempty(SET(no).Scar)
    return;
    %viability('viabilityreset_Callback');
end;

%Check if menu is enabled
if isequal(get(DATA.Handles.hideoverlayicon,'state'),'on')
    manualinteraction = false;
else
    manualinteraction = true;
end;

if isequal(get(DATA.Handles.viabilityshowinfarctaswhitemenu,'checked'),'on')
    showaswhite = true;
else
    showaswhite = false;
end;

if ~isempty(SET(no).Scar.GreyZone.map) && ...
        isequal(get(DATA.Handles.viabilityshowgrayzonemenu,'checked'),'on')
    showgreyzone = true;
    manualinteraction = false;
else
    showgreyzone = false;
end

%If neither then we can just safely exit.
if not(manualinteraction || showaswhite || showgreyzone)
    return;
end;

tempnos=no;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
    return;
end

panelstodo = find(DATA.ViewPanels==no);

for panel=panelstodo
    
    %Determine if it is a RGB image or not.
    if not(DATA.Silent)
        temp = get(DATA.Handles.imagehandle(panel),'CData');
        if ndims(temp)>2
            isrgbimage = true;
        else
            isrgbimage = false;
        end;
    else
        isrgbimage=false;
    end
    
    clear temp;
    
    %Ok lets draw it
    switch DATA.ViewPanelsType{panel}
        case {'one','mmodespatial','ortho'}
            temp = SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
            if showgreyzone
                temp = imresize(temp,size(SET(no).Scar.GreyZone.map(:,:,1)),'bilinear');
            end
            if isrgbimage  %RGB iamge
                colmap = SET(no).Colormap;
                if isempty(colmap)
                    colmap = colormap('gray');
                end
                tempuint8 = min(uint8(255),uint8(1)+calcfunctions('remapuint8',temp,no,colmap));
            else  %Default grayscale image
                tempuint8 = min(uint8(255),uint8(1)+calcfunctions('remapuint8',temp));
                tempuint8 = repmat(tempuint8,[1 1 3]); %EH:
            end;
            clear temp;
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
                    
                    
                end;
            end;
            
            if not(isfield(SET(no).Scar,'NoReflow'))
                SET(no).Scar.NoReflow = SET(no).Scar.Auto;
            end;
            
            if manualinteraction
                tmp=SET(no).Scar.Manual(:,:,SET(no).CurrentSlice);
                tmp=tmp(:);
                scarimrgb((tmp==int8( 1)),2) = uint8(255);
                scarimrgb((tmp==int8(-1)),3) = uint8(255);
                tmp=SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice);
                tmp=tmp(:);
                scarimrgb(logical(tmp),1) = uint8(255);
            end
            scarimrgb=reshape(scarimrgb,[sz(1:2) 3]); %EH: added 3
            
        case {'montage','montagerow','montagefit','sax3'}
            % Convert to 2D and layout
            
            if isrgbimage  %RGB iamge
                colmap = SET(no).Colormap;
                if isempty(colmap)
                    colmap = colormap('gray');
                end
                tempuint8 = calcfunctions('remapuint8',...
                    segment('reshape2layout',squeeze(SET(no).IM(:,:,SET(no).CurrentTimeFrame,:)),no,panel),...
                    no,colmap);
            else  %Default grayscale image
                tempuint8 = calcfunctions('remapuint8',...
                    segment('reshape2layout',squeeze(SET(no).IM(:,:,SET(no).CurrentTimeFrame,:)),no,panel),...
                    no,calcfunctions('returnmapping',no,true));
            end
            tempuint8 = min(uint8(255),uint8(1)+tempuint8);
            sz = size(tempuint8);
            scarimrgb = reshape(tempuint8,[prod(sz(1:2)) 3]);
            
            tmp = segment('reshape2layout',SET(NO).Scar.Result,no,panel);
            tmp = logical(tmp(:));
            
            %Show infarct
            if DATA.Pref.LineWidth>0
                if showaswhite
                    scarimrgb(tmp,:) = uint8(255);
                    %Greyzone
                elseif showgreyzone
                    temp = imresize(SET(no).Scar.GreyZone.map,[SET(no).XSize SET(no).YSize],'bilinear');
                    greytmp = segment('reshape2layout',temp,no,panel);
                    gztmp = greytmp(:) == 1;
                    coretmp = greytmp(:) == 2;
                    
                    scarimrgb(coretmp,:) = repmat(uint8([127 0 0]),sum(coretmp),1);
                    scarimrgb(gztmp(:),:) = repmat(uint8([127 127 0]),sum(gztmp(:)),1);
                end;
            end;
            
            if not(isfield(SET(no).Scar,'NoReflow'))
                SET(no).Scar.NoReflow = repmat(uint8(1),size(SET(no).Scar.Auto));
            end;
            
            if manualinteraction
                tmp=segment('reshape2layout',SET(no).Scar.Manual,no,panel);
                tmp=tmp(:);
                scarimrgb((tmp==int8( 1)),2) = uint8(255);
                scarimrgb((tmp==int8(-1)),3) = uint8(255);
                tmp=segment('reshape2layout',SET(no).Scar.NoReflow,no,panel);
                tmp=tmp(:);
                scarimrgb(logical(tmp),1) = uint8(255);
            end
            scarimrgb=reshape(scarimrgb,sz);
    end
    if not(DATA.Silent)
        set(DATA.Handles.imagehandle(panel),'CData',scarimrgb);
    end
end

%Update menu. Seems unnecessary.
viability('viabilitymenu');

%-------------------------------
function drawmarhelper(no,panel)
%-------------------------------
%Function to draw MaR contour on screeen used from drawimageslice,
%drawimagemontage

global DATA SET

delete(DATA.Handles.marcontour{panel});

%Get mask
switch DATA.ViewPanelsType{panel}
    case {'montage','montagerow','montagefit','sax3'}
        result = segment('reshape2layout',squeeze(SET(no).MaR.Result(:,:,SET(no).CurrentTimeFrame,:)),no,panel);
    case {'one','ortho'}
        result = SET(no).MaR.Result(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
    otherwise
        myfailed('Unknown viewtype in mar-drawhelper.');
end;

if (sum(result(:))>0)
    
    %Ensure that we delete stuff
    set(DATA.Handles.imageaxes(panel),'NextPlot','add');
    
    %Normal MaR outline
    [c,DATA.Handles.marcontour{panel}] = contour(DATA.Handles.imageaxes(panel),double(result),[0.5 0.5]);
    set(DATA.Handles.marcontour{panel},'linecolor',[1 1 1]);
    if DATA.Pref.LineWidth>0
        set(DATA.Handles.marcontour{panel},'linewidth',DATA.Pref.LineWidth);
    else
        set(DATA.Handles.marcontour{panel},'visible','off');
    end;
else
    DATA.Handles.marcontour{panel} = [];
end;

%------------------------
function showmaredits(no)
%------------------------
%Show MaR edits on screen as a temporary overlay.
global DATA SET NO

if nargin==0
    no = NO;
end;

if isempty(SET(no).MaR)
    return;
end;

%Check if menu is enabled
if isequal(get(DATA.Handles.hideoverlayicon,'state'),'on')
    manualinteraction = false;
else
    manualinteraction = true;
end;

if not(manualinteraction)
    return;
end;

tempnos=no;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
    return;
end

panelstodo = find(DATA.ViewPanels==no);

for panel=panelstodo
    
    %Determine if it is a RGB image or not.
    if not(DATA.Silent)
        temp = get(DATA.Handles.imagehandle(panel),'CData');
        if ndims(temp)>2
            isrgbimage = true;
        else
            isrgbimage = false;
        end;
    else
        isrgbimage=false;
    end
    if isempty(SET(no).Colormap)
        isrgbimage = false;
    end
    
    clear temp;
    
    %Ok lets draw it
    switch DATA.ViewPanelsType{panel}
        case {'one','mmodespatial','ortho'}
            temp = SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
            
            if isrgbimage  %RGB iamge
                tempuint8 = min(uint8(255),uint8(1)+calcfunctions('remapuint8',temp,no,SET(no).Colormap));
            else  %Default grayscale image
                tempuint8 = min(uint8(255),uint8(1)+calcfunctions('remapuint8',temp));
                tempuint8 = repmat(tempuint8,[1 1 3]); %EH:
            end;
            clear temp;
            sz=size(tempuint8);
            marimrgb = reshape(tempuint8,[prod(sz(1:2)) 3]);
            
            tmp=SET(no).MaR.Manual(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
            tmp=tmp(:);
            marimrgb((tmp==int8( 1)),2) = uint8(255);
            marimrgb((tmp==int8(-1)),3) = uint8(255);
            marimrgb=reshape(marimrgb,[sz(1:2) 3]); %EH: added 3
            
        case {'montage','montagerow','montagefit','sax3'}
            % Convert to 2D and layout
            
            if isrgbimage  %RGB iamge
                tempuint8 = calcfunctions('remapuint8',...
                    segment('reshape2layout',squeeze(SET(no).IM(:,:,SET(no).CurrentTimeFrame,:)),no,panel),...
                    no,SET(no).Colormap);
            else  %Default grayscale image
                tempuint8 = calcfunctions('remapuint8',...
                    segment('reshape2layout',squeeze(SET(no).IM(:,:,SET(no).CurrentTimeFrame,:)),no,panel),...
                    no,calcfunctions('returnmapping',no,true));
            end
            tempuint8 = min(uint8(255),uint8(1)+tempuint8);
            sz = size(tempuint8);
            marimrgb = reshape(tempuint8,[prod(sz(1:2)) 3]);
            
            tmp=segment('reshape2layout',SET(no).MaR.Manual(:,:,SET(no).CurrentTimeFrame,:),no,panel);
            tmp=tmp(:);
            marimrgb((tmp==int8( 1)),2) = uint8(255);
            marimrgb((tmp==int8(-1)),3) = uint8(255);
            marimrgb=reshape(marimrgb,sz);
    end
    if not(DATA.Silent)
        set(DATA.Handles.imagehandle(panel),'CData',marimrgb);
    end
end

%---------------------
function showedits(no)
%---------------------
%Show edits (overlays) of scar and mar if current mode
global DATA SET

if isequal(DATA.CurrentTheme,'scar') || ...
        ~isempty(SET(no).Scar) && ~isempty(SET(no).Scar.GreyZone.map) && ...
        max(SET(no).Scar.GreyZone.map(:)) > 0
    if isequal(get(DATA.Handles.hidescaricon,'state'),'off')
        showviabilityedits(no);
    end
elseif isequal(DATA.CurrentTheme,'mar')
    if isequal(get(DATA.Handles.hidemaricon,'state'),'off')
        showmaredits(no);
    end
end

%-------------------------------------
function viewupdatetextposition(panel)
%-------------------------------------
% Updates the location of the corner text, to always stay still.
% From drawx(), it is called panelwise, since the linked panels aren't
% ready when the first panel is drawn. From other situations, it is called
% without arguements, and handles the linkage by itself. /JU

global DATA SET NO

%Input to force black background around corner text to improve visibility
if DATA.Pref.BackgroundColor
    bgcolor = 'blue';
else
    bgcolor = 'none';
end
myset(DATA.Handles.imagetypetext,'BackgroundColor',bgcolor); %Use myset to handle NaNs in vector
myset(DATA.Handles.seriesdescriptiontext,'BackgroundColor',bgcolor)
myset(DATA.Handles.dicomimagetypetext,'BackgroundColor',bgcolor)
myset(DATA.Handles.slicetimetext,'BackgroundColor',bgcolor)

if nargin==1
    ind = panel;
else
    no = NO;
    nos=SET(no).Linked;
    ind = [];
    for loop=1:length(nos)
        ind = [ind find(DATA.ViewPanels==nos(loop))]; %#ok<AGROW>
    end;
    ind = unique(ind); %remove duplicates.
    ind = ind(~strcmp(DATA.ViewPanelsType(ind),'mmodetemporal'));
end

for panel=ind
    xlim=get(DATA.Handles.imageaxes(panel),'xlim');
    ylim=get(DATA.Handles.imageaxes(panel),'ylim');
    
    extent = get(DATA.Handles.imagetypetext(panel),'extent');
    set(DATA.Handles.imagetypetext(panel),'color','white','position',...
        [xlim(2)-extent(3) ...
        ylim(2)-extent(4) 0]);
    
    extent = get(DATA.Handles.dicomimagetypetext(panel),'extent');
    set(DATA.Handles.dicomimagetypetext(panel),'color','white','position',...
        [xlim(2)-extent(3) ...
        ylim(2)-2*extent(4) 0]);
    
    extent = get(DATA.Handles.seriesdescriptiontext(panel),'extent');
    set(DATA.Handles.seriesdescriptiontext(panel),'color','white','position',...
        [xlim(2)-extent(3) ...
        ylim(2)-3*extent(4) 0]);
    
    extent = get(DATA.Handles.slicetimetext(panel),'extent');
    set(DATA.Handles.slicetimetext(panel),'color','white','position',...
        [xlim(2)-extent(3) ...
        ylim(2)-4*extent(4) 0]);
end

%------------------------------
function viewupdateannotext(panel)
%------------------------------
% Updates the visibility of point/measurement text.
% Respects hideX settings.
%
% If point is inbound, pointtext is visible.
% If measurementtext is inbound, measurementtext is visible.
% If whole ROI is inbound, place XXXX
% If only some of ROI is inbound, place at any inbound point
% If whole ROI is outside, make ROItext not visible.

global DATA SET NO

if nargin==1
    ind = panel;
    no = DATA.ViewPanels(panel(1));
else
    no = NO;
    nos=SET(no).Linked;
    ind = [];
    for loop=1:length(nos)
        ind = [ind find(DATA.ViewPanels==nos(loop))]; %#ok<AGROW>
    end;
    ind = unique(ind); %remove duplicates.
end

nom = no;
if ~isempty(SET(no).Parent)
    nom=SET(no).Parent;
end

for panel=ind
    if ~strcmp(DATA.ViewPanelsType(panel),'mmodetemporal')
        xlim=get(DATA.Handles.imageaxes(panel),'xlim');
        ylim=get(DATA.Handles.imageaxes(panel),'ylim');
        
        if ~isempty(SET(nom).Point)
            for loop=1:length(DATA.Handles.pointo{panel})
                px=get(DATA.Handles.pointo{panel}(loop),'xdata');
                py=get(DATA.Handles.pointo{panel}(loop),'ydata');
                if ~isnan(px)
                    if (px>xlim(1))&&(py>ylim(1))&&(px<xlim(2))&&(py<ylim(2))
                        set(DATA.Handles.pointtext{panel}(loop),'Position',[px+2 py]);
                    else
                        set(DATA.Handles.pointtext{panel}(loop),'Position',[nan nan]);
                    end
                end
            end
        end
        if ~isempty(SET(nom).Measure)
            for loop=1:length(DATA.Handles.measureline{panel})
                mx=get(DATA.Handles.measureline{panel}{loop}(end),'xdata');
                my=get(DATA.Handles.measureline{panel}{loop}(end),'ydata');
                if ~isnan(mx(1))
                    [mx,ix] = max(mx);
                    my = my(ix);
                    if (mx>xlim(1))&&(my>ylim(1))&&(mx<xlim(2))&&(my<ylim(2))
                        set(DATA.Handles.measuretext{panel}(loop),'Position',[mx+1 my]);
                    else
                        set(DATA.Handles.measuretext{panel}(loop),'Position',[nan nan]);
                    end
                end
            end
        end
        
        if ~ismember(DATA.ViewPanelsType{panel},{'montage','montagerow','montagefit','sax3'})
            for loop=1:SET(nom).RoiN
                rx=get(DATA.Handles.roicontour{panel}(loop),'xdata');
                ry=get(DATA.Handles.roicontour{panel}(loop),'ydata');
                if ~isnan(rx)
                    rbool=(rx>xlim(1))&(ry>ylim(1))&(rx<xlim(2))&(ry<ylim(2));
                    if all(rbool)
                        [ymin,ix] = min(SET(nom).Roi(loop).Y(:,SET(nom).CurrentTimeFrame));
                        set(DATA.Handles.roitext{panel}(loop),'position',...
                            [ymin-1 SET(nom).Roi(loop).X(ix,SET(no).CurrentTimeFrame)]);
                    elseif any(rbool)
                        inpoint=find(rbool);
                        %Take inbound point somewhere in the middle, put text there
                        inpoint=inpoint(round(length(inpoint)/2));
                        %inpoint = min(rx(rbool)) == rx;
                        set(DATA.Handles.roitext{panel}(loop),'position',...
                            [SET(nom).Roi(loop).Y(inpoint,SET(nom).CurrentTimeFrame) ...
                            SET(nom).Roi(loop).X(inpoint,SET(nom).CurrentTimeFrame)]);
                    else
                        set(DATA.Handles.roitext{panel}(loop),'Position',[nan nan]);
                    end
                end
            end
        end
    end
end

%------------------------
function updatevisibility
%------------------------
%Make sure visibility of handles correspond to status of hide/show icons
global DATA

%Make everything visible first

%Center cross
set([DATA.Handles.center{:}],'Visible','on');

%Pins
pinh = [DATA.Handles.endopin ...
    DATA.Handles.epipin ...
    DATA.Handles.rvendopin ...
    DATA.Handles.rvepipin];
myset(pinh,'Visible','on');

%Plane intersections
set([DATA.Handles.planeintersectionline{:}],'Visible','on');

%Contour intersections
ish = mycat(2,DATA.Handles.endointersectionline{:}, ...
    DATA.Handles.endointersectionpoints{:}, ...
    DATA.Handles.epiintersection{:}, ...
    DATA.Handles.rvendointersection{:}, ...
    DATA.Handles.rvepiintersection{:});
myset(ish,'Visible','on');

%Interpolation points
iph = mycat(2,DATA.Handles.endointerp, ...
    DATA.Handles.epiinterp, ...
    DATA.Handles.rvendointerp, ...
    DATA.Handles.rvepiinterp);
  
myset(iph,'Visible','on');

%LV segmentation
lvh = mycat(2, DATA.Handles.endocontour, ...
  DATA.Handles.epicontour,...
  DATA.Handles.endointerp, ...
  DATA.Handles.epiinterp, ...
  DATA.Handles.endopin, ...
  DATA.Handles.epipin, ...
  DATA.Handles.endointersectionline{:}, ...
  DATA.Handles.endointersectionpoints{:}, ...
  DATA.Handles.epiintersection{:});

myset(lvh,'Visible','on');

%RV segmentation
rvh = mycat(2,DATA.Handles.rvendocontour, ...
    DATA.Handles.rvepicontour, ...
    DATA.Handles.rvendointerp, ...
    DATA.Handles.rvepiinterp, ...
    DATA.Handles.rvendopin, ...
    DATA.Handles.rvepipin, ...
    DATA.Handles.rvendointersection{:}, ...
    DATA.Handles.rvepiintersection{:});

myset(rvh,'Visible','on');

%ROI
roih = [DATA.Handles.roicontour{:} ...
    DATA.Handles.roitext{:}];
set(roih,'Visible','on');

%Measurements
mlh = [DATA.Handles.measureline{:}];
mh = [mlh{:} ...
    DATA.Handles.measuretext{:}];
set(mh,'Visible','on');

%Annotation points
ph = [DATA.Handles.pointo{:} ...
    DATA.Handles.pointp{:} ...
    DATA.Handles.pointtext{:}];
set(ph,'Visible','on');

%Text
txth1 = [DATA.Handles.imagetypetext ...
    DATA.Handles.seriesdescriptiontext ...
    DATA.Handles.dicomimagetypetext ...
    DATA.Handles.slicetimetext];
txth2 = [DATA.Handles.roitext{:} ...
    DATA.Handles.measuretext{:} ...
    DATA.Handles.pointtext{:}];
myset([txth1 txth2],'Visible','on');

myset([DATA.Handles.marcontour{:}],'Visible','on');
myset([DATA.Handles.scarcontour{:}],'Visible','on');   %Scar contour
myset([DATA.Handles.weightedscarcontour{:}],'Visible','on'); %Weighted scar contour
myset([DATA.Handles.mocontour{:}],'Visible','on'); %microvascular obstruction contour
myset([DATA.Handles.moextentcontour{:}],'Visible','on'); %microvascular obstruction contour

%Then hide what is to hide

%Center cross
stateandicon=segment('iconson',{'hideplus','hideintersections','hideothercontour',...
  'hideinterp','hidelv','hiderv','hideroi','hidemeasure','hidepoint',...
  'hidetext','hidemar','hidescar'});
state=stateandicon{1,1};
if state%strcmp(get(DATA.Handles.hideplusicon,'state'),'on')
    set([DATA.Handles.center{:}],'Visible','off');
end

% stateandicon=segment('iconson','hidepins');
% state=stateandicon{1};
% %Pins
% if state%strcmp(get(DATA.Handles.hidepinsicon,'state'),'on')
%     myset(pinh,'Visible','off');
% end

%stateandicon=segment('iconson','hideintersections');
state=stateandicon{2,1};
%Plane intersections
if state%strcmp(get(DATA.Handles.hideintersectionsicon,'state'),'on')
    set([DATA.Handles.planeintersectionline{:}],'Visible','off');
end

%stateandicon=segment('iconson','hideothercontour');
state=stateandicon{3,1};
%Contour intersections
if state%strcmp(get(DATA.Handles.hideothercontouricon,'state'),'on')
    myset(ish,'Visible','off');
end

%stateandicon=segment('iconson','hideinterp');
state=stateandicon{4,1};
%Interpolation points
if state%strcmp(get(DATA.Handles.hideinterpicon,'state'),'on')
    myset(iph,'Visible','off');
end

%stateandicon=segment('iconson','hidelv');
state=stateandicon{5,1};
%LV segmentation
if state%strcmp(get(DATA.Handles.hidelvicon,'state'),'on')
    myset(lvh,'Visible','off');
end

%stateandicon=segment('iconson','hiderv');
state=stateandicon{6,1};
%RV segmentation
if state%strcmp(get(DATA.Handles.hidervicon,'state'),'on')
    myset(rvh,'Visible','off');
end

%stateandicon=segment('iconson','hideroi');
state=stateandicon{7,1};
%ROI
if state%strcmp(get(DATA.Handles.hideroiicon,'state'),'on')
    set(roih,'Visible','off');
end

%stateandicon=segment('iconson','hidemeasure');
state=stateandicon{8,1};
%Measurements
if state%strcmp(get(DATA.Handles.hidemeasuresicon,'state'),'on')
    set(mh,'Visible','off');
end

%stateandicon=segment('iconson','hidepoint');
state=stateandicon{9,1};
%Annotation points
if state%strcmp(get(DATA.Handles.hidepointsicon,'state'),'on')
    set(ph,'Visible','off');
end

%stateandicon=segment('iconson','hidetext');
state=stateandicon{10,1};
%Text
if state%strcmp(get(DATA.Handles.hidetexticon,'state'),'on')
    myset([txth1 txth2],'Visible','off');
end

%stateandicon=segment('iconson','hidemar');
state=stateandicon{11,1};

if state%strcmp(get(DATA.Handles.hidetexticon,'state'),'on')
    myset([DATA.Handles.marcontour{:}],'Visible','off');
    %myset(DATA.Handles.moextentcontour,'Visible','off'); %microvascular obstruction contour
end

%stateandicon=segment('iconson','hidescar');
state=stateandicon{12,1};

if state%strcmp(get(DATA.Handles.hidetexticon,'state'),'on')
  myset([DATA.Handles.scarcontour{:}],'Visible','off');   %Scar contour
  myset([DATA.Handles.mocontour{:}],'Visible','off');
  myset([DATA.Handles.moextentcontour{:}],'Visible','off'); %microvascular obstruction contour
  myset([DATA.Handles.weightedscarcontour{:}],'Visible','off'); %Weighted scar contour  
end
%------------------------------------
function updatemodeldisplay(no,panel)
%------------------------------------
%Update data used to correctly display segmentation in montage view. This
%fcn needs to be called when the segmentation has changed.
global DATA SET NO

if DATA.Run
    %No need to update if movie is playing
    return
end

if nargin < 2
    panel = DATA.CurrentPanel;
end
if nargin==0
    no = NO;
end;

%--- Check for size of structures
%EndoX
if ~isempty(SET(no).EndoX)
    SET(no).EndoXView = nan((length(SET(no).EndoX)+1)*SET(no).ZSize,SET(no).TSize);
    SET(no).EndoYView = SET(no).EndoXView;
else
    SET(no).EndoXView = NaN;
    SET(no).EndoYView = NaN;
end;
%RVEndo
if ~isempty(SET(no).RVEndoX)
    SET(no).RVEndoXView = nan((length(SET(no).RVEndoX)+1)*SET(no).ZSize,SET(no).TSize);
    SET(no).RVEndoYView = SET(no).RVEndoXView;
else
    SET(no).RVEndoXView = NaN;
    SET(no).RVEndoYView = NaN;
end;
%RVEpi
if ~isempty(SET(no).RVEpiX)
    SET(no).RVEpiXView = nan((length(SET(no).RVEpiX)+1)*SET(no).ZSize,SET(no).TSize);
    SET(no).RVEpiYView = SET(no).RVEpiXView;
else
    SET(no).RVEpiXView = NaN;
    SET(no).RVEpiYView = NaN;
end;
%Epi
if ~isempty(SET(no).EpiX)
    SET(no).EpiXView = nan((length(SET(no).EpiX)+1)*SET(no).ZSize,SET(no).TSize);
    SET(no).EpiYView = SET(no).EpiXView;
else
    SET(no).EpiXView = NaN;
    SET(no).EpiYView = NaN;
end;

temp = cell(SET(no).TSize,1);

if isempty(SET(no).EndoPinX)
    SET(no).EndoPinXView = [];
    SET(no).EndoPinYView = [];
else
    SET(no).EndoPinXView = temp;
    SET(no).EndoPinYView = temp;
end;
if isempty(SET(no).EpiPinX)
    SET(no).EpiPinXView = [];
    SET(no).EpiPinYView = [];
else
    SET(no).EpiPinXView = temp;
    SET(no).EpiPinYView = temp;
end;
if isempty(SET(no).RVEndoPinX)
    SET(no).RVEndoPinXView = [];
    SET(no).RVEndoPinYView = [];
else
    SET(no).RVEndoPinXView = temp;
    SET(no).RVEndoPinYView = temp;
end;
if isempty(SET(no).RVEpiPinX)
    SET(no).RVEpiPinXView = [];
    SET(no).RVEpiPinYView = [];
else
    SET(no).RVEpiPinXView = temp;
    SET(no).RVEpiPinYView = temp;
end;

if isempty(SET(no).EndoInterpX)
    SET(no).EndoInterpXView = [];
    SET(no).EndoInterpYView = [];
else
    SET(no).EndoInterpXView = temp;
    SET(no).EndoInterpYView = temp;
end;
if isempty(SET(no).EpiInterpX)
    SET(no).EpiInterpXView = [];
    SET(no).EpiInterpYView = [];
else
    SET(no).EpiInterpXView = temp;
    SET(no).EpiInterpYView = temp;
end;
if isempty(SET(no).RVEndoInterpX)
    SET(no).RVEndoInterpXView = [];
    SET(no).RVEndoInterpYView = [];
else
    SET(no).RVEndoInterpXView = temp;
    SET(no).RVEndoInterpYView = temp;
end;
if isempty(SET(no).RVEpiInterpX)
    SET(no).RVEpiInterpXView = [];
    SET(no).RVEpiInterpYView = [];
else
    SET(no).RVEpiInterpXView = temp;
    SET(no).RVEpiInterpYView = temp;
end;

%TODO: Semi-vectorize the contours so they can fit into the optimized
% single-loop below this double-loop. /JU

for zloop=1:SET(no).ZSize
    [xofsall,yofsall] = calcfunctions('calcoffset',zloop,[],no,panel); %[] was 'montage'
    xofs = xofsall;
    yofs = yofsall;
    for tloop=1:SET(no).TSize
        if numel(xofsall) > 1
            xofs = xofsall(tloop);
            yofs = yofsall(tloop);
        end
        %zloop = SET(no).SAX3.slices{loop,tloop};
        
        %Endocontour
        if ~isempty(SET(no).EndoX)
            SET(no).EndoXView((1+(zloop-1)*(length(SET(no).EndoX)+1)):(zloop*(length(SET(no).EndoX)+1)-1),tloop) = ...
                SET(no).EndoX(:,tloop,zloop)+xofs;
            SET(no).EndoYView((1+(zloop-1)*(length(SET(no).EndoX)+1)):(zloop*(length(SET(no).EndoX)+1)-1),tloop) = ...
                SET(no).EndoY(:,tloop,zloop)+yofs;
        end;
        %Epicontour
        if ~isempty(SET(no).EpiX)
            SET(no).EpiXView((1+(zloop-1)*(length(SET(no).EpiX)+1)):(zloop*(length(SET(no).EpiX)+1)-1),tloop) = ...
                SET(no).EpiX(:,tloop,zloop)+xofs;
            SET(no).EpiYView((1+(zloop-1)*(length(SET(no).EpiX)+1)):(zloop*(length(SET(no).EpiX)+1)-1),tloop) = ...
                SET(no).EpiY(:,tloop,zloop)+yofs;
        end;
        %RV Endocontour
        if ~isempty(SET(no).RVEndoX)
            SET(no).RVEndoXView((1+(zloop-1)*(length(SET(no).RVEndoX)+1)):(zloop*(length(SET(no).RVEndoX)+1)-1),tloop) = ...
                SET(no).RVEndoX(:,tloop,zloop)+xofs;
            SET(no).RVEndoYView((1+(zloop-1)*(length(SET(no).RVEndoX)+1)):(zloop*(length(SET(no).RVEndoX)+1)-1),tloop) = ...
                SET(no).RVEndoY(:,tloop,zloop)+yofs;
        end;
        %RV Epicontour
        if ~isempty(SET(no).RVEpiX)
            SET(no).RVEpiXView((1+(zloop-1)*(length(SET(no).RVEpiX)+1)):(zloop*(length(SET(no).RVEpiX)+1)-1),tloop) = ...
                SET(no).RVEpiX(:,tloop,zloop)+xofs;
            SET(no).RVEpiYView((1+(zloop-1)*(length(SET(no).RVEpiX)+1)):(zloop*(length(SET(no).RVEpiX)+1)-1),tloop) = ...
                SET(no).RVEpiY(:,tloop,zloop)+yofs;
        end;
        
        %     %Endopin
        %     if ~isempty(SET(no).EndoPinX)
        %       SET(no).EndoPinXView{tloop} = [...
        %         SET(no).EndoPinXView{tloop};...
        %         SET(no).EndoPinX{tloop,zloop}+xofs];
        %       SET(no).EndoPinYView{tloop} = [...
        %         SET(no).EndoPinYView{tloop};...
        %         SET(no).EndoPinY{tloop,zloop}+yofs];
        %     end;
        %     %Epipin
        %     if ~isempty(SET(no).EpiPinX)
        %       SET(no).EpiPinXView{tloop} = [...
        %         SET(no).EpiPinXView{tloop};...
        %         SET(no).EpiPinX{tloop,zloop}+xofs];
        %       SET(no).EpiPinYView{tloop} = [...
        %         SET(no).EpiPinYView{tloop};...
        %         SET(no).EpiPinY{tloop,zloop}+yofs];
        %     end;
        %     %RVEndopin
        %     if ~isempty(SET(no).RVEndoPinX)
        %       SET(no).RVEndoPinXView{tloop} = [...
        %         SET(no).RVEndoPinXView{tloop};...
        %         SET(no).RVEndoPinX{tloop,zloop}+xofs];
        %       SET(no).RVEndoPinYView{tloop} = [...
        %         SET(no).RVEndoPinYView{tloop};...
        %         SET(no).RVEndoPinY{tloop,zloop}+yofs];
        %     end;
        %     %RVEpipin
        %     if ~isempty(SET(no).RVEpiPinX)
        %       SET(no).RVEpiPinXView{tloop} = [...
        %         SET(no).RVEpiPinXView{tloop};...
        %         SET(no).RVEpiPinX{tloop,zloop}+xofs];
        %       SET(no).RVEpiPinYView{tloop} = [...
        %         SET(no).RVEpiPinYView{tloop};...
        %         SET(no).RVEpiPinY{tloop,zloop}+yofs];
        %     end;
        %
        %     %Endo Interp Pts
        %     if ~isempty(SET(no).EndoInterpX)
        %       SET(no).EndoInterpXView{tloop} = [...
        %         SET(no).EndoInterpXView{tloop};...
        %         SET(no).EndoInterpX{tloop,zloop}+xofs];
        %       SET(no).EndoInterpYView{tloop} = [...
        %         SET(no).EndoInterpYView{tloop};...
        %         SET(no).EndoInterpY{tloop,zloop}+yofs];
        %     end;
        %     %Epi Interp Pts
        %     if ~isempty(SET(no).EpiInterpX)
        %       SET(no).EpiInterpXView{tloop} = [...
        %         SET(no).EpiInterpXView{tloop};...
        %         SET(no).EpiInterpX{tloop,zloop}+xofs];
        %       SET(no).EpiInterpYView{tloop} = [...
        %         SET(no).EpiInterpYView{tloop};...
        %         SET(no).EpiInterpY{tloop,zloop}+yofs];
        %     end;
        %     %RVEndo Interp Pts
        %     if ~isempty(SET(no).RVEndoInterpX)
        %       SET(no).RVEndoInterpXView{tloop} = [...
        %         SET(no).RVEndoInterpXView{tloop};...
        %         SET(no).RVEndoInterpX{tloop,zloop}+xofs];
        %       SET(no).RVEndoInterpYView{tloop} = [...
        %         SET(no).RVEndoInterpYView{tloop};...
        %         SET(no).RVEndoInterpY{tloop,zloop}+yofs];
        %     end;
        %     %RVEpi Interp Pts
        %     if ~isempty(SET(no).RVEpiInterpX)
        %       SET(no).RVEpiInterpXView{tloop} = [...
        %         SET(no).RVEpiInterpXView{tloop};...
        %         SET(no).RVEpiInterpX{tloop,zloop}+xofs];
        %       SET(no).RVEpiInterpYView{tloop} = [...
        %         SET(no).RVEpiInterpYView{tloop};...
        %         SET(no).RVEpiInterpY{tloop,zloop}+yofs];
        %     end;
        
    end; %tloop
end; %zloop

%check to eliminate 1*0 matrix in Interpcells
if ~isempty(SET(no).EndoInterpX)
    [SET(no).EndoInterpX{cellfun('isempty',SET(no).EndoInterpX)}]=deal([]);
    [SET(no).EndoInterpY{cellfun('isempty',SET(no).EndoInterpY)}]=deal([]);
end
if ~isempty(SET(no).EpiInterpX)
    [SET(no).EpiInterpX{cellfun('isempty',SET(no).EpiInterpX)}]=deal([]);
    [SET(no).EpiInterpY{cellfun('isempty',SET(no).EpiInterpY)}]=deal([]);
end
if ~isempty(SET(no).RVEndoInterpX)
    [SET(no).RVEndoInterpX{cellfun('isempty',SET(no).RVEndoInterpX)}]=deal([]);
    [SET(no).RVEndoInterpY{cellfun('isempty',SET(no).RVEndoInterpY)}]=deal([]);
end
if ~isempty(SET(no).RVEpiInterpX)
    [SET(no).RVEpiInterpX{cellfun('isempty',SET(no).RVEpiInterpX)}]=deal([]);
    [SET(no).RVEpiInterpY{cellfun('isempty',SET(no).RVEpiInterpY)}]=deal([]);
end

[cellx,celly]=calcfunctions('calcoffsetcells',no,panel);
for tloop=1:SET(no).TSize
    %Interp Pts
    if ~isempty(SET(no).EndoInterpX) && tloop<=size(SET(no).EndoInterpX,1)%&&~isempty(cell2mat(SET(no).EndoInterpX(tloop,:)))
        SET(no).EndoInterpXView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).EndoInterpX(tloop,:)',cellx,'UniformOutput',0));
        SET(no).EndoInterpYView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).EndoInterpY(tloop,:)',celly,'UniformOutput',0));
    end
    if ~isempty(SET(no).EpiInterpX) && tloop<=size(SET(no).EpiInterpX,1)
        SET(no).EpiInterpXView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).EpiInterpX(tloop,:)',cellx,'UniformOutput',0));
        SET(no).EpiInterpYView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).EpiInterpY(tloop,:)',celly,'UniformOutput',0));
    end
    if ~isempty(SET(no).RVEndoInterpX)&& tloop<=size(SET(no).RVEndoInterpX,1)
        SET(no).RVEndoInterpXView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).RVEndoInterpX(tloop,:)',cellx,'UniformOutput',0));
        SET(no).RVEndoInterpYView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).RVEndoInterpY(tloop,:)',celly,'UniformOutput',0));
    end
    if ~isempty(SET(no).RVEpiInterpX) && tloop<=size(SET(no).RVEpiInterpX,1)
        SET(no).RVEpiInterpXView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).RVEpiInterpX(tloop,:)',cellx,'UniformOutput',0));
        SET(no).RVEpiInterpYView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).RVEpiInterpY(tloop,:)',celly,'UniformOutput',0));
    end
    
    %Pins
    if ~isempty(SET(no).EndoPinX)
        SET(no).EndoPinXView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).EndoPinX(tloop,:)',cellx,'UniformOutput',0));
        SET(no).EndoPinYView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).EndoPinY(tloop,:)',celly,'UniformOutput',0));
    end
    if ~isempty(SET(no).EpiPinX)
        SET(no).EpiPinXView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).EpiPinX(tloop,:)',cellx,'UniformOutput',0));
        SET(no).EpiPinYView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).EpiPinY(tloop,:)',celly,'UniformOutput',0));
    end
    if ~isempty(SET(no).RVEndoPinX)
        SET(no).RVEndoPinXView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).RVEndoPinX(tloop,:)',cellx,'UniformOutput',0));
        SET(no).RVEndoPinYView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).RVEndoPinY(tloop,:)',celly,'UniformOutput',0));
    end
    if ~isempty(SET(no).RVEpiPinX)
        SET(no).RVEpiPinXView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).RVEpiPinX(tloop,:)',cellx,'UniformOutput',0));
        SET(no).RVEpiPinYView{tloop}=cell2mat(cellfun(@plus,...
            SET(no).RVEpiPinY(tloop,:)',celly,'UniformOutput',0));
    end
end

%-------------------------------
function s = cellref(a,varargin)
%-------------------------------
%Returns vector of content of e(varargin), for all elements e of cell a
if ~iscell(a{1})
    s = nan(1,numel(a)*sum(varargin{1}~=0));
else
    s = cell(1,numel(a)*sum(varargin{1}~=0));
end
n = numel(s)/numel(a);
for i = 1:numel(a)
    s((i-1)*n+(1:n)) = a{i}(varargin{:});
end

%----------------------------------
function updateglazoomstate(no,ysz)
%----------------------------------
global SET
xzoomsz = max(SET(no).XSize,SET(no).YSize);
repos = (ysz-xzoomsz)/2;
SET(no).GLA.ZoomState = 0.5+[repos+[0;xzoomsz];0;SET(no).ZSize];