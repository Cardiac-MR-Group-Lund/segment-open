function varargout = atriumpenfunctions(varargin)
%This contains all callbacks related to atrium pen objects.

% Justine Le Douaron, Medviso, 2023-2024

%#ok<*GVMIS>

if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%---------------------------------------------------------
function createnewobject(atriumtype,no)
%---------------------------------------------------------
%Create a new Atrium object
arguments
  atriumtype {mustBeMember(atriumtype,{'la','ra'})}
  no = []
end
global DATA

if isempty(no)
  global NO %#ok<TLEV> we don't want to call the global variable unless no NO is passed as arguemnt
  no = NO;
end

%Instantiate a new object
DATA.GeneralPenSettings.createatriumobject(atriumtype,no);

%---------------------------------------------------------
function deleteobject_Callback(atriumtype)
%---------------------------------------------------------
%Delete an Atrium object. Callback for icons.
arguments
  atriumtype {mustBeMember(atriumtype,{'la','ra'})}
end

deleteobject(atriumtype);

%---------------------------------------------------------
function deleteobject(atriumtype,no)
%---------------------------------------------------------
%Delete an Atrium object

global DATA SET NO
% Argument check
if nargin < 1
  error('atriumtype is required');
end

% Validate atriumtype
if ~ismember(atriumtype, {'la','ra'})
  error('atriumtype must be ''la'' or ''ra''');
end
% Optional stack number
if nargin < 2 || isempty(no)
  no = NO;
end
tools('enableundo');

%Remove object from SET
SET(no).(upper(atriumtype)) = [];

%Graphical update
handles = DATA.Handles;
panel = DATA.CurrentPanel;
set([...
  handles.([atriumtype,'contour'])(panel), ...
  handles.([atriumtype,'interp'])(panel)],...
  'XData',nan,'YData',nan);

%Graphics and results panel update
update(atriumtype,no);

%---------------------------------------------------------
function clearslice_Callback(atriumtype)
%---------------------------------------------------------
%Delete selected slices in an Atrium object. Callback for icons.
arguments
  atriumtype {mustBeMember(atriumtype,{'la','ra'})}
end

clearslice(atriumtype);

%---------------------------------------------------------
function clearslice(atriumtype)
%---------------------------------------------------------
%Delete selected slices in an Atrium object
arguments
  atriumtype {mustBeMember(atriumtype,{'la','ra'})}
end
global DATA SET NO
no = NO;

tools('enableundo');
tf = SET(no).CurrentTimeFrame;

%Update SET
slices = SET(no).StartSlice:SET(no).EndSlice;
if ~isempty(SET(no).(upper(atriumtype))) && ~isempty(SET(no).(upper(atriumtype)).X)
  SET(no).(upper(atriumtype)).X(:,tf,slices) = NaN;
  SET(no).(upper(atriumtype)).Y(:,tf,slices) = NaN;
end
if ~isempty(SET(no).(upper(atriumtype))) && ~isempty(SET(no).(upper(atriumtype)).InterpX)
  SET(no).(upper(atriumtype)).InterpX(tf,slices) = {zeros(0,0)};
  SET(no).(upper(atriumtype)).InterpY(tf,slices) = {zeros(0,0)};
end

%Graphical update
handles = DATA.Handles;
panel = DATA.CurrentPanel;
set([...
  handles.([atriumtype,'contour'])(panel), ...
  handles.([atriumtype,'interp'])(panel)],...
  'XData',nan,'YData',nan);

%Graphics and results panel update
update(atriumtype,no);

%---------------------------------------------------------
function deleteallobjects(no)
%---------------------------------------------------------
%Delete all Atrium objects. Used in "clear all" function.
global NO
if nargin ==  0
  no = NO;
end

atriumtypes = {'la','ra'};
for loop = 1:numel(atriumtypes)
  deleteobject(atriumtypes{loop},no);
end

%---------------------------------------------------------
function deletethispoint_Callback
%---------------------------------------------------------
%Delete interpolation point in Atrium objects.
global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

if isempty(DATA.LastObject)
  return
end

pointind = DATA.LastObject(1);
slice = DATA.LastObject(2);
tf = DATA.LastObject(3);
pointtype = DATA.LastObjectType;
pointtype = pointtype(1:end-6);

SET(no).(pointtype).InterpX{tf,slice}(pointind) = [];
SET(no).(pointtype).InterpY{tf,slice}(pointind) = [];

%do resampling of curve when deleting single point.
X = SET(no).(pointtype).InterpX{tf,slice};
Y = SET(no).(pointtype).InterpY{tf,slice};
%removes duplicate points and resamples the contour
opencontour = true;
datanumpoints = SET(no).(pointtype).getnumpoints;
[x,y] = calcfunctions('resamplecurve',X,Y,datanumpoints,opencontour);
SET(no).(pointtype).Y(:,tf,slice) = y;
SET(no).(pointtype).X(:,tf,slice) = x;

%redo strain only if in ED
redo = true;
if ~isequal(tf,SET(no).EDT)
  redo = false;
end

%Update flag for manual contour 
if isequal(tf,SET(no).EST)
  SET(no).(pointtype).isManualContourinES = true;
end

%Graphics and results panel update
update(lower(pointtype),no,tf,redo);

%---------------------------------------------------------
function update(atriumtype,no,tf,redostrain)
%---------------------------------------------------------
%Helper function to perform all necessary updates after a change in contour
%in an Atrium object
arguments
  atriumtype {mustBeMember(atriumtype,{'la','ra'})}
  no = []
  tf = []
  redostrain = true
end
global SET

if isempty(no)
  global NO %#ok<TLEV> we don't want to call the global variable unless no NO is passed as arguemnt
  no = NO;
end

if isempty(tf)
  tf = SET(no).CurrentTimeFrame;
end

%Update draw list by calling drawno
drawfunctions('drawno',no);

%Update redo flag for strain
if redostrain && ~isempty(SET(no).StrainMitt)
  linkedlaxnos = findfunctions('findlaxnowithheartpart',atriumtype,SET(no).StrainMitt.LAXGroup);
  for loop = 1:numel(linkedlaxnos)
    if ~isempty(SET(linkedlaxnos(loop)).StrainMitt)
      updateredo(SET(linkedlaxnos(loop)).StrainMitt,upper(atriumtype));
    end
  end
end

%Update result panel
if (strcmp(atriumtype,'la')||strcmp(atriumtype,'ra'))  && ismember(tf,[SET(no).EST, SET(no).EDT])
  segment('updatemeasurement');
end

%--------------------------------
function duplicate(atriumtype,no)
%--------------------------------
% Helper function to duplicate/create a complete new instance of the object
global SET
atriumtype = upper(atriumtype);
if ~isempty(SET(no).(atriumtype))
  if isa(SET(no).(atriumtype),'generalpen.clatriumpen')
    newobj = copy(SET(no).(atriumtype));
    SET(no).(atriumtype) = newobj;
  end
end

%-------------------------------------
function answ = isemptycontour(atriumtype,no)
%-------------------------------------
% Helper function to check if LA/RA is empty or not
global SET
answ = true;
atriumtype = upper(atriumtype);
if ~isempty(SET(no).(atriumtype))
  if isa(SET(no).(atriumtype),'generalpen.clatriumpen')
    if ~isempty(SET(no).(atriumtype).X) &&  ~isempty(SET(no).(atriumtype).Y)
      answ = false;
    end
  end
end

%-------------------------------------
function [startx,endx,starty,endy] = getstartendxy(atriumtype,no)
%-------------------------------------
% Helper function to check if LA/RA is empty or not
global SET

atriumtype = upper(atriumtype);
startx = floor(min(SET(no).(atriumtype).X(:)));
endx = ceil(max(SET(no).(atriumtype).X(:)));
starty = floor(min(SET(no).(atriumtype).Y(:)));
endy = ceil(max(SET(no).(atriumtype).Y(:)));


