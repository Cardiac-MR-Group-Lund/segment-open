function varargout = plugin_CSPAMM(fcn,varargin)
%Reconstructs two CSPAMM images (vertical stripes and horisontal) by
%multiplication of the two images.
%
%Einar Heiberg, 2018-08-28

% Main plugin function. Acts as switchyard, also
% reports dependencys and register menu items.
global DATA

if nargin==0
  myfailed('Expects at least one input argument.');
  return;
end;

switch fcn
  case 'getname'
    varargout = cell(1,1);

    %Register callbacks
    uimenu(varargin{1},'Label','CSPAMM','Callback','plugin_CSPAMM(''cspamminit'')');
    set(varargin{1},'Callback','');

    %Return title
    varargout{1} = 'Reconstruct CSPAMM images';
    
  case 'getdependencies'
    %Here: List all depending files. This is required if your plugin should
    %be possible to compile to the stand-alone version of Segment.
    varargout = cell(1,4);

    %M-files
    varargout{1} = {};

    %Fig-files
    varargout{2} = {'mergestacks.fig'};

    %Mat-files
    varargout{3} = {};

    %Mex-files
    varargout{4} = {};

  otherwise
    macro_helper(fcn,varargin{:});
    [varargout{1:nargout}] = feval(fcn,varargin{:}); % FEVAL switchyard
end;

%-----------------------------------
function cspamminit
%-----------------------------------
global SET DATA

DATA.GUI.cspammreconstruct = mygui('mergestacks.fig');
handles = DATA.GUI.cspammreconstruct.handles;
stackc = cell(1,numel(SET));
for no = 1:numel(SET)
  stackc{no} = sprintf('%d. %s, %s',no,SET(no).ImageType,SET(no).ImageViewPlane);
end

set(DATA.GUI.cspammreconstruct.fig,'Name','Reconstruct CSPAMM');
set(handles.imagestackslistbox,'String',stackc);
set(handles.mergestackspushbutton,'String','Reconstruct');
set(handles.mergestackspushbutton,'Callback','plugin_CSPAMM(''cspammadjustment'')');


%--------------------------------
function cspammadjustment
%----------------------------------
global SET DATA

handles = DATA.GUI.cspammreconstruct.handles;

str = get(handles.imagestackslistbox,'String');
vals = get(handles.imagestackslistbox,'Value');

if length(vals)>2 || length(vals)<2
  mywarning('Need to select 2 image stacks for reconstruction')
  return
end

no1=str2num(str{vals(1)}(1));
no2=str2num(str{vals(2)}(1));
nonew = length(SET)+1;

SET(nonew) = SET(no1);
SET(nonew).Linked = nonew;
SET(nonew).IM = SET(no1).IM.*SET(no2).IM;
SET(nonew).IM = SET(nonew).IM/max(SET(nonew).IM(:));
SET(nonew).IntensityScaling = 1;
SET(nonew).IntensityOffset = 0;

%close the gui
close(DATA.GUI.cspammreconstruct.fig)
segment('viewrefreshall_Callback')%mymsgbox('Click refresh!');