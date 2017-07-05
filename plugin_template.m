function varargout = plugin_template(fcn,varargin)
%-------------------------------------------------
%Demo code to produce own plug-ins to segment. Use this code as a template
%to write own functions.

%Einar Heiberg

if nargin==0
  myfailed('Expects at least one input argument.');
  return;
end;

switch fcn
  case 'getname'
    varargout = cell(1,1);
    
		%Segment with versions >1.636 calls with two input arguments where
		%the second input argument is the handle to the menu item.
    
		%Register submenus. You need to change label and callback here!!!
		uimenu(varargin{1},'Label','Test A','Callback','plugin_template(''funa_Callback'')');
		uimenu(varargin{1},'Label','Test B','Callback','plugin_template(''funb_Callback'')');
    
    %Here write your own title that will be shown on menu.
    varargout{1} = 'Example plug-in template'; 
    
    %The above code is an example of a plug-in with subfunctions. To create
    %a simple function comment the above code and uncomment the code below.
    
    %set(varargin{1},'Label','Demo plugin','Callback','plugin_template(''funa'')');
    %varargout{1} = 'Plug-in template';     
    set(varargin{1},'Callback',''); 
  case 'getdependencies'
    %Here: List all depending files. This is required if your plugin should
    %be possible to compile to the stand-alone version of Segment.
    varargout = cell(1,4);
    
    %M-files, list as {'hello.m' ...};
    varargout{1} = {};

    %Fig-files, list as {'hello.fig' ... };
    varargout{2} = {};
    
    %Mat-files, list as {'hello.mat' ... };
    varargout{3} = {};
    
    %Mex-files, list as {'hello' ...}; %Note i.e no extension!!!
    varargout{4} = {};
    
  otherwise
    macro_helper(fcn,varargin{:}); %Future use to record macros
		[varargout{1:nargout}] = feval(fcn,varargin{:}); % FEVAL switchyard    
end;

%---------------------
function funa_Callback %#ok<DEFNU>
%---------------------
%Here your write all your code, you may use sub-functions.
%
%The global variables are:
%-DATA contains GUI information and edgedetected images
%-NO   current image stack, this is a scalar
%-SET  contains each image set info such as
%      IM (image data)
%      XSize (size in x etc)
%      ..
%      TIncr (time increment)
%      Resolution (in mm)
%      SliceThickness (in mm)
%      ...
%      EndoX (x-points of segmentation endocardie)

global DATA SET NO

%Usually good to check!
if not(DATA.DataLoaded)
  myfailed('No data loaded.');
  return;
end;
  
mymsgbox(dprintf('You have %d image stacks loaded.',length(SET)));

%---------------------
function funb_Callback %#ok<DEFNU>
%---------------------
myfailed('Not yet implemented.');