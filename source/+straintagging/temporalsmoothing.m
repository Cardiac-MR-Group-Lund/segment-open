function varargout = temporalsmoothing(varargin)
%overwrites all existing temporal smoothing settings in current mat file
%according to 

macro_helper(varargin{:});
if nargin < 1 || isempty(varargin{1})
  varargin{1} = 'init';
end

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init
%------------
global DATA

% defaultpresets=0;
% 
% if ~isfield(DATA.GUI,'temporalsmoothing') 
%   defaultpresets=1;
% elseif isempty(DATA.GUI.temporalsmoothing) 
%   defaultpresets=1;
% elseif ~isfield(DATA.GUI.temporalsmoothing,'ts2ch')
%   defaultpresets=1;
% end

gui = mygui(fullfile('+straintagging','temporalsmoothing.fig'));
DATA.GUI.temporalsmoothing=gui;

%check that DATA.Pref.tsX is initialized
if ~isfield(DATA.Pref,'ts2ch')
  DATA.Pref.ts2ch=0;
end

if ~isfield(DATA.Pref,'ts3ch')
  DATA.Pref.ts3ch=0;
end

if ~isfield(DATA.Pref,'ts4ch')
  DATA.Pref.ts4ch=0;
end

if ~isfield(DATA.Pref,'tssax')
  DATA.Pref.tssax=0;
end

set(gui.handles.twochamberpopup,'value',3+DATA.Pref.ts2ch);
set(gui.handles.threechamberpopup,'value',3+DATA.Pref.ts3ch);
set(gui.handles.fourchamberpopup,'value',3+DATA.Pref.ts4ch);
set(gui.handles.saxpopup,'value',3+DATA.Pref.tssax);

%---------------------
function ok_callback
%---------------------

global SET DATA

%do warning as this overrides all current strain data.
answer=yesno('Resetting temporalsmoothing clears all strain data. Are you sure?');

if ~answer
  close(gcbf)
  return
end

gui = DATA.GUI.temporalsmoothing;

DATA.Pref.ts2ch=get(gui.handles.twochamberpopup,'value')-3;
DATA.Pref.ts3ch=get(gui.handles.threechamberpopup,'value')-3;
DATA.Pref.ts4ch=get(gui.handles.fourchamberpopup,'value')-3;
DATA.Pref.tssax=get(gui.handles.saxpopup,'value')-3;

for i=1:length(SET)
  redo=0;
  if strcmp(SET(i).ImageViewPlane,'2CH')
    if ~(isfield(SET(i).StrainTagging,'temporalsmoothing') && SET(i).StrainTagging.temporalsmoothing==DATA.Pref.ts2ch)
      if isfield(SET(i).StrainTagging,'transformparameters')
        redo=1;
      end
      SET(i).StrainTagging=[];
      SET(i).StrainTagging.temporalsmoothing=DATA.Pref.ts2ch;
      SET(i).StrainTagging.redo=redo;
    end
  elseif strcmp(SET(i).ImageViewPlane,'3CH')
    if ~(isfield(SET(i).StrainTagging,'temporalsmoothing') && SET(i).StrainTagging.temporalsmoothing==DATA.Pref.ts3ch)
      if isfield(SET(i).StrainTagging,'transformparameters')
        redo=1;
      end
      SET(i).StrainTagging=[];
      SET(i).StrainTagging.temporalsmoothing=DATA.Pref.ts3ch;
      SET(i).StrainTagging.redo=redo;
    end
  elseif strcmp(SET(i).ImageViewPlane,'4CH')
    if ~(isfield(SET(i).StrainTagging,'temporalsmoothing') && SET(i).StrainTagging.temporalsmoothing==DATA.Pref.ts4ch)
      if isfield(SET(i).StrainTagging,'transformparameters')
        redo=1;
      end
      SET(i).StrainTagging=[];
      SET(i).StrainTagging.temporalsmoothing=DATA.Pref.ts4ch;
      SET(i).StrainTagging.redo=redo;
    end
  else
    if ~(isfield(SET(i).StrainTagging,'temporalsmoothing') && SET(i).StrainTagging.temporalsmoothing==DATA.Pref.tssax)
      if isfield(SET(i).StrainTagging,'transformparameters')
        redo=1;
      end
      SET(i).StrainTagging=[];
      SET(i).StrainTagging.temporalsmoothing=DATA.Pref.tssax;
      SET(i).StrainTagging.redo=redo;
    end
  end
end

segpref('save_Callback');

close(gcbf)
  