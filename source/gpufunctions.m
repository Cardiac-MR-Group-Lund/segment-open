function varargout = gpufunctions(varargin)
%Functions related to GPU settings

%#ok<*GVMIS>

%Justine Le Douaron, Medviso, 2024

if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%----------------------------------------------------------------
function execenv = getexecutionenvironement
%----------------------------------------------------------------
% get execution environement based on GPU compatibility
global DATA
persistent executionenv

if isempty(executionenv)
  if DATA.isGPUCompatible
    executionenv  = 'gpu';
  else
    executionenv  = 'cpu';
  end
end
execenv = executionenv;

%----------------------------------------------------------------
function minibatchsize = getminibatchsize(minibatchsize)
%----------------------------------------------------------------
% get "minibatchsize" DL functions based on user preferences
global DATA

if nargin < 1
  minibatchsize = 10; % default mini batch size
end

switch lower(DATA.Pref.GPUPerformance)
  case 'low'
    minibatchsize = floor(0.25*minibatchsize);
  case 'high'
    minibatchsize = 2*minibatchsize;
  case 'medium'
    return % take the default
  otherwise
    return % take the default
end