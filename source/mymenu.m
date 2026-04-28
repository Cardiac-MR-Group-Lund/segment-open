function [varargout] = mymenu(header,varargin)
%m = MYMENU(header,items,defaultstring)
%m = MYMENU(header,items,defaultstring,maxselections)
%m = MYMENU(header,item1,item2,item3,...)
%
%Same as menu, but with the following additions:
%- more elegant! uses keyboard to select items.
%- modal display
%- can use default string
%- last optional argument is figure handle
%
%See also MYINPUTSTRUCT, YESNO

%Einar Heiberg
global DATA %#ok<*GVMIS>

if nargout > 0
  varargout = cell(1,1);
end

if isa(DATA, 'maingui') && DATA.Testing
  testing = DATA.Testing;
else
  testing = false;
end

if testing
  %If testing then also DATA.Buffer should exist
  v = popfrombuffer('Mymenu');
  if isempty(v)
    myfailed('Menu buffer is empty.');
    return;
  elseif isnan(v)
    %Test of warning and nbr of options (one less than cell size)
    pushtobuffer('Warnings',sprintf('%d %s',numel(varargin{1})-1,header));
    v = 2;
  end
  varargout{1} = v;
  return;
end

% rewrite of mymenu to use myinputstruct
n = 1;
fieldstr = 'Input';
instruct(n).Field = fieldstr;

instruct(n).Label = dprintf('Select an option');

inputargs = varargin{:};
if ~iscell(inputargs)
  % do not unpack, if is no longer a cell array after the unpacking
  inputargs = varargin(:);
end
% check if all arguments are chars and translate only them
charind = cell2mat(cellfun(@(x) ischar(x),inputargs ,'UniformOutput', false));
defaultlabels = cellfun(@dprintf, inputargs(charind),'UniformOutput', false);
instruct(n).Default = defaultlabels;

[outs,ok] = myinputstruct(instruct,dprintf(header),10);
if ok
  varargout{1} = outs.(fieldstr);
else
  varargout{1} = 0;
end
