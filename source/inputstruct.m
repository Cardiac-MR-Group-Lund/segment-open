function [res,varargout] = inputstruct(s,tit)
%INPUTSTRUCT creates a popup input dialog box from a struct
%
% INPUTSTRUCT(struct_var,title)
%
% Input is a struct with exploraty names, and a title for the dialogbox
% Output is the modified struct, and optional an ok variable.
%
% Example
%   s.FirstString = 'First string';
%   s.DoubleVector = [1 2]; %A double vector
%   s.TrueOrFalse = true; %A logical
%   [res,ok] = inputstruct(s,'mytitel') %Brings up dialog box.
%
% Limitations:
%   The function behaviours poorly when there a very many fields. In
%   these cases split it up by the usage of COPYFIELDS. Example
%   s.A = 1;
%   s.B = 2;
%   %etc
%   s = inputstruct(s,'First parameter set.')
%   temp = [];
%   temp.C = 3;
%   temp.D = 4;
%   s = copyfields(s,temp);
%
%   The function does only currently handles the following datatype
%   - double(s)
%   - chars or string arrays
%   - logicals
%
%   See also INPUTDLG, COPYFIELDS.

%Einar Heiberg 2003-12-04

if nargin==0
  error('Expected at least one input argument.');
end;

if nargin==1
  tit = inputname(1);
end;

if nargin>2
  error('Expected no more than two input arguments.');
end;

if length(s)>1
  error('Expected a struct, not an array of structs.');
end;

if nargout>2
  error('Expected only two output arguments.');
end;

%Get the fieldsnames
fields = fieldnames(s);
numfields = length(fields);

if numfields<1
  error('Input is not a struct.');
end;

%Set default answer as failed
if nargout==2
  varargout = cell(1,1);
  varargout{1} = false;
end;

%Create variables for dialogbox
prompt = cell(size(fields));
def = cell(size(fields));
lineno = zeros(numfields,1);
for loop=1:numfields
  prompt{loop} = [fields{loop} ':'];
  switch class(getfield(s,fields{loop}))
    case 'logical'
      temp = getfield(s,fields{loop});
      if temp
        def{loop} = 'true';
      else
        def{loop} = 'false';
      end;
    otherwise
      def{loop} = num2str(getfield(s,fields{loop}));
  end;
  lineno(loop) = max(size(def{loop},1),1);
end;

%Bring up the dialogbox
keystroke = popfrombuffer('KeyStroke');
if isempty(keystroke)
  answ = inputdlg(prompt,tit,lineno,def);
else
  if isequal(keystroke,'ok')
    answ = def;
  else
    error('Expected ''ok'' as keystroke.');
  end;
end;

if isempty(answ)
  %User pressed cancel, return the same.
  res = s;  
  return;
else
  %Parse the input
  res = [];
  ok = 1;
  for loop=1:numfields
    switch class(getfield(s,fields{loop}))
      case 'char'
        res = setfield(res,fields{loop},answ{loop});
      case 'double'
        %Convert to number
        [temp,tempok] = str2num(answ{loop});
        if tempok
          res = setfield(res,fields{loop},temp);
        else
          %Could not convert
          myfailed(dprintf('Could not convert %s to number.',fields{loop}));
          res = setfield(res,fields{loop},getfield(s,fields{loop}));
          ok = false;
        end;
      case 'logical'
        %Convert to number
        switch lower(answ{loop})
          case {'false','no'}
            res = setfield(res,fields{loop},false);
          case {'true','yes'}
            res = setfield(res,fields{loop},true);            
          otherwise
            %Could not convert
            myfailed(dprintf('Could not convert %s to logical.',fields{loop}));
            ok = false;
            res = setfield(res,fields{loop},getfield(s,fields{loop}));
        end;        
      otherwise
        error(sprintf('Input type %s not supported.',class(getfield(s,fields{loop}))));
    end;    
  end;
end;
varargout{1} = ok;