function macro_helper(varargin)
%Helper function to track function history and to generate macros
%
%A dream would be to make this function actually really record macros.
%Left to fix:
%- clicked coords should be recorded.
%- generate keystroke command for mymsgbox,...
%- folder selection,... => implement myuigetdir,...
%
%Problems for true recording
%1) non accessible GUI handles, maybe solved by mygui class??
%2) folder/filename selections

%Einar Heiberg

global DATA

%return; %temporarily disable

%--- Check if start/stop
if nargin>0
  
  %start clause
  if isequal(varargin{1},'_start')
    DATA.RecordMacro = true;

    DATA.Macro = cell(1,1000);
    DATA.MacroN = 1;
    DATA.Macro{DATA.MacroN} = dprintf('%% Starting on %s',datestr(now));
    DATA.MacroN = DATA.MacroN+1;
    
    mydisp('--- Starting macro recording ---');
    return;
  end;
  
  %stop clause
  if isequal(varargin{1},'_stop')
    
    %End file nicely
    outfid = fopen('macro.txt','wt');
    if isequal(outfid,-1)
      myfailed('Could not open macro file for writing.');
      return;
    end;
    for loop=1:(DATA.MacroN-1)
      fprintf(outfid,'%s\n',DATA.Macro{loop});
    end;
    fprintf(outfid,'%% Ending\n');
    fclose(outfid);
    
    mydisp('--- Stop macro recording ---');        
    DATA.RecordMacro = false;
    DATA.Macro = [];
    DATA.MacroN = 0;
    return;
  end;  
  
  if isequal(varargin{1},'put')
    if DATA.RecordMacro
      DATA.Macro{DATA.MacroN} = varargin{2};
      DATA.MacroN = DATA.MacroN+1;
      return;
    end;
  end;
  
  if isequal(varargin{1},'switchorder')
    if DATA.RecordMacro
      switchorder;
    end;
  end;
end;

try
  if isempty(DATA) || ~DATA.RecordMacro
    return; %done
  end;
catch
  return
end;

st = dbstack;

if length(st)==2
  %Only when stack not deeper than two since then internal calls
  stri = '';
  fname = varargin{1};
  for loop=2:nargin
    switch class(varargin{loop})
      case 'char'
        tempstri = sprintf('%s',varargin{loop});
      case {'double','single'}
        if numel(varargin{loop})>1
          tempstri = ['[' sprintf('%d ',size(varargin{loop})) ']'];
        else
          tempstri = sprintf('%0.5g',varargin{loop});
        end;
    end;

    if loop>2
      stri = [stri ',' tempstri]; %#ok<AGROW>
    else
      stri = tempstri;
    end;

  end;

  switch nargin
    case 0
      DATA.Macro{DATA.MacroN} = removeext(st(2).file);
      DATA.MacroN = DATA.MacroN+1;
    case 1
      DATA.Macro{DATA.MacroN} = sprintf('%s(''%s'');',removeext(st(2).file),fname);      
      DATA.MacroN = DATA.MacroN+1;      
    otherwise
      DATA.Macro{DATA.MacroN} = sprintf('%s(''%s'',''%s'');',removeext(st(2).file),fname,stri);
      DATA.MacroN = DATA.MacroN+1;
  end;

  %For handles
  h = gco;
  tag = get(h,'tag');
  if not(isempty(tag))
    switch tag
      case 'uicontrol'
        switch get(h,'style')
          case 'listbox'
            v = get(h,'value');
            DATA.Macro{DATA.MacroN} = sprintf('%%%s set to value %d\n',tag,v);
            switchorder;
            DATA.MacroN = DATA.MacroN+1;
          case 'slider'
            v = get(h,'value');
            DATA.Macro{DATA.MacroN} = sprintf('%%%s set to value %d\n',tag,v);
            switchorder;
            DATA.MacroN = DATA.MacroN+1;
        end;
      otherwise
    end;
  end;
    
end;

%------------------------------
function stri = removeext(stri)
%------------------------------
stri = stri(1:(end-2));

%-------------------
function switchorder
%-------------------
%Switch order for two elements on macro stack

global DATA

if DATA.MacroN<3
  myfailed('Cannot switch place on top element.');
end;

temp = DATA.Macro{DATA.MacroN-2};
DATA.Macro{DATA.MacroN-2} = DATA.Macro{DATA.MacroN-1};
DATA.Macro{DATA.MacroN-1} = temp;
