function [filename,pathname,filterindex,ok] = myuigetfile(varargin)
%[FILENAME, PATHNAME, FILTERINDEX,OK] = MYUIGETFILE(FILTERSPEC, TITLE)
%Corresponding to uigetfile, but also fixes macro recording and test
%scripts.

%Nils Lundahl

global DATA

if isa(DATA,'maingui')
  if DATA.Testing
    testing = DATA.Testing;
  else
    testing = false;
  end
  if DATA.RecordMacro
    recordmacro = DATA.RecordMacro;
  else
    recordmacro = false;
  end
else
  testing = false;
  recordmacro = false;
end

ok = false;
if testing
  %Take from buffer
  s = popfrombuffer('File');
  if ~isempty(s)
    ok = true;
    [pathname,file,ext] = fileparts(s);
    filename = [file ext];
    filterindex = 1;
  else
    error('No file selection in buffer.');
  end
else
  %Ask user
  [filename,pathname,filterindex] = uigetfile(varargin{:});
  if isequal(filename,0)
    return
  end
  ok = true;
end

if recordmacro
  file = fullfile(pathname,filename);
  macro_helper('put',sprintf('pushtobuffer(''File'',''%s'') ; %%add to buffer',file));
  macro_helper('switchorder');
end