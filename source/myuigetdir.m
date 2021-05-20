function [pathname,ok] = myuigetdir(pathname, titlestr)
%[PATHNAME,OK] = MYUIGETDIR(PATHNAME,TITLESTRI)
%Corresponding to uigetdir, but also fixes macro recording and test
%scripts.

%Einar Heiberg

global DATA

if nargin<1
  error('Expected at least one input arguments (prompt).');
end

if nargin<2
  titlestr = '';
end
%translation
titlestr = translation.dictionary(titlestr);

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
  pathname = popfrombuffer('Dir');
  if ~isempty(pathname)
    ok = true;
  else
    error('No path selection in buffer.');
  end
else
  %Ask user
  if DATA.isSiemensVersion
    % in Open apps version the path is fixed
    pathname = '\\tsclient\c';
  end
  pathname = uigetdir(pathname, titlestr);
  if isequal(pathname,0)
    return;
  end
  ok = true;
end

if recordmacro
  macro_helper('put',sprintf('pushtobuffer(''Dir'',''%s'') ; %%add to buffer',pathname));
  macro_helper('switchorder'); %We need to store data in buffer before the callback
end