function h = mymsgbox(stri,title,fighandle)
%MYMSGBOX Helper function, works as MSGBOX
%
%MYMSGBOX(STRI,TITLE,[ALIGNMENT])
%  STRI String to display
%  TITLE Title of messagebox
%  FIGHANDLE Optional figure handle indicating alignment
%
%See also MYWARNING, MYMSGBOX, MYADJUST, MYWAITBARSTART.

%Einar Heiberg

global DATA
if nargin<3
  fighandle=[];
end
if nargin<2
  title = '';
  fighandle=[];
end
  
stri = translation.dictionary(stri);
title = translation.dictionary(title);
try 
  if DATA.Pref.DoNotAsk
    %If do not ask then also do not display
    mydisp(dprintf('Message: %s\n',stri));
    return;
  end
catch %#ok<CTCH>
end

try
  mydisp(dprintf('Message: %s\n',stri));
catch %#ok<CTCH>
  dispstri = tools('maskpatientstrings',stri);
  fprintf('Message: %s\n\n',dispstri);
end

keystroke = popfrombuffer('KeyStroke');

h = [];
if isempty(keystroke)
  if contains(stri,'Your license does not include')
    createmode = struct('WindowStyle','modal','Interpreter', 'tex');
  else
    createmode = struct('WindowStyle','modal','Interpreter', 'none');
  end
  h = msgbox(stri,title,createmode);
  try
    htext = findobj(h, 'Type', 'Text');  %find text control in dialog
    set(htext,'Color',DATA.GUISettings.ForegroundColor,'FontSize',9,'Units','normalized');
  catch me
    mydispexception(me)
  end
  
   try
    set(h,'Color',DATA.GUISettings.BackgroundColor);
    kids = h.Children;
    for i=1:length(kids)
      try set(kids(i),'BackgroundColor',DATA.GUISettings.BackgroundColor);catch, end
      try set(kids(i),'Color',DATA.GUISettings.BackgroundColor);catch, end
      try set(kids(i),'ForegroundColor',DATA.GUISettings.ForegroundColor);catch, end
      try set(kids(i),'FontSize',12);catch, end
      try
        kidstyle = kids(i).Style;
        if strcmp(kidstyle,'pushbutton')
          set(kids(i),'Units','normalized','FontSize',9)
        end
      catch
      end
    end
  catch me
    mydispexception(me)
  end
  if length(stri) > 40
    ps = get(h,'Position');
    set(h,'Position',[ps(1) ps(2) ps(3)*1.2 ps(4)])
  end
  
  myadjust(h,fighandle);
  %   set(h,'windowstyle','modal');
  setupicon(h); %set up Segment icon
  uiwait(h);
  try
    macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from mymsgbox');
    macro_helper('switchorder'); %We need to store data in buffer before the callback
  catch %#ok<CTCH>
    %Do nothing if this fails. Most likely DATA variable is not
    %initialized.
  end
else
  if isequal(lower(keystroke),'ok')
    macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from mymsgbox');
    macro_helper('switchorder'); %We need to store data in buffer before the callback
    return;
  else
    error('Expected ''ok'' as keystroke');
  end
end
flushlog;
