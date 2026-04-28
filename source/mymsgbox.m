function h = mymsgbox(stri,title,fighandle,force)
%MYMSGBOX Helper function, works as MSGBOX
%
%MYMSGBOX(STRI,TITLE,[ALIGNMENT],[FORCE])
%  STRI String to display
%  TITLE Title of messagebox
%  FIGHANDLE Optional figure handle indicating alignment
%  FORCE Optional if true, then display regardless of DoNotAsk status in preferences
%
%See also MYWARNING, MYMSGBOX, MYADJUST, MYWAITBARSTART.

%Einar Heiberg

global DATA %#ok<*GVMIS> 

if nargin < 4
  force = false;
end

if nargin < 3
  fighandle=[];
end

if nargin < 2
  title = '';
  fighandle=[];
end
  
stri = translation.dictionary(stri);
title = translation.dictionary(title);
try 
  if (not(force) && DATA.Pref.DoNotAsk) || DATA.Autoloader
    %If do not ask then also do not display
    logdisp(dprintf('Message: %s',stri),false); %false = mask patient string
    return;
  end
catch %#ok<CTCH>
end

try
  logdisp(dprintf('Message: %s',stri),false); %false = mask patient string
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
  if ~isempty(DATA)
    try
      htext = findobj(h, 'Type', 'Text');  %find text control in dialog
      set(htext,'Color',DATA.GUISettings.ForegroundColor,'FontSize',9,'Units','normalized');
    catch me
      mydispexception(me)
    end

     try
      set(h,'Color',DATA.GUISettings.BackgroundColor);
      kids = h.Children;

      try set(kids,'BackgroundColor',DATA.GUISettings.BackgroundColor);catch, end
        try set(kids,'Color',DATA.GUISettings.BackgroundColor);catch, end
        try set(kids,'ForegroundColor',DATA.GUISettings.ForegroundColor);catch, end
        try set(kids,'FontSize',9);catch, end
        try
          for ind = 1:length(kids)
            kidstyle = {kids.Style};
            if strcmp(kidstyle,'pushbutton')
              set(kids(ind),'Units','normalized','FontSize',9)
            end
          end
        catch
        end
      
    catch me
      mydispexception(me)
     end
  end
  if length(stri) > 35
    ps = get(h,'Position');
    set(h,'Position',[ps(1) ps(2) ps(3)*1.1 ps(4)])
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
