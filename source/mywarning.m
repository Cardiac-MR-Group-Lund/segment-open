function mywarning(stri,fighandle)
%MYWARNING(STRI,HANDLE) Displays a warning message
%  User needs to click OK to discard message.
%
%See also  MYWARNINGNOBLOCK
%
%  If DATA.Silent then message is just displayed in the console
%  and exectution continues as if user had pressed ok.
%

%%%% If changing this file, please also change mywarning no block, since
%%%% they should be identical except for the uiwait(h); line.
global DATA %#ok<*GVMIS> 

%Einar Heiberg
% make sure that the string we display is in English to be able to read and
% understand our log file
try
  logstr = translation.dictionary(stri,'English',DATA.Pref.Language);
catch
  logstr = stri;
end

%For autoloader, add warning to warninglist instead of showing it
try
  if ~isempty(DATA) && DATA.Autoloader
    fprintf('Error: %s\n',logstr);
    warningfunctions('addwarning', logstr); 
    return;
  end
catch
  %Do nothing as if error then depends on DATA not initialised
end
%If silent then just write to log and return
try
  if ~isempty(DATA) && DATA.Silent
    fprintf('Attention: %s\n',logstr); 
    return;
  end
catch
  %Do nothing as if error then depends on DATA not initialised
end

stri = translation.dictionary(stri);

try
  if ~isempty(DATA) && DATA.Pref.DoNotAsk
    %Just print warning message in window if DoNotAsk mode.
    fprintf('%s Attention: %s\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),logstr);
    return;
  end
catch %#ok<CTCH>
end

keystroke = popfrombuffer('KeyStroke');

if isempty(keystroke)
  fprintf('%s Attention: %s\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),logstr);
  h = warndlg(stri,dprintf('Attention'),'invisible');
  
  %set color
  if ~isempty(DATA)
    try
      htext = findobj(h, 'Type', 'Text');  %find text control in dialog
      set(htext,'Color',DATA.GUISettings.ForegroundColor,'FontSize',9,'Units','normalized');
    catch me
      mydispexception(me)
    end
  end

  %set icon
  setupicon(h);
  
  if ~isempty(DATA)
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
  end
  if length(stri) > 40
    ps = get(h,'Position');
    set(h,'Position',[ps(1) ps(2) ps(3)*1.2 ps(4)])
  end
  %Adjust horizontally
  if nargin > 1
    myadjust(h,fighandle);
  else
    myadjust(h);
  end
  set(h,'visible','on');
  flushlog;
  uiwait(h); %This is the only difference to mywarningnoblock.
  try
    macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from mywarning');
    macro_helper('switchorder'); %We need to store data in buffer before the callback
  catch %#ok<CTCH>
    %Ignore if failed, likely due to DATA is cleared or not present
  end
else
  fprintf('Attention: %s\n\n',logstr);  
  if isequal(lower(keystroke),'ok')
    pushtobuffer('Warnings',stri);
    macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from mywarning');
    macro_helper('switchorder'); %We need to store data in buffer before the callback
    return;
  else
    error(sprintf('Expected ''ok'' as keystroke, got %s',keystroke)); %#ok<SPERR>
  end
end
flushlog;
