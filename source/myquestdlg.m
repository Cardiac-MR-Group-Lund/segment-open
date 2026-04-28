function buttonname = myquestdlg(dialogstr,buttonstrings,defaultstr,fighandle)
% dialog function with 3 answer's possibilities
arguments
  dialogstr {mustBeText} %dialog string presented to the user
  buttonstrings cell % cell array with all strings for the answer buttons
  defaultstr {mustBeText} = ''; % string for the default answer button
  fighandle = []; % figure the dialog will be placed over, if empty ''Segment GUI'' is used
end

global DATA %#ok<*GVMIS>

%color settings
if ~isempty(DATA)
  if DATA.Testing
    buttonname = popfrombuffer('KeyStroke');
    return
  end
  backgroundcolor = DATA.GUISettings.BackgroundColor;
  foregroundcolor = DATA.GUISettings.ForegroundColor;
  highlightcolor = 'blue';
else
  % use default colors
  backgroundcolor = get(0,'defaultUicontrolBackgroundColor');
  foregroundcolor = get(0,'defaultUicontrolForegroundColor');
  highlightcolor = foregroundcolor;
end

%settings for margins, font size etc.
fontsize = 13;
marginleftright = 10;
margintopbottom = 20;
distvertical = 15; %distance between lines
disthorizontal = 15; % distance between buttons

% setup invisible figure
fig = figure(...
  'Color',backgroundcolor,...
  'MenuBar','none',...
  'Visible','off',....
  'KeyPressFcn',@keypress,...
  'WindowStyle','modal',...% ensures that figure reacts to uiwait/uiresume
  'NumberTitle','off');

%Set up Segment icon
setupicon(fig);

% set uicontrol to be used to get the extent
tmpuicontrol = uicontrol( ...
  'parent',fig   , ...
  'Style','text', ...
  'FontSize',  fontsize,...
  'Visible','off'...
  );
% get necessary extent of the buttons
defaultstr = dprintf(defaultstr);
numbuttons = length(buttonstrings);
buttonwidth = zeros(numbuttons,1);

isvaliddefaultbutton = false;
for btn = 1:numbuttons
  btnstr = buttonstrings{btn};
  set(tmpuicontrol,'String',btnstr);
  buttonextent = tmpuicontrol.Extent;
  buttonwidth(btn) = buttonextent(3);% get minimal width for each button
  buttonheight = buttonextent(4);
  if strcmpi(btnstr,defaultstr)
    isvaliddefaultbutton = true;
    defaultbuttonind = btn;
  end
end
buttonwidthsingle = max(buttonwidth)+marginleftright; % button width + margin
buttonwidthall = numbuttons*buttonwidthsingle + ...
  (numbuttons-1)*disthorizontal;

% get necessary extent of the dialog text
if length(dialogstr) > 90 && ~contains(dialogstr, newline)
  tmpstr = mysplitstring(dialogstr);
  dialogstr = sprintf('%s\n%s',tmpstr{1},tmpstr{2});
end
dialogstr = dprintf(dialogstr);
set(tmpuicontrol,'String',dialogstr);
txtextent = tmpuicontrol.Extent;
txtwidth = txtextent(3);
txtheight = txtextent(4);

% calculate minimal figure width and height
figwidth = max(txtwidth,buttonwidthall) + 2*marginleftright;
figheight = txtheight + buttonheight + distvertical +2*margintopbottom;
[px, py, figwidth, figheight] = helperfunctions('getfigposition',figwidth,figheight,fighandle);
fig.Position = [px py figwidth figheight];

% place question text in the figure
txtx = marginleftright + (figwidth-2*marginleftright-txtwidth)/2;
txty = margintopbottom + distvertical + buttonheight;
uicontrol(...
  'parent',fig,...
  'Style','text',...
  'Position',[txtx txty txtwidth txtheight],...
  'String',dialogstr,...
  'HorizontalAlignment','left',...
  'ForegroundColor',foregroundcolor,...
  'BackgroundColor',backgroundcolor,...
  'Fontsize',fontsize...
  );

% place buttons in the figure
buttony = margintopbottom;
xposbuttons = zeros(1,3);
for btn = 1:numbuttons
  if btn == 1
    xposbuttons(1) = marginleftright + (figwidth-2*marginleftright-buttonwidthall)/2;
  else
    xposbuttons(btn) = xposbuttons(btn-1)+buttonwidthsingle+disthorizontal;
  end
  uicontrol(fig,...
    'Style','pushbutton',...
    'Position',[xposbuttons(btn) buttony buttonwidthsingle buttonheight],...
    'Stri',buttonstrings{btn},...
    'ForegroundColor',foregroundcolor,...
    'BackgroundColor',backgroundcolor,...
    'Callback',@button_Callback,...
    'Fontsize',fontsize...
    );
end
% place default frame over default button
if isvaliddefaultbutton
  xyoffset = 2;
  sizeoffset = 2*xyoffset;
  defaultx = xposbuttons(defaultbuttonind) - xyoffset;
  defaulty = buttony - xyoffset;
  defaultwidth = buttonwidthsingle+sizeoffset;
  defaultheight = buttonheight+sizeoffset;

  buttondefault = uipanel(fig, ...
    'HighlightColor', highlightcolor, ...
    'BorderType', 'etchedin', ...
    'units', 'pixels', ...
    'Position', [defaultx defaulty defaultwidth defaultheight],...
    'UserData',buttonstrings{defaultbuttonind});
  uistack(buttondefault, 'bottom');
end

fig.Visible = 'on';
uiwait;
if fig.isvalid
  buttonname = fig.UserData;
  close(fig);
else
  % figure was closed by user
  buttonname = '';
end

%------------------------------
function button_Callback(obj,~)
%------------------------------
% get name of the button object
fig = obj.Parent;
fig.UserData = obj.String;
uiresume(fig);

%------------------------------
function keypress(obj, evnt)
%------------------------------
% check pressed key
switch evnt.Key
  case {'return','space'}
    defaultbuttonobj = findobj(obj.Children,'Type','uipanel');
    if ~isempty(defaultbuttonobj)
      obj.UserData = defaultbuttonobj.UserData;
      uiresume(obj)
    end
end
