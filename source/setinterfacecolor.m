function setinterfacecolor(handle)
%set background color and text color for all object in all interfaces according to selected background color
global DATA

if isa(DATA,'maingui')
  setcolorall(handle);
%   setcolormaingui(handle);
end

%------------------------
function setcolorall(handle)
global DATA

switch handle.Type
  case 'figure'
    set(handle,'Color',DATA.GUISettings.BackgroundColor);
    
  case 'axes'
    if ~(strcmp(handle.Tag,'printuipanel'))
      set(handle,'XColorMode','manual','GridColorMode','manual',...
        'GridColor',DATA.GUISettings.BackgroundColor,...
        'XColor',DATA.GUISettings.ForegroundColor,...
        'YColor',DATA.GUISettings.ForegroundColor,...
        'ZColor',DATA.GUISettings.ForegroundColor);
    end
    
  case {'uicontrol','uipanel','uibuttongroup'}
    %this handle is specific to Segment 3DP, should have ligth gray background
    if (strcmp(handle.Tag,'printuipanel'))
      h = findall(handle.printuipanel);
      for loop = 1:length(h)
        if any(contains({'figure','panel','listbox','edit','text','button'},handle.Tag))
          set(h(loop),'BackgroundColor',[0.94,0.94,0.94]);
          set(h(loop),'ForegroundColor',[0, 0, 0]);
        end
      end
      
    %this handle is specific to Segment 3DP, should have ligth gray background
    elseif (strcmp(handle.Tag,'ribbonuipanel'))
      set(handle,'BackgroundColor',[0.94,0.94,0.94]);
      set(handle,'ForegroundColor',[0, 0, 0]);
      
    %other panels can have user defined background color
    else
      set(handle,'BackgroundColor',DATA.GUISettings.BackgroundColor,...
        'ForegroundColor',DATA.GUISettings.ForegroundColor);
    end
end

%go through all objects
if any(contains({'figure','uipanel','uibuttongroup'},handle.Type))
  kids = handle.Children;
  
  for i = 1:numel(kids)
    try
      setcolorall(kids(i));
    catch me %#ok<NASGU>
    end
  end
end

%------------------------

%------------------------
function setcolor(handle)
%------------------------
global DATA

%find all objects in all interfaces
hstruct = get(handle);
fnames = fieldnames(hstruct);

if strcmp(hstruct.Type,'uimenu') || strcmp(hstruct.Type,'uicontextmenu')
  %Do nothing
else
  try
    if not(strcmp(hstruct.Type,'axes'))
      %set background color
      if not(isempty(intersect('Color',fnames)))
        disp([handle.Type,': ','Color'])
        set(handle,'Color',DATA.GUISettings.BackgroundColor);
      end
      if not(isempty(intersect('BackgroundColor',fnames)))
        disp([handle.Type,': ','BackgroundColor'])
        set(handle,'BackgroundColor',DATA.GUISettings.BackgroundColor);
      end
    end

    %set text color
    if not(isempty(intersect('ForegroundColor',fnames)))
      disp([handle.Type,': ','ForegroundColor'])
      set(handle,'ForegroundColor',DATA.GUISettings.ForegroundColor);
    end
    if strcmp(hstruct.Type,'axes')
      disp([handle.Type,': ','XColorMode'])
      disp([handle.Type,': ','GridColorMode'])
      disp([handle.Type,': ','GridColor'])
      set(handle,'XColorMode','manual');
      set(handle,'GridColorMode','manual');
      set(handle,'GridColor',DATA.GUISettings.BackgroundColor);
    end
    if not(isempty(intersect('XColor',fnames)))
      disp([handle.Type,': ','XColor'])
      set(handle,'XColor',DATA.GUISettings.ForegroundColor);
    end
    if not(isempty(intersect('YColor',fnames)))
      disp([handle.Type,': ','YColor'])
      set(handle,'YColor',DATA.GUISettings.ForegroundColor);
    end
    if not(isempty(intersect('ZColor',fnames)))
      disp([handle.Type,': ','ZColor'])
      set(handle,'ZColor',DATA.GUISettings.ForegroundColor);
    end
    
  catch
    disp('Error in set color');
  end
end

%go through all objects
kids = hstruct.Children;

for i = 1:numel(kids)
  try
    if not(strcmp(kids(i).Parent.Type,'axes'))
      if not(strcmp(kids(i).Tag,'printuipanel'))
        setcolor(kids(i));
      end
    end
  catch me %#ok<NASGU>
  end
end

%change background color and text color for specific objects
% lgcolor = [0.94 0.94 0.94];
% try
%    set(DATA.Handles.iconuipanel,'BackgroundColor',lgcolor);
% catch
% end
% try
%   set(DATA.Handles.permanentaxes,'Color',lgcolor);
% catch
% end
% try
%   set(DATA.Handles.ribbonaxes,'Color',lgcolor);
% catch
% end
% try
%   set(DATA.Handles.configaxes,'Color',lgcolor);
% catch
% end
