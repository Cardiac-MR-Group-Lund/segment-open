function translatealllabels(handle,fromlanguage)
%TRANSLATEALLLABELS Recursively translate all labels of handle
global DATA
if nargin < 2
  fromlanguage = 'English';
end
if isa(DATA,'maingui') && isfield(DATA.Pref,'Language') && ~strcmp(DATA.Pref.Language,fromlanguage)
  translateall(handle,DATA.Pref.Language,fromlanguage);
end
%-------------------------------------------------
function translateall(handle,language,fromlanguage)
%-------------------------------------------------

switch handle.Type
  case 'figure'
    propertylist = {'Name'};
  case 'uimenu'
    propertylist = {'Label','Text'};
  case 'uipanel'
    propertylist = {'Title'};
  case 'axes'
    propertylist = {'Title','Xlabel','YLabel','Label'};
  case 'uicontrol'
    propertylist = {'String','Tooltip'};
  case 'uibuttongroup'
    propertylist = {'Title','Tooltip'};
  otherwise
    propertylist = {};
end
if ~isempty(propertylist)
  alltext = get(handle,propertylist);
  for ind = 1:numel(alltext)
    % get and translate only the translatable properties
    txt = alltext{ind};
    if ~isempty(txt)
      if ischar(txt)
        set(handle,propertylist{ind},translation.dictionary(txt,language,fromlanguage));
      elseif iscell(txt)
        for j = 1:numel(txt)
          if ischar(txt{j}) && ~isempty(txt{j})
            txt{j} = translation.dictionary(txt{j},language,fromlanguage);
          end
        end
        set(handle,propertylist{ind},txt);
      end
    end
  end
end
kids = get(handle, 'Children');
for i = 1:numel(kids)
  translateall(kids(i),language,fromlanguage);
end

%-------------------------------------------------
function dotranslate(handle,language,fromlanguage)
%-------------------------------------------------
propertylist = {'Name','String','Label','Title','XLabel','YLabel','TooltipString','Text'};
for i = 1:numel(propertylist)
  p = propertylist{i};
  try
    
    txt = get(handle,p);
    if ~isempty(txt)
      if ischar(txt)
        set(handle,p,translation.dictionary(txt,language,fromlanguage));
      elseif iscell(txt)
        for j = 1:numel(txt)
          if ischar(txt{j}) && ~isempty(txt{j})
            txt{j} = translation.dictionary(txt{j},language,fromlanguage);
          end
        end
        set(handle,p,txt);
      end
    end
  catch
  end
end

kids = handle.Children;
for i = 1:numel(kids)
  dotranslate(kids(i),language,fromlanguage);
end


% hstruct = get(handle);
% fnames = fieldnames(hstruct);
% propertylist = {'Name','String','Label','Title','XLabel','YLabel','TooltipString','Text'};
% plist = intersect(fnames,propertylist);
% for i = 1:numel(plist)
%   p = plist{i};
%   txt = hstruct.(p);
%   if ~isempty(txt)
%     if ischar(txt)
%       set(handle,p,translation.dictionary(txt,language,fromlanguage));
%     elseif iscell(txt)
%       for j = 1:numel(txt)
%         if ischar(txt{j}) && ~isempty(txt{j})
%           txt{j} = translation.dictionary(txt{j},language,fromlanguage);
%         end
%       end
%       set(handle,p,txt);
%     end
%   end
% end
% 
% kids = hstruct.Children;
% for i = 1:numel(kids)
%   dotranslate(kids(i),language,fromlanguage);
% end

