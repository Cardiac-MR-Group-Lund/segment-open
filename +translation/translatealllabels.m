function translatealllabels(handle,fromlanguage)
%TRANSLATEALLLABELS Recursively translate all labels of handle
global DATA
if nargin < 2
  fromlanguage = 'English';
end
if isa(DATA,'maingui') && isfield(DATA.Pref,'Language') && ~strcmp(DATA.Pref.Language,fromlanguage)
  dotranslate(handle,DATA.Pref.Language,fromlanguage);
end

%-------------------------------------------------
function dotranslate(handle,language,fromlanguage)
%-------------------------------------------------
hstruct = get(handle);
fnames = fieldnames(hstruct);
propertylist = {'Name','String','Label','Title','XLabel','YLabel','TooltipString'};
plist = intersect(fnames,propertylist);
for i = 1:numel(plist)
  p = plist{i};
  txt = hstruct.(p);
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
end

kids = hstruct.Children;
for i = 1:numel(kids)
  dotranslate(kids(i),language,fromlanguage);
end