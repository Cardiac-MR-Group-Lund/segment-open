function tooltip = maketooltipstring(tooltip,varargin)
hotkeystr = '';
numkeys = numel(varargin);
for n = 1:numkeys
  if n > 1 
    separator = '-';
  else
    separator = '';
  end
  hotkeystr = [hotkeystr,separator,varargin{n}]; %#ok<AGROW>
end
if numkeys > 0
  tooltip = sprintf('%s [%s]',tooltip,hotkeystr);
else
  tooltip = sprintf('%s',tooltip);
end
