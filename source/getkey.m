function key = getkey(evnt)
%--------------------------
%Decodes matlab evnt structure into a key string.
%
%Example of output strings:
%'A' (shift-a)
%'ctrl-shift-a'
%'ctrl-alt-f'

%Einar Heiberg

modifier = '';
if ~isempty(evnt(1).Modifier)
  for loop=1:length(evnt.Modifier)
    switch evnt(1).Modifier{loop}
      case 'control'
        modifier = [modifier 'ctrl']; %#ok<AGROW>
      case 'shift'
        modifier = [modifier 'shift']; %#ok<AGROW>
      case 'alt'
        modifier = [modifier 'alt']; %#ok<AGROW>
    end
    if loop<length(evnt.Modifier)
      modifier = [modifier '-']; %#ok<AGROW>
    end
  end
end

key = evnt.Key;

switch key
  case {'control','alt','shift'}
    key = '';
end

if ~isempty(modifier)
  if ~isempty(key)
    ch = evnt.Character;
    if ~isempty(key) && strcmp(ch,'+')
      key = ch;
    end
      key = [modifier '-' key];
  else
    key = modifier;
  end
end