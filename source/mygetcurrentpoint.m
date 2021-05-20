function [x,y] = mygetcurrentpoint(h)
%Returns CurrentPoint of an axes h. In case of macro or testing, returns
%current point from buffer.

%Einar Heiberg
global DATA

x = 0;
y = 0;

try
  if ~DATA.Testing
    p = get(h,'CurrentPoint');
    x = p(1);
    if numel(p)<3
      y = p(2);
    else
      y = p(3);
    end
  else
    if length(DATA.Buffer.CurrentPoint)<2
      p = get(h,'CurrentPoint');
      x = p(1);
      if numel(p)<3
        y = p(2);
      else
        y = p(3);
      end
    else
      x = popfrombuffer('CurrentPoint');
      y = popfrombuffer('CurrentPoint');
    end
  end
catch
  return;
end

if DATA.RecordMacro
  stri = get(h,'tag'); 
  namestri = inputname(1);
  macro_helper('put',sprintf('pushtobuffer(''CurrentPoint'',[%0.5g %0.5g]); %% %s %s',x,y,stri,namestri));
  macro_helper('switchorder'); %We need to store data in buffer before the callback
end
