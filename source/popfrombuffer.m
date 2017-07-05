function v = popfrombuffer(f)
%V = POPFROMBUFFER(F). Where F is the field to pop from. When buffer is
%empty an empty value is returned.
%
%See also CLEARBUFFER, PUSHTOBUFFER.

%Einar Heiberg

global DATA

try
  buffer = DATA.Buffer.(f);
  isacell = iscell(buffer);
  if isempty(buffer)
    if isacell
      v = {}; %Empty => return empty
    else
      v = [];
    end;
  else
    if isacell
      v = buffer{1}; %Return first element
    else
      v = buffer(1); %Return first element
    end;

    %Update stored buffer
    if isacell
      %Cell buffers
      if length(buffer)>1
        if length(buffer)==2
          buffer = buffer(2);
        else
          buffer = buffer(2:end);
        end
      else
        buffer = {};
      end;
    else
      %vector buffers
      buffer = buffer(2:end);      
      if isempty(buffer)
        buffer = []; %make is nice format
      end;        
    end;
    
    %Store it.
    DATA.Buffer.(f) = buffer;
  end;
catch
  v = [];
end;