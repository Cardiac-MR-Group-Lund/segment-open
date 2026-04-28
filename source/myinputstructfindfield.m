function n = myinputstructfindfield(field,s)
%Helper function that goes together with myinputstruct. 
%Find which n has the name field, this is useful for @callback usage in
%myinputstruct
%
%N = MYINPUTSTRUCTFINDFIELD(FIELD,S)

%Einar Heiberg

%This can be written faster but fast enough for now
n = [];
for loop = 1:length(s)
  if isequal(field,s(loop).Field)
    n = loop;
  end
end
