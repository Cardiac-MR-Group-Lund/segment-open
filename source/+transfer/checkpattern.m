function [ok,msg] = checkpattern(pattern,newname)
%Check if the new name is valid according to the pattern.

ok = true;
msg = '';

if isequal(newname,pattern)
  ok = false;
  msg = 'You need to change the default value. Upload aborted.';
  return;
end;

if ~isequal(length(newname),length(pattern))
  ok = false;
  msg = sprintf('Length of new name and pattern need to be same. Length of pattern is %d and length of new name is %d',...
    length(pattern),length(newname)); 
  return;
end;

for loop = 1:length(pattern)

  %If it is a dash then it should be a dash
  if isequal(pattern(loop),'-') && ~isequal(newname(loop),'-') 
    ok = false;
    msg = sprintf('Need to have a dash in position %d.',loop); 
    return;
  end;
  
  %If it is a colon then it should be a dash
  if isequal(pattern(loop),':') && ~isequal(newname(loop),':') 
    ok = false;
    msg = sprintf('Need to have a colon in position %d.',loop); 
    return;
  end;  
  
  %If it is a semi-colon then it should be a 
  if isequal(pattern(loop),';') && ~isequal(newname(loop),';') 
    ok = false;
    msg = sprintf('Need to have a semi-colon in position %d.',loop); 
    return;
  end;  
  
  %If it is a point then it should be a 
  if isequal(pattern(loop),'.') && ~isequal(newname(loop),'.') 
    ok = false;
    msg = sprintf('Need to have a point in position %d.',loop); 
    return;
  end;      

  %If it is a space then it should be a 
  if isequal(pattern(loop),' ') && ~isequal(newname(loop),' ') 
    ok = false;
    msg = sprintf('Need to have a space in position %d.',loop); 
    return;
  end;      
  
  %If it is an underscore then it should be a 
  if isequal(pattern(loop),'_') && ~isequal(newname(loop),'_') 
    ok = false;
    msg = sprintf('Need to have an underscorein position %d.',loop); 
    return;
  end;      
  
  %It is a upper case letter => should be equal
  if isequal(pattern(loop),upper(pattern(loop))) 
    if ~isequal(pattern(loop),newname(loop)) 
      ok = false;
      msg = sprintf('Letter number %d needs to be a ''%s''',loop,pattern(loop)); 
    end;
  end;

  %If it is a lower case n then it should be numeric
  if isequal(pattern(loop),'n') %#ok<PROP>
    if ~ismember(newname(loop),'0123456789')
      ok = false;
      msg = sprintf('New name in position %d needs to be numeric.',loop);
      return;
    end;
  end;

  %If it is a lower case x then it should be non numeric
  if isequal(pattern(loop),'x') %#ok<PROP>
    if ismember(newname(loop),'0123456789')
      ok = false;
      msg = sprintf('New name in position %d needs to be non-numeric.',loop); 
      return;
    end;
  end;

end;