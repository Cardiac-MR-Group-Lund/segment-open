function z = isapproxequal(a,b,tol) 
%ISAPPROXEQUAL Same as isequal, but with tolerance for numeric arrays.

%Einar Heiberg

if nargin<3
  tol = 1e-4;
end
if strcmp(mexext,'mexw32')
  if nargin >= 3
    tol = max(tol,b*0.02); 
  else
    tol = b*0.02;
  end
  tol = abs(tol);
end
if ~isnumeric(a)
  z = isequal(a,b);
  return;
end
if isequal(size(a),size(b))
  if all(abs((a(:)-b(:)))<=tol(:))
    z = true;
  else
    z = false;
  end
else
  z = false;
end
end