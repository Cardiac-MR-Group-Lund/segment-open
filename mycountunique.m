function n = mycountunique(v,tol)
%MYCOUNTUNIQUE counts the number unique elements in the input vector.
%
%N = MYCOUNTUNIQUE(V,<tol>)
%
%N is the number of unique elements
%V is numeric input vector
%tol is numeric tolerance, if omitted then set 1e-6.
%
%See also matlab builtin unique.

%Einar Heiberg, 2008-04-02

if nargin<2
  tol = 1e-6;
end;

if isempty(v)
  n = 0;
  return;
end;

switch class(v)
  case {'double','single'}
    %Ensure dimension
    v = v(:);
    
    %Normalize vector
    maxv = max(v);
    minv = min(v);
    v = v-minv;
    %if abs(maxv-minv)>1e-6 EH:2009-04-27
    %  v = v/(maxv-minv);
    %end;
    
    %Multiply with 1/tol and round
    v = round(v*(1/tol));
    
    %Count elements
    n = length(unique(v));
  case 'cell'
    %check if the elements are numeric
    if nargin>1
      warning('Tolerance not used for alphanumeric string.');
    end;
    n = length(unique(v));
  otherwise
    error('Not implemented for other than numeric cells, double, and single.');
end;

