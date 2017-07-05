function a = copyfields(a,b)
%COPYFIELDS Copy fields from a struct to another struct.
%  A = COPYFIELDS(A,B) Copies the fields from the struct 
%  B to the struct A.
%
%  Limitations: The function assumes that A & B are 1-by-1 
%  structures.
%
%  See also FIELDNAMES, GETFIELD, ISFIELD, SETFIELD, RMFIELD.

%Einar Heiberg 2003-12-08

temp = fieldnames(b);
for loop=1:length(temp)
  a = setfield(a,temp{loop},getfield(b,temp{loop}));
end;

