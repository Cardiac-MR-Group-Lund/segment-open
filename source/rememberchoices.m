function rememberchoices(varargin)
%Function to help remember choices in user interfaces.
%
%Usage 
% 
% rememberchoices('Key1',value,'Key2',value,...) 
% 
% or
%
% rememberchoices({'Key1','Key2',...},outs) %where outs is a myinputstructstruct
%
%To recall a pre-stored choice use
%
%v = getchoice('Key1',default); %return stored choices if available, else returns default
%
%The keys must be valid fieldnames in Matlab.
%
%See also getchoice

%Einar Heiberg

global DATA %#ok<GVMIS> 

if nargin<2
  error('Need at least two input arguments');
end

arg1 = varargin{1};

%Get it
C = DATA.GUIChoices;

if isa(arg1,'cell')
  %--- Cell and outs

  if nargin~=2
    error('When called with cell array first only two arguments are expected')
  end

  arg2 = varargin{2};

  %Loop over and store
  for loop = 1:length(arg1)
    key = arg1{loop};
    if isfield(arg2,key)
      C.(key) = arg2.(key); %Store it
    else
      %Do nothing, should not occur unless bad call by user
    end
  end

else
  %--- Key value pairs

  n = nargin;
  if mod(n,2) > 0
    error('Expected key-value pairs')
  end

  %Loop over arguments
  for loop = 1:(n/2)

    %Extract
    keyind = (loop-1)*2+1;
    valueind = keyind+1;
    key = varargin{keyind};
    value = varargin{valueind};

    %Store
    C.(key) = value;

  end
end

%Store it back
DATA.GUIChoices = C;

