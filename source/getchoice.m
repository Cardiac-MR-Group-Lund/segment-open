function c = getchoice(key,default)
%Function to extract choice from 
%
%To recall a pre-stored choice use
%
%v = getchoice('Key1',default); %return stored choices if available, else returns default
%
%The keys must be valid fieldnames in Matlab.
%
%See also rememberchoices

%Einar Heiberg

global DATA %#ok<GVMIS> 

%Get it
C = DATA.GUIChoices;

if isfield(C,key)
  c = C.(key); %Return it as it is stored
else
  c = default;
end
