function mainhandles = appendemptyhandles(mainfig)
% append all handles that are present in all other figures but not present
% in the current one. Used whan loading a GUI

%Jelena Bock

mainhandles = guihandles(mainfig);
% get current handle names
fnames = fieldnames(mainhandles);
% get all handle names
getallhandlenames; % this places variable "allhandlenames" into woking space

missingnames = setdiff(allhandlenames,fnames);

if ~isempty(missingnames)
  for h = 1:length(missingnames)
    field = missingnames{h}; 
    mainhandles.(field) = [];  
  end
end