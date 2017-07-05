function [study, user, pw, url,pattern] = readcredentials()
%Read the credentials file.

%Jonatan Wulcan. Modified by Einar Heiberg
name = transfer.getcredsfile;
if isempty(name)
  error('SEGMENT:ERROR','Could not find any credentials file.');
end;

if iscell(name)
  %Found multiple
  
  %Loop over them to extract site name
  sitenames = cell(size(name)); 
  for loop = 1:length(sitenames)
    fid = fopen(name{loop},'rt');
    if isequal(fid,-1)
      error('SEGMENT:ERROR','Could not open credentials file.');
    end;
    stri = fgetl(fid);
    sitenames{loop} = stri;
    fclose(fid);
  end;
  m = mymenu('Select trial.',sitenames{:});
  
  if isequal(m,0)
    error('SEGMENT:ERROR','User aborted.');
  end;
  
  name = name{m};
end;

fid = fopen(name,'r');
if fid == -1
  error('SEGMENT:ERROR', sprintf('Couldn''t open %s. Please contact core lab for support.',name)); %#ok<SPERR>
end

%Read study
stri = fgetl(fid);
if isempty(stri)
  error('SEGMENT:ERROR','Could not read study name.');
end;
study = stri;

%Read username
stri = fgetl(fid);
if isempty(stri)
  error('SEGMENT:ERROR','Could not read user name.');
end;
user = stri;

%Read password
stri = fgetl(fid);
if isempty(stri)
  error('SEGMENT:ERROR','Could not read password.');
end;
pw = stri;

%Read url
stri = fgetl(fid);
if isempty(stri)
  error('SEGMENT:ERROR','Could not read URL.');
end;
url = stri;

%Read pattern
stri = fgetl(fid);
pattern = stri;

fclose(fid);