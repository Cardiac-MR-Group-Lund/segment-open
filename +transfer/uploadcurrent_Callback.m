function uploadcurrent_Callback
%uploadcurrent_Callback Uploads the current study (i.e matfile) to core
%lab features.

%Jonatan Wulcan, modified by Einar Heiberg

global DATA SET

if isempty(SET)
  myfailed('No data to upload');
  return
end

% Get the study id
if (not(isfield(SET(1), 'transferid'))) || isempty(SET(1).transferid)
  myfailed('Couldn''t find Transfer id. The current stack must originate from the Transfer system');
  return
end

idstr = SET(1).transferid;
spaceloc = regexp(idstr,'[ ]');
if length(spaceloc) ~= 1
  myfailed('Bad format in transferid');
  return
end
study_id = idstr(1:spaceloc-1);
project_url = idstr(spaceloc+1:end);

% Save matfile
filename = 'transfer-temp.mat';
matfile = [getpreferencespath filesep filename];
topatientdatabase = false;
h = msgbox('Storing file to disk. Please wait.');
fail = filemenu('saveallas_helper',getpreferencespath,filename,topatientdatabase);
try
  close(h);
catch    %#ok<CTCH>
end;

if fail
  myfailed('Couldn''t save. Aborting.');
  return;
end;

try
  
  try
    [~, user, pw, url] = transfer.readcredentials();
  catch me
    mydispexception(me);
    myfailed('Could not find any credential files.');
    return;
  end;
 
  match_url = url;
  
  %Remove http://
  if isequal(findstr(match_url,'http://'),1)  %#ok<*FSTR>
    match_url = match_url(8:end);
  end;
  
  %Remove www.
  if isequal(findstr(match_url,'www.'),1) 
    match_url = match_url(5:end);
  end;
    
  %Remove http://
  if isequal(findstr(project_url,'http://'),1)
    project_url = project_url(8:end);
  end;
  
  %Remove www.
  if isequal(findstr(project_url,'www.'),1)
    project_url = project_url(5:end);
  end;  
    
  while ~strcmp(match_url,project_url)
    myfailed(dprintf(...
      'Project URL does not match origin, You are trying to upload to wrong trial! \n\nUploading to %s and trial originated from %s.',...
      url,project_url));

    [~, user, pw, url, ~] = transfer.readcredentials();
    
    match_url = url;
    
    %Remove http://
    if isequal(findstr(match_url,'http://'),1)
      match_url = match_url(8:end);
    end;
  
    %Remove www.
    if isequal(findstr(match_url,'www.'),1)
      match_url = match_url(5:end);
    end;  
    
  end
catch %#ok<CTCH>
  return
end

try
  s = transfer.server(user, pw, url);
  transfer.sendmatfile(s, study_id, matfile);
catch e
  if strcmp(e.identifier, 'SEGMENT:ERROR')
    myfailed(e.message);
    return
  else
    rethrow(e);
  end
end

% Delete matfile
delete(matfile);

% Display message
old_dna = DATA.Pref.DoNotAsk;
DATA.Pref.DoNotAsk = false;
mymsgbox('Upload completed successfully!');
DATA.Pref.DoNotAsk = old_dna;
