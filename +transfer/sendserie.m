function sendserie(server, study, desc, files, wbfac)

serie_id = server.new_serie(desc, study);
wb = wbfac(numel(files), 'Uploading files...');
  
%First upload DICOMS
senddicoms(server, files, wb, serie_id);

%-----------------------------------------------
function senddicoms(server, files, wb, serie_id)
%-----------------------------------------------
%Send DICOMs (stored in files) to server

% Create fileinfos
fileinfos = struct('name', {}, 'hash', {});
for i=1:numel(files)
  fileinfos(i) = struct(...
    'name', files{i}, ...
    'hash', transfer.hashfile(files{i}));
end

% new dicomfiles
fileinfos = server.new_dicomfiles(fileinfos, serie_id);

% send files
sendfile = transfer.sendfile(@(data) server.update_dicomfile(data));
for i=1:numel(fileinfos)
  sendfile.addfile(fileinfos(i));
  wb.update();
end
sendfile.finalize();

% finalize files
server.finalize_dicomfile({fileinfos.id});