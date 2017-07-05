function sendmatfile(server, study, matfile)
  % Get the hash and data
  f = fopen(matfile, 'r');
  data = fread(f, inf, '*uint8');
  fclose(f);
  h = transfer.hashfile(matfile);
  
  % Create the file
  file_id = server.new_matfile(study, h);
  
  % Send the acctual file
  pos = 1;
  chunksize = 524288;
  hw = waitbar(0,'Please wait uploading .mat file.');
  while(pos <= numel(data))
    server.update_matfile(file_id, pos-1, ...
      transfer.base64encode( ...
      data( ...
      pos:min([pos+chunksize-1 numel(data)]) ...
      ) ...
      ));
    pos = pos + chunksize;
    waitbar(pos/numel(data),hw);
  end
  close(hw);
  server.finalize_matfile(file_id);
end