function adminrequirement

if isdeployed
  %check if can write file
  fid = fopen('delete.me','wt');
  if fid>0
    fclose(fid);
    delete('delete.me');
  else
    mymsgbox('This function may require that you are running the software as Administrator (not only be logged in with Administrator rights). To do so right click when starting software and select "Run as administrator".');
  end
end
