function didfail = sendstudy(server, series, name, wbfac)
%Send a study, series is a cell array struct of files

if(numel(series) == 0)
  error('SEGMENT:ERROR', 'No series to send');
end

%Set up study object
study = server.new_study('N/A', name);
  
%Loop over series to send
didfail = false;
for i=1:numel(series)
  
  %Send serie
  try
    transfer.sendserie(...
      server, study, 'N/A', series{i}, wbfac);
  catch me
    disp(sprintf('Could not send serie %d',i)); %#ok<DSPS>
    mydispexception(me);
    didfail = true;
  end;
    
  %Update progress bar
  transfer.progressbar(0.3+0.7*i/numel(series),...
    sprintf('Sending serie %d/%d',i,numel(series)));
end

try
  server.finalize_study(study);
catch me
  disp('Could not finalize study');
  mydispexception(me);
  didfail = true;
end;
