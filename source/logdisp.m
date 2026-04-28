function logdisp(message, keeppatientinfo, filename)
%LOGDISP(MESSAGE,KEEPPATIENTINFO,FILENAME) 
% Disp function to log with datestr 
% if keeppatientinfo=true, patient info will be kept in log message. 
% Specify filename if a special logfile should be used.
%
% If keeppatientinfo is not supplied, then it is set to true. Thus if the
% string could contain patientinfo and you do not want that, then call with
% keeppatientinfo = false

%Fanny Månefjord, Medviso, 2024

arguments
  message 
  keeppatientinfo = true; %this is default so we don't call mydisp unnecessarily often
  filename = ''
end

message = strrep(message, '\', '\\'); %Make it safe for sprintf
message = sprintf("%s %s \n", datestr(now,'yyyy-mm-dd HH:MM:SS '), message); 

if ~isempty(filename) %specific file chosen
  logfid = fopen(filename,'at');
  if keeppatientinfo
    %disp with date and potential patient info
    fprintf(logfid, message);
  else
    %disp with date and without patient info
    mydisp(message, logfid,false); %false = nowrap
  end
  fclose(logfid);
else %write to ordinary log file
  if keeppatientinfo
    %disp with date and potential patient info
    fprintf(message);
  else
    %disp with date and without patient info
    mydisp(message,[],false); %[] = logfid, false = nowrap
  end
end