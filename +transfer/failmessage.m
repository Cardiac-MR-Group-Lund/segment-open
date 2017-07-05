function failmessage
%This throws a fail message to instruct the site on what to do.

%Einar Heiberg
global DATA

% Display message
old_dna = DATA.Pref.DoNotAsk;
DATA.Pref.DoNotAsk = false; %Ensure that message box is displayed.
DATA.Buffer.KeyStroke = {}; %Ensure that message box is displayed.

mymsgbox(dprintf([...
  'Something went wrong! \n\n' ...
  'An error report will be generated after you have clicked OK on \n' ...
  'this message. Please copy errors from the error report and send \n' ...
  'them in an email to support@medviso.com.\n\n' ...
  'We will then look into the source of the problem and assist you as soon as possible.']));

flushlog;
if ispc
  dos(sprintf('notepad.exe "%s" &',DATA.LogFile));
else
  mymsgbox(dprintf('Log file for this run is %s',DATA.LogFile),'',DATA.GUI.Segment);
end;


mybrowser('mailto:support@medviso.com');

DATA.Pref.DoNotAsk = old_dna;