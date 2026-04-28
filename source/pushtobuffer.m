function pushtobuffer(f,o)
%PUSHTOBUFFER(F). Where F is the field to push to.
%
%Used to push keystrokes and values to user interfaces, used during
%testing. Push to the buffer before you call the interface in the
%testscript. Note that it acts as stack and the latest pushed is used
%first.
%
%Examples
% pushtobuffer('KeyStroke','yes'); %yes in yesno questions
% pushtobuffer('KeyStroke','ok'); %Ok in myinputstruct
% pushtobuffer('Dir',pathname); %push a foldername
% pushtobuffer('ListboxValue',2); % %push a listbox setting in myinputstruct
% pushtobuffer('EditString','Hello'); %Push a string to an editstring in myinputstring
% pushtobuffer('SliderValue',2); %Push slidervalue in uicontrol executes

%See also CLEARBUFFER, POPFROMBUFFER.

%Einar Heiberg

%#ok<*GVMIS>
global DATA

if nargin<2
  error('Expected two input arguments.');
end

buffer = DATA.Buffer.(f);
isacell = iscell(buffer);
if isacell
  buffer = [{o} buffer];
else
  buffer = [o buffer];
end

%Store it.
DATA.Buffer.(f) = buffer;
