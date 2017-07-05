function Buffer = clearbuffer
%Clears buffers used for macropurposes
%
%See also POPFROMBUFFER, PUSHTOBUFFER.

%Einar Heiberg

Buffer = []; %Used by maketest and macro recordning
Buffer.KeyStroke = {}; 
Buffer.CurrentPoint = [];
Buffer.Dir = {}; %Reserved for folder selections.
Buffer.EditString = {}; %Edit boxes.
Buffer.File = {}; %Reserved for File selections.
Buffer.ListboxValue = []; %Reserved for valus from listboxes.
Buffer.Mymenu = []; %Menu selections from mymenu
Buffer.Number = []; %Entered numbers using mygetnumber
Buffer.SliderValue = []; %Get values from sliders.
Buffer.String = {}; %Reserved for entered strings. Replace inputdlg
Buffer.Structs = {}; %Reserved for values from inputstruct.
Buffer.Warnings = {}; %Reserved for warnings from security test
