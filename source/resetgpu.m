function resetgpu
%Resets GPU device and throughs an error message
%
%See also usegpu()

%Einar Heiberg

e = lasterror; %#ok<LERR>

if strfind(e.message,'Out of memory')
  stri1 = 'Out of memory on GPU. Reset GPU. If problem persists try to close 3D view.';
else
  stri1 = 'Unknown GPU error.';
end

gpuDevice([]); %deselects and clears its
g = gpuDevice;
stri2 = dprintf('Total GPU memory %0.5g GB. Available GPU memory %0.5g GB',g.TotalMemory/1e9,g.AvailableMemory/1e9);

stri = [stri1 ' ' stri2];

h = msgbox(stri,'GPU error');

mydisp(stri);

try
  pause(2);
  delete(h);  
catch
end
