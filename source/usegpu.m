function use = usegpu()
%Checks if to use GPU, displays message while GPU is "warming up".
%
%See also resetgpu()

%Einar Heiberg

global DATA

use = DATA.Pref.GPU.Use;

if ~use
  return;
end

if isempty(DATA.GPU) || (~DATA.GPU.Hot)  
  
  %Fire it up
  DATA.GPU.Hot = true;
  h = waitbar(0,'Initalizing GPU engine.');  
  temp = gpuArray(single(0)); %#ok<NASGU>
  clear temp;
  close(h);
  
end