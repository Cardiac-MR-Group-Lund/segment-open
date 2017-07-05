function mystubfailed
%Displays generic error message from stubs
global DATA

myfailed(dprintf('This module is not available in this version of %s.',DATA.ProgramName));
