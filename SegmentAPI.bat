@echo off
cd "C:\Program Files\Segment CMR"
set stri=SegmentAPI.exe
:buildstri
if \{%1\}==\{\} goto :call
set stri=%stri% %1
shift
goto buildstri
:call
%stri%