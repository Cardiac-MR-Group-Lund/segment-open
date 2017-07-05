function [s,ind]=removechars(stri)
%REMOVECHARS Remove numeric chars from string stri

%Einar Heiberg

ind = (stri=='0');
ind = ind|(stri=='1');
ind = ind|(stri=='2');
ind = ind|(stri=='3');
ind = ind|(stri=='4');
ind = ind|(stri=='5');
ind = ind|(stri=='6');
ind = ind|(stri=='7');
ind = ind|(stri=='8');
ind = ind|(stri=='9');
ind = ind|(stri=='.');
s = stri(ind);