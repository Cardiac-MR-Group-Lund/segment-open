function blockfig(fig)
%Adds figure to list of figures that needs to be closed before image stack
%can be switched.

%Einar Heiberg
global DATA

DATA.BlockingFigs = union(DATA.BlockingFigs,fig); %Also removes