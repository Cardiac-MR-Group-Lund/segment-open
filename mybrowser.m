function mybrowser(url)
%MYBROWSER(URL) Opens an URL platform independent.

%Jane Sj�gren
global DATA

command=['"' DATA.Pref.WebBrowser '" "' url '"'];
system(command);
flushlog;