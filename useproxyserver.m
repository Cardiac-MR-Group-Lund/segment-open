function useproxyserver(use)
%Turn on/off proxy server usage

%Einar Heiberg
global DATA

if nargin==0
  use = false;
end;

if use
  %Start and configure proxy server
  com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(true);
  com.mathworks.mlwidgets.html.HTMLPrefs.setProxyHost(DATA.Pref.ProxySettings.HostName);
  com.mathworks.mlwidgets.html.HTMLPrefs.setProxyPort(DATA.Pref.ProxySettings.Port);
  if isequal(DATA.Pref.ProxySettings.UserName,'') || isempty(DATA.Pref.ProxySettings.UserName)
    disp('Using proxy server without authentication.');
    com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxyAuthentication(false)
  else
    disp('Using proxy server with authentication.');
    com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxyAuthentication(true)
    com.mathworks.mlwidgets.html.HTMLPrefs.setProxyUsername(DATA.Pref.ProxySettings.UserName)
    com.mathworks.mlwidgets.html.HTMLPrefs.setProxyPassword(DATA.Pref.ProxySettings.Password)      
  end;
else
  %End proxy server
  
  com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(false);
end;
