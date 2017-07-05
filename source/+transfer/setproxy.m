function setproxy()
  [h, p] = transfer.findproxy();
  if not(isequal(h, ''))
    com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(true);
    com.mathworks.mlwidgets.html.HTMLPrefs.setProxyHost(h);
    com.mathworks.mlwidgets.html.HTMLPrefs.setProxyPort(p);
    com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxyAuthentication(false);
  else
    com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(false);
  end
end