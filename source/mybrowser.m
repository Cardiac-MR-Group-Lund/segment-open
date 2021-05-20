function mybrowser(url,browser)
%MYBROWSER(URL) Opens an URL platform independent.

%Jane Sjögren
global DATA
if nargin == 1
  try
    browser = DATA.Pref.WebBrowser;
  catch
    mywarning(dprintf('Could not open default browser.\nPlease open your browser and go to\n%s',url))
    return
  end
end

try
  command=['"' browser '" "' url '"'];
  system(command);
  flushlog;
catch
  mywarning(dprintf('Could not open default browser.\nPlease open your browser and go to\n%s',url))
end
