function isuptodate = checkversion() %skiptranslation

global DATA %#ok<GVMIS> 

isuptodate = false;
if isdeployed()
  platform = mexext();
else
  platform = 'source';
end

[versionnbr, subversionnbr, releaseid] = helperfunctions('parseversionstr',changelog);
if strcmp(subversionnbr,'a')
  versionstr = sprintf('%d',versionnbr);
else
  versionstr = sprintf('%d%s',versionnbr,subversionnbr);
end
try
  address = java.net.InetAddress.getLocalHost;
  IPaddress = char(address.getHostAddress);
  %hostname = getenv('COMPUTERNAME'); % for windows
  [online_version] = myurlread(sprintf( ...
    'https://www.medviso.com/getversion.php?platform=%s&releaseid=%s&version=%d&versionstr=%s&ip=%s', ...
    platform,releaseid,versionnbr,versionstr,IPaddress),750); %750 is 750ms timeout
catch %#ok<CTCH>
  return
end

try
  [online_version, online_subversion] = parseonlineversionstr(online_version);
  if online_version < versionnbr || online_version == versionnbr && online_subversion <= subversionnbr
    isuptodate = true;
    return
  end
catch %#ok<CTCH>
  return
end

c = mymenu('A new version is available', ...
  'Update now', ...
  'Ask again later');
if c == 1
  if yesno('This will close Segment and open installation website. Continue?')
    DATA.CurrentTool = 'updating version';
    browser = DATA.Pref.WebBrowser;
    filemenu('quit_Callback');  %close Segment
    mybrowser('https://www.medviso.com/download2/',browser);
  end
end
end


function [ver,subver] = parseonlineversionstr(verstr)
tt = sscanf(verstr, '%d%s');
if ~(numel(tt) == 1 || numel(tt) == 2)
  error('SEGMENT:ERROR', 'Couldn''t parse online version string');
end
ver = tt(1);
if numel(tt) == 2
  subver = char(tt(2));
else
  subver = 'a';
end
end