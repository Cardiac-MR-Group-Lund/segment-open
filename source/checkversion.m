function checkversion()
  if isdeployed()
    platform = mexext();
  else
    platform = 'source';
  end
  
  versionnbr = parseversionstr(changelog);
  try
    address = java.net.InetAddress.getLocalHost;
    IPaddress = char(address.getHostAddress);
    [online_version] = myurlread(sprintf( ...
      'http://www.medviso.com/getversion.php?platform=%s&version=%d&ip=%s', ...
      platform,versionnbr,IPaddress),750); %750 is 750ms timeout
  catch %#ok<CTCH>
    return
  end
  
  try 
    if parseversionstr(online_version) <= versionnbr
      return
    end
  catch %#ok<CTCH>
    return
  end
  
  c = mymenu('A new version is available', ...
    'Update now', ...
    'Ask again later');
  if c == 1
    mybrowser('http://www.medviso.com/products/segment/download/');
  end
end

function r = parseversionstr(vstr)
  tt = sscanf(vstr, '%d.%d R%d');
  if numel(tt) ~= 3
    error('SEGMENT:ERROR', 'Couldn''t parse version string');
  end
  r = tt(3);
end