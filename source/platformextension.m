function ext = platformextension(mext)
%Returns platform extension for segment / dcmtk binaries.

%Einar Heiberg

%case {'mexmaci','mexmaci64'}
%  ext = '_maci';
    
if nargin == 0
  mext = mexext();
end;

switch mext
  case {'mexw32','mexw64'}
    ext = '.exe';
  case 'mexa64'
    ext = '_linux';
  otherwise
    myfailed('This platform is not supported.');
    return;
end;
    