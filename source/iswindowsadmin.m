function out = iswindowsadmin
% returns true if this user is in admin role

% Written by Jelena
persistent outvalue

if isempty(outvalue)
  if ispc
    %Original code
    wi = System.Security.Principal.WindowsIdentity.GetCurrent();
    wp = System.Security.Principal.WindowsPrincipal(wi);
    outvalue = wp.IsInRole(System.Security.Principal.WindowsBuiltInRole.Administrator);
  else
    outvalue = false;
  end
end
out = outvalue;