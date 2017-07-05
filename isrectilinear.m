function out = isrectilinear(vec)
% Function to determine whether vector vec is rectilinear or not.
% Used to see origin of SET field TimeVector

out = std(diff(vec)) < 0.005;