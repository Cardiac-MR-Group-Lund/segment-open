function [a,b,c] = getangles(v1, v2, v3)
% Find angles of an image block for ParaView.
%
% Let Rx(a) be a rotation of a radians around the x-axis
%     Ry(b) be a rotation of b radians around the y-axis
%     Rz(c) be a rotation of c radians around the z-axis
%
% Let T(a,b,c) = Rz(c)*Ry(b)*Rx(a).
%
%
% Input:
% v1, v2, v3: Three orthonormal vectors in R^3.
%
% Output:
% a, b, c: Angles in radians such that
%   T(a,b,c)*[1 0 0]' = v1
%   T(a,b,c)*[0 1 0]' = v2
%   T(a,b,c)*[0 0 1]' = v3
%
% To convert to degrees, multiply by 180/pi.
%
% Johannes Töger 2012-09-28

T = [v1 v2 v3];

assert(norm(T*T' - eye(3)) < 1e-3, 'Vectors must form orthogonal basis.');

% Find possible a,b,c
a0 = atan(T(3,2)/T(3,3));
b0 = asin(-T(3,1));
c0 = atan(T(2,1)/T(1,1));

asol = [a0, a0+pi];
bsol = [b0, pi-b0];
csol = [c0, c0+pi];

if T(3,3) == 0
  % cos(b)*cos(a) = 0 => cos(a) = 0 or cos(b) = 0
  asol = [asol, pi/2, -pi/2];
  bsol = [bsol, pi/2, -pi/2];
end

if T(1,1) == 0
  % cos(b)*cos(c) = 0 => cos(b) = 0 or cos(c) = 0
  bsol = [bsol, pi/2, 3/2*pi];
  csol = [csol, pi/2, 3/2*pi];
end

% Wrap candidates to [-pi, pi]
asol = mod(asol+pi, 2*pi)-pi;
bsol = mod(bsol+pi, 2*pi)-pi;
csol = mod(csol+pi, 2*pi)-pi;

% Test all combinations of possible a, b, c
N = length(asol)*length(bsol)*length(csol);
avec = zeros(1,N);
bvec = zeros(size(avec));
cvec = zeros(size(avec));
errs = zeros(size(avec));

ctr = 1;
for ii = 1:length(asol)
  for jj = 1:length(bsol)
    for kk = 1:length(csol)
      a = asol(ii);
      b = bsol(jj);
      c = csol(kk);
      avec(ctr) = a;
      bvec(ctr) = b;
      cvec(ctr) = c;
      errs(ctr) = norm(R([a b c])-T, 2);

      ctr = ctr+1;
    end
  end
end

% Find one of the good solutions (there are at least 2)
ok = find(errs < 1e-3);
if isempty(ok)
  error('Could not find angles.');
end

abc(1) = avec(ok(1));
abc(2) = bvec(ok(1));
abc(3) = cvec(ok(1));

end

function R_ = R(abc)
% Evaluate composite rotation matrix Rz(c)*Ry(b)*Rx(a).

a = abc(1);
b = abc(2);
c = abc(3);

R_ = [cos(c)*cos(b), -sin(c)*cos(a)+cos(c)*sin(b)*sin(a), sin(c)*sin(a)+cos(c)*sin(b)*cos(a); sin(c)*cos(b), cos(c)*cos(a)+sin(c)*sin(b)*sin(a), -cos(c)*sin(a)+sin(c)*sin(b)*cos(a); -sin(b), cos(b)*sin(a), cos(b)*cos(a);];

end
