% -----------------------
function area = stablepolyarea(x,y)
% -----------------------
% Computes the area of the polygon described by coordinates x and y in
% square-pixels in a numerically stable way.
%
% Algorithm: Discretized version Green's formula for area.
% Persson, Böiers: Analys i flera variabler. Studentlitteratur, Lund,
% Sweden 1988.
%
% Johannes Töger 2015-03-23
%

% Error checks
errmsgonecurve = 'stablepolyarea supports only one curve at a time.';

assert(ndims(x) <= 2, errmsgonecurve);
assert(ndims(y) <= 2, errmsgonecurve);

assert(size(x,1) == 1 || size(x,2) == 1, errmsgonecurve);
assert(size(y,1) == 1 || size(y,2) == 1, errmsgonecurve);

x = x(:);
y = y(:);

assert(length(x) == length(y), 'length(x) must equal length(y)');

% closes curve if not closed, no-op on already closed curves
x = [x; x(1)];
y = [y; y(1)];

% Compute area
x = x - mean(x);
y = y - mean(y);

dy = diff(y);
xm = conv(x, [1 1]/2, 'valid'); % dy is defined *between* points, so we need to approximate x between points also

dA = xm.*dy; % integrand in Green's formula

area = abs(stablesum(dA));


% -----------------------
function S = stablesum(x)
% -----------------------
% For a vector x, compute its sum S in a slightly more stable way.
% 1. Separate positive and negative parts
% 2. Sort on magnitude, ascending
% 3. Do the sum

x = x(:);

xpos = x(x >= 0);
xneg = x(x < 0);

S = kahansum(sort(xpos,1,'ascend')) + kahansum(sort(xneg,1,'descend')); % descend, since xneg < 0


% -----------------------
function S = kahansum(x)
% -----------------------
% Perform Kahan summation of the vector x.
% The intermediate variable c keeps track of the rounding error.

S = 0;
c = 0;

for ii = 1:length(x)
  y = x(ii) - c;
  t = S + y;
  c = (t - S) - y;
  S = t;
end