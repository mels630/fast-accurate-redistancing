function d = distanceToPoints(cpx,cpy,m,n,dx,dy)
% function d = distanceToPoints(cpx,cpy,m,n,dx,dy)

if(nargin < 3)
  disp(' Must enter cpx, cpy, and m=n as minimum number of inputs.')
  d = NaN;
  return;
end
if(nargin < 4)
  n = m;
end
if(nargin < 5)
  dx = 1/n;
end
if(nargin < 6)
  dy = dx;
end

d = ones(m,n) * max(dx*n,dy*m); d = d.^2;
x = (0:(n-1)) * dx;
y = (0:(m-1)) * dy;

[X,Y] = meshgrid(x,y);

for ii = 1:numel(cpx)
  % do squared distance, first
  d = min(d,(X-cpx(ii)).^2 + (Y-cpy(ii)).^2);
end
d = sqrt(d);
