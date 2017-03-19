function [v,cp] = redist(u,width,flag,dx,dy,dz)
% function [v,cp] = redist(u,width,flag,dx,dy,dz)
%   calls redistC mex-file to redistance u out to width using method
%   flag.
%  If nargin < 2, width = 9
%  If nargin < 3, flag = 1
%  If nargin < 4, dx = 1/n
%  If nargin < 5, dy = 1/m
%  If nargin < 6, dz = 1/k
%  
%  Flag:  0 == taxicab metric
%         1 == First-order accurate redistancing (e.g. fast marching method of Sethian; 
%                equivalently, method of Tsitsiklis)
%         2 == Directional optimization (biquadratic interpolation)
%         3 == Directional optimization (bicubic interpolation)
%  Returns the updated signed distance function (v) 
%  If flag > 1, also returns closest point functions (cpx,cpy).

if(nargin < 6)
  dz = 1/size(u,3);
end
if(nargin < 5)
  dy = 1/size(u,1);
end
if(nargin < 4)
  dx = 1/size(u,2);
end
if(nargin < 3)
  flag = 1;
end
if(nargin < 2)
  width = 9;
end
if(nargin < 1)
  disp('Not a proper function call!');
  v = NaN;
  cp = [];
  return;
end
size3 = size(u,3);

width = min(width,max(size(u)));
[v,cpx,cpy,cpz] = redistMex(u,width,flag,size3,dx,dy,dz);

if(flag < 2)
  cp = [];
elseif(size3 == 1)
  cp = {cpx,cpy};
else
  cp = {cpx,cpy,cpz};
end
v = reshape(v,size(u));
