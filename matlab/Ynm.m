function [Yn, Hn, Wn] = Ynm(n, m, u, v)
%YNM_MEX Fortran-backed spherical harmonics Yn using legacy associated-Legendre convention.

if nargin == 0
  error('Ynm_mex requires inputs (n,m,u,v)');
end

u = u(:);
v = v(:);
np = numel(u);

nm = 2*n + 1;
Yall = ynm_all_mex(n, np, nm, u, v, complex(zeros(np, nm)));
if isempty(m)
  Yn = Yall;
else
  if m > n || m < -n
    Yn = zeros(np, 1);
  else
    Yn = Yall(:, n+1+m);
  end
end

if nargout >= 2
  if exist('Ynm', 'file') == 2
    [~, Hn, Wn] = Ynm(n, m, u, v);
  else
    error('Ynm_mex:NoDerivatives', 'Hn/Wn are not implemented in mex path and Ynm.m is not on path.');
  end
end
