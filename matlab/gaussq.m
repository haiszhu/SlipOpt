function [t,w] = gaussq(kind,n,alpha,beta,kpts,endpts)

b = zeros(n,1);
t = zeros(n,1);
w = zeros(n,1);

mex_id_ = 'gaussq(c i int[x], c i int[x], c i double[x], c i double[x], c i int[x], c i double[x], c io double[], c io double[], c io double[])';
[b, t, w] = rotgrid_r2014a(mex_id_, kind, n, alpha, beta, kpts, endpts, b, t, w, 1, 1, 1, 1, 1, 2);


