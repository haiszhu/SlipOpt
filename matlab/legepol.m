function [pol,der]=legepol(x,n)

pol=0;
der=0;
mex_id_ = 'legepol(c i double[x], c i int[x], c io double[x], c io double[x])';
[pol, der] = rotgrid_r2014a(mex_id_, x, n, pol, der, 1, 1, 1, 1);



