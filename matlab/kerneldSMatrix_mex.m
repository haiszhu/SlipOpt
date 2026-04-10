function [A] = kerneldSMatrix_mex(p, n, n3, rxyz, nxyz, w, A)
mex_id_ = 'kerneldSMatrix_mex(c i int[x], c i int[x], c i int[x], c i double[xx], c i double[xx], c i double[x], c io double[xx])';
[A] = stokes_mex(mex_id_, p, n, n3, rxyz, nxyz, w, A, 1, 1, 1, 3, n, 3, n, n, n3, n3);
