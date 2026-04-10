function [A] = stokesKernelMat_mex(p, n, n3, rxyz, w, A)
mex_id_ = 'stokesKernelMat_mex(c i int[x], c i int[x], c i int[x], c i double[xx], c i double[x], c io double[xx])';
[A] = stokes_mex(mex_id_, p, n, n3, rxyz, w, A, 1, 1, 1, 3, n, n, n3, n3);
