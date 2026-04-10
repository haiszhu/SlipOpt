function [yn] = ynm_all_mex(n, np, nm, u, v, yn)
mex_id_ = 'ynm_all_mex(c i int[x], c i int[x], c i int[x], c i double[x], c i double[x], c io dcomplex[xx])';
[yn] = stokes_mex(mex_id_, n, np, nm, u, v, yn, 1, 1, 1, np, np, np, nm);
