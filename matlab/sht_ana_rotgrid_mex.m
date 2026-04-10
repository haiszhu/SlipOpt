function [shc] = sht_ana_rotgrid_mex(p, np, nc, f, nshc, shc)
mex_id_ = 'sht_ana_rotgrid_mex(c i int[x], c i int[x], c i int[x], c i double[xx], c i int[x], c io dcomplex[xx])';
[shc] = stokes_mex(mex_id_, p, np, nc, f, nshc, shc, 1, 1, 1, np, nc, 1, nshc, nc);
