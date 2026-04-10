function [A] = Sto3dSLPmat_mex(nt, ns, nt3, ns3, xt, xs, ws, if_self, A)
mex_id_ = 'Sto3dSLPmat_mex(c i int[x], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[xx], c i double[x], c i int[x], c io double[xx])';
[A] = stokes_mex(mex_id_, nt, ns, nt3, ns3, xt, xs, ws, if_self, A, 1, 1, 1, 1, 3, nt, 3, ns, ns, 1, nt3, ns3);
