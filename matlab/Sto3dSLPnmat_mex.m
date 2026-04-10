function [T] = Sto3dSLPnmat_mex(nt, ns, nt3, ns3, xt, nxt, xs, ws, if_self, T)
mex_id_ = 'Sto3dSLPnmat_mex(c i int[x], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[xx], c i double[xx], c i double[x], c i int[x], c io double[xx])';
[T] = stokes_mex(mex_id_, nt, ns, nt3, ns3, xt, nxt, xs, ws, if_self, T, 1, 1, 1, 1, 3, nt, 3, nt, 3, ns, ns, 1, nt3, ns3);
