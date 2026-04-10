subroutine stokesKernelMat_mex(p, n, n3, rxyz, w, a)
  use stokes_mod, only: stokesKernelMat
  implicit none
  integer, intent(in) :: p, n, n3
  real(8), intent(in) :: rxyz(3, n)
  real(8), intent(in) :: w(n)
  real(8), intent(inout) :: a(n3, n3)
  if (n3 /= 3*n) stop 'stokesKernelMat_mex: n3 mismatch'
  call stokesKernelMat(p, n, rxyz, w, a)
end subroutine stokesKernelMat_mex

subroutine Sto3dSLPmat_mex(nt, ns, nt3, ns3, xt, xs, ws, if_self, a)
  use stokes_mod, only: Sto3dSLPmat
  implicit none
  integer, intent(in) :: nt, ns, nt3, ns3, if_self
  real(8), intent(in) :: xt(3, nt), xs(3, ns), ws(ns)
  real(8), intent(inout) :: a(nt3, ns3)
  if (nt3 /= 3*nt) stop 'Sto3dSLPmat_mex: nt3 mismatch'
  if (ns3 /= 3*ns) stop 'Sto3dSLPmat_mex: ns3 mismatch'
  call Sto3dSLPmat(nt, ns, xt, xs, ws, if_self, a)
end subroutine Sto3dSLPmat_mex

subroutine Sto3dSLPnmat_mex(nt, ns, nt3, ns3, xt, nxt, xs, ws, if_self, t)
  use stokes_mod, only: Sto3dSLPnmat
  implicit none
  integer, intent(in) :: nt, ns, nt3, ns3, if_self
  real(8), intent(in) :: xt(3, nt), nxt(3, nt), xs(3, ns), ws(ns)
  real(8), intent(inout) :: t(nt3, ns3)
  if (nt3 /= 3*nt) stop 'Sto3dSLPnmat_mex: nt3 mismatch'
  if (ns3 /= 3*ns) stop 'Sto3dSLPnmat_mex: ns3 mismatch'
  call Sto3dSLPnmat(nt, ns, xt, nxt, xs, ws, if_self, t)
end subroutine Sto3dSLPnmat_mex

subroutine kerneldSMatrix_mex(p, n, n3, rxyz, nxyz, w, a)
  use stokes_mod, only: kerneldSMatrix
  implicit none
  integer, intent(in) :: p, n, n3
  real(8), intent(in) :: rxyz(3, n)
  real(8), intent(in) :: nxyz(3, n)
  real(8), intent(in) :: w(n)
  real(8), intent(inout) :: a(n3, n3)
  if (n3 /= 3*n) stop 'kerneldSMatrix_mex: n3 mismatch'
  call kerneldSMatrix(p, n, rxyz, nxyz, w, a)
end subroutine kerneldSMatrix_mex

subroutine sht_ana_rotgrid_mex(p, np, nc, f, nshc, shc)
  use stokes_mod, only: sht_ana_rotgrid
  implicit none
  integer, intent(in) :: p, np, nc, nshc
  real(8), intent(in) :: f(np, nc)
  complex(8), intent(inout) :: shc(nshc, nc)
  call sht_ana_rotgrid(p, np, nc, f, nshc, shc)
end subroutine sht_ana_rotgrid_mex

subroutine ynm_all_mex(n, np, nm, u, v, yn)
  use stokes_mod, only: ynm_all_rotgrid
  implicit none
  integer, intent(in) :: n, np, nm
  real(8), intent(in) :: u(np), v(np)
  complex(8), intent(inout) :: yn(np, nm)
  if (nm /= 2*n+1) stop 'ynm_all_mex: nm mismatch'
  call ynm_all_rotgrid(n, np, u, v, yn)
end subroutine ynm_all_mex
