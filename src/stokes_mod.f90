module stokes_mod
  use sht_mod, only: gl_nodes, shc_expand, shc_shrink
  implicit none
  private

  public :: rot_mat
  public :: sht_ana_rotgrid
  public :: ynm_all_rotgrid
  public :: Sto3dSLPmat
  public :: Sto3dSLPnmat
  public :: stokesKernelMat
  public :: kerneldSMatrix

  integer, save :: sphtrans_p_cached = -1
  real(8), save, allocatable :: ctheta_sph_cache(:)
  real(8), save, allocatable :: wht_sph_cache(:)
  real(8), save, allocatable :: ynms_sph_cache(:,:,:)
  complex(8), save, allocatable :: wsave_sph_cache(:)

  interface
    subroutine rotviaprojf90(beta, nterms, m1, m2, mpole, lmp, mpout, lmpn)
      real(8), intent(in) :: beta
      integer, intent(in) :: nterms, m1, m2, lmp, lmpn
      complex(8), intent(in) :: mpole(0:lmp,-lmp:lmp)
      complex(8), intent(inout) :: mpout(0:lmpn,-lmpn:lmpn)
    end subroutine rotviaprojf90

    subroutine legepols(x, n, pols)
      real(8), intent(in) :: x
      integer, intent(in) :: n
      real(8), intent(out) :: pols(n+1)
    end subroutine legepols

    subroutine sphtrans_cmpl_lege_init(nterms, nphi, ntheta, ctheta, whts, ynms, wsave)
      integer, intent(in) :: nterms, nphi, ntheta
      real(8), intent(out) :: ctheta(ntheta), whts(ntheta)
      real(8), intent(out) :: ynms(0:nterms,0:nterms,ntheta/2+1)
      complex(8), intent(out) :: wsave(4*nphi+15)
    end subroutine sphtrans_cmpl_lege_init

    subroutine sphtrans_fwd_cmpl(nterms, mpole, nphi, ntheta, fgrid, ctheta, whts, ynms, wsave)
      integer, intent(in) :: nterms, nphi, ntheta
      complex(8), intent(out) :: mpole(0:nterms,-nterms:nterms)
      complex(8), intent(inout) :: fgrid(nphi,ntheta)
      real(8), intent(in) :: ctheta(ntheta), whts(ntheta)
      real(8), intent(in) :: ynms(0:nterms,0:nterms,ntheta/2+1)
      complex(8), intent(in) :: wsave(4*nphi+15)
    end subroutine sphtrans_fwd_cmpl

    subroutine sphtrans_cmpl(nterms, mpole, nphi, ntheta, fgrid, ctheta, ynms, wsave)
      integer, intent(in) :: nterms, nphi, ntheta
      complex(8), intent(in) :: mpole(0:nterms,-nterms:nterms)
      complex(8), intent(out) :: fgrid(nphi,ntheta)
      real(8), intent(in) :: ctheta(ntheta)
      real(8), intent(in) :: ynms(0:nterms,0:nterms,ntheta/2+1)
      complex(8), intent(in) :: wsave(4*nphi+15)
    end subroutine sphtrans_cmpl
  end interface

contains

  subroutine rot_mat(p, np, beta, r)
    integer, intent(in) :: p, np
    real(8), intent(in) :: beta
    real(8), intent(out) :: r(np, np)

    integer :: nu, nv, nshc, nf
    integer :: i, j, k, n, m, col1, col2, mm, c, ind, sgn, iu, iv
    real(8) :: sqrt4pi
    real(8), allocatable :: f(:,:)
    real(8), allocatable :: ratio(:)
    complex(8), allocatable :: shc(:,:), shc_rot(:,:)
    complex(8), allocatable :: shc_col(:,:), shc_full(:,:), shc_out(:,:)
    complex(8), allocatable :: mmat(:,:), mflip(:,:), mout(:,:), mtmp(:,:)
    complex(8), allocatable :: frot(:,:)
    complex(8), allocatable :: fgrid(:,:), mpole(:,:)

    nu = p + 1
    nv = 2 * p
    nshc = nu * nu
    nf = nu * (2*p + 1)

    if (np /= nu * (2*p)) stop 'rot_mat: np mismatch'

    allocate(f(np, np), ratio(2*p+1))
    allocate(shc(nshc, np), shc_rot(nshc, np))
    allocate(shc_col(nshc,1), shc_full(nf,1), shc_out(nshc,1))
    allocate(mmat(nu,2*p+1), mflip(nu,2*p+1), mout(nu,2*p+1), mtmp(nu,2*p+1))
    allocate(frot(np, np))
    allocate(fgrid(nv, nu), mpole(0:p,-p:p))
    sqrt4pi = sqrt(4.0d0 * acos(-1.0d0))

    f = 0.0d0
    do i = 1, np
      f(i,i) = 1.0d0
    end do

    if (p /= sphtrans_p_cached) then
      if (allocated(ctheta_sph_cache)) deallocate(ctheta_sph_cache)
      if (allocated(wht_sph_cache)) deallocate(wht_sph_cache)
      if (allocated(ynms_sph_cache)) deallocate(ynms_sph_cache)
      if (allocated(wsave_sph_cache)) deallocate(wsave_sph_cache)
      allocate(ctheta_sph_cache(nu), wht_sph_cache(nu))
      allocate(ynms_sph_cache(0:p,0:p,nu/2+1))
      allocate(wsave_sph_cache(4*nv+15))
      call sphtrans_cmpl_lege_init(p, nv, nu, ctheta_sph_cache, wht_sph_cache, ynms_sph_cache, wsave_sph_cache)
      sphtrans_p_cached = p
    end if
    do c = 1, np
      do iu = 1, nu
        do iv = 1, nv
          fgrid(iv, iu) = dcmplx(f(iu + (iv-1)*nu, c), 0.0d0)
        end do
      end do

      mpole = (0.0d0, 0.0d0)
      call sphtrans_fwd_cmpl(p, mpole, nv, nu, fgrid, ctheta_sph_cache, wht_sph_cache, ynms_sph_cache, wsave_sph_cache)

      shc_full(:,1) = (0.0d0, 0.0d0)
      do m = -p, p
        ind = (m + p) * nu + 1
        do n = 0, p
          if (abs(m) <= n) then
            if (m <= 0) then
              sgn = (-1)**n
            else
              sgn = (-1)**(n+m)
            end if
            shc_full(ind+n,1) = dble(sgn) * sqrt4pi * mpole(n, m)
            if (n == p .and. abs(m) == p) then
              shc_full(ind+n,1) = 0.5d0 * shc_full(ind+n,1)
            end if
          end if
        end do
      end do
      call shc_shrink(p, 1, shc_full, shc_out)
      shc(:,c) = shc_out(:,1)
    end do

    do i = 1, 2*p+1
      mm = i - (p+1)
      if (mm <= 0) then
        if (mod(abs(mm),2) == 0) then
          ratio(i) = 1.0d0
        else
          ratio(i) = -1.0d0
        end if
      else
        ratio(i) = 1.0d0
      end if
    end do

    do k = 1, np
      shc_col(:,1) = shc(:,k)
      call shc_expand(p, 1, shc_col, shc_full)
      mmat = reshape(shc_full(:,1), [nu, 2*p+1])

      do n = 0, p
        col1 = p + 1 - n
        col2 = p + 1 + n
        do j = col1, col2
          mmat(n+1,j) = mmat(n+1,j) * ratio(j)
        end do
      end do

      do j = 1, 2*p+1
        mflip(:,j) = mmat(:, 2*p+2-j)
      end do

      call rotviaprojf90(beta, p, p, p, mflip, p, mout, p)

      do j = 1, 2*p+1
        mtmp(:,j) = mout(:, 2*p+2-j)
      end do

      do n = 0, p
        col1 = p + 1 - n
        col2 = p + 1 + n
        do j = col1, col2
          mtmp(n+1,j) = mtmp(n+1,j) / ratio(j)
        end do
      end do

      shc_full(:,1) = reshape(mtmp, [nf])
      call shc_shrink(p, 1, shc_full, shc_out)
      shc_rot(:,k) = shc_out(:,1)
    end do

    do c = 1, np
      shc_col(:,1) = shc_rot(:,c)
      call shc_expand(p, 1, shc_col, shc_full)

      mpole = (0.0d0, 0.0d0)
      do m = -p, p
        ind = (m + p) * nu + 1
        do n = 0, p
          if (abs(m) <= n) then
            if (m <= 0) then
              sgn = (-1)**n
            else
              sgn = (-1)**(n+m)
            end if
            mpole(n, m) = shc_full(ind+n,1) / (dble(sgn) * sqrt4pi)
          end if
        end do
      end do

      call sphtrans_cmpl(p, mpole, nv, nu, fgrid, ctheta_sph_cache, ynms_sph_cache, wsave_sph_cache)
      do iu = 1, nu
        do iv = 1, nv
          frot(iu + (iv-1)*nu, c) = dcmplx(dble(fgrid(iv, iu)), 0.0d0)
        end do
      end do
    end do
    r = dble(frot)

    deallocate(f, ratio, shc, shc_rot, shc_col, shc_full, shc_out)
    deallocate(mmat, mflip, mout, mtmp, frot)
    deallocate(fgrid, mpole)
  end subroutine rot_mat


  subroutine sht_ana_rotgrid(p, np, nc, f, nshc, shc)
    integer, intent(in) :: p, np, nc, nshc
    real(8), intent(in) :: f(np, nc)
    complex(8), intent(out) :: shc(nshc, nc)

    integer :: nu, nv, nf
    integer :: c, iu, iv, m, n, ind, sgn
    real(8) :: sqrt4pi
    complex(8), allocatable :: shc_full(:,:), shc_out(:,:), fgrid(:,:), mpole(:,:)

    nu = p + 1
    nv = 2 * p
    nf = (2*p + 1) * nu

    if (np /= nu*nv) stop 'sht_ana_rotgrid: np mismatch'
    if (nshc /= nu*nu) stop 'sht_ana_rotgrid: nshc mismatch'

    if (p /= sphtrans_p_cached) then
      if (allocated(ctheta_sph_cache)) deallocate(ctheta_sph_cache)
      if (allocated(wht_sph_cache)) deallocate(wht_sph_cache)
      if (allocated(ynms_sph_cache)) deallocate(ynms_sph_cache)
      if (allocated(wsave_sph_cache)) deallocate(wsave_sph_cache)
      allocate(ctheta_sph_cache(nu), wht_sph_cache(nu))
      allocate(ynms_sph_cache(0:p,0:p,nu/2+1))
      allocate(wsave_sph_cache(4*nv+15))
      call sphtrans_cmpl_lege_init(p, nv, nu, ctheta_sph_cache, wht_sph_cache, ynms_sph_cache, wsave_sph_cache)
      sphtrans_p_cached = p
    end if

    allocate(shc_full(nf, nc), shc_out(nshc, nc), fgrid(nv, nu), mpole(0:p,-p:p))
    shc_full = (0.0d0, 0.0d0)
    sqrt4pi = sqrt(4.0d0 * acos(-1.0d0))

    do c = 1, nc
      do iu = 1, nu
        do iv = 1, nv
          fgrid(iv, iu) = dcmplx(f(iu + (iv-1)*nu, c), 0.0d0)
        end do
      end do

      mpole = (0.0d0, 0.0d0)
      call sphtrans_fwd_cmpl(p, mpole, nv, nu, fgrid, ctheta_sph_cache, wht_sph_cache, ynms_sph_cache, wsave_sph_cache)

      do m = -p, p
        ind = (m + p) * nu + 1
        do n = 0, p
          if (abs(m) <= n) then
            if (m <= 0) then
              sgn = (-1)**n
            else
              sgn = (-1)**(n+m)
            end if
            shc_full(ind+n, c) = dble(sgn) * sqrt4pi * mpole(n, m)
            if (n == p .and. abs(m) == p) then
              shc_full(ind+n, c) = 0.5d0 * shc_full(ind+n, c)
            end if
          else
            shc_full(ind+n, c) = (0.0d0, 0.0d0)
          end if
        end do
      end do
    end do

    call shc_shrink(p, nc, shc_full, shc_out)
    shc = shc_out
    deallocate(shc_full, shc_out, fgrid, mpole)
  end subroutine sht_ana_rotgrid

  subroutine ynm_all_rotgrid(n, np, u, v, yn)
    integer, intent(in) :: n, np
    real(8), intent(in) :: u(np), v(np)
    complex(8), intent(out) :: yn(np, 2*n+1)

    integer :: i, m
    real(8), allocatable :: y(:,:)
    real(8) :: phase

    allocate(y(np, 2*n+1))
    call assoc_lege(n, np, u, y)

    do m = -n, n
      do i = 1, np
        phase = dble(m) * v(i)
        yn(i, n+1+m) = dcmplx(y(i, n+1+m) * cos(phase), y(i, n+1+m) * sin(phase))
      end do
    end do

    deallocate(y)
  end subroutine ynm_all_rotgrid

  subroutine assoc_lege(n, ntheta, theta, y)
    integer, intent(in) :: n, ntheta
    real(8), intent(in) :: theta(ntheta)
    real(8), intent(out) :: y(ntheta, 2*n+1)

    integer :: i, m, l
    real(8) :: x, sx, pmm, pmmp1, pll, pnm, pkm2, pkm1
    real(8) :: nmfac, pi

    pi = 4.0d0 * atan(1.0d0)
    y = 0.0d0

    do i = 1, ntheta
      x = cos(theta(i))
      sx = sqrt(max(0.0d0, 1.0d0 - x*x))

      pmm = 1.0d0
      do m = 0, n
        if (m > 0) pmm = dble(2*m-1) * sx * pmm

        if (n == m) then
          pnm = pmm
        else
          pmmp1 = dble(2*m+1) * x * pmm
          if (n == m+1) then
            pnm = pmmp1
          else
            pkm2 = pmm
            pkm1 = pmmp1
            pll = pkm1
            do l = m+2, n
              pll = (dble(2*l-1) * x * pkm1 - dble(l+m-1) * pkm2) / dble(l-m)
              pkm2 = pkm1
              pkm1 = pll
            end do
            pnm = pll
          end if
        end if

        nmfac = sqrt((dble(2*n+1)/(4.0d0*pi)) * exp(log_gamma_ratio(n-m, n+m)))
        y(i, n+1+m) = dble((-1)**m) * nmfac * pnm
        if (m > 0) then
          y(i, n+1-m) = dble((-1)**m) * y(i, n+1+m)
        end if
      end do
    end do
  end subroutine assoc_lege

  real(8) function log_gamma_ratio(a, b)
    integer, intent(in) :: a, b
    log_gamma_ratio = log_gamma(dble(a+1)) - log_gamma(dble(b+1))
  end function log_gamma_ratio


  subroutine stokesKernelMat(p, n, rxyz, w, a)
    integer, intent(in) :: p, n
    real(8), intent(in) :: rxyz(3, n)
    real(8), intent(in) :: w(n)
    real(8), intent(out) :: a(3*n, 3*n)

    integer :: nu, nv, np
    integer :: i, j, k, c, src, t, ii, jj, indg
    real(8) :: pi, theta_i, d
    real(8), allocatable :: xg(:), gwt(:), pols(:), wt_theta(:), ywt(:)
    real(8), allocatable :: utheta(:), wsph(:), w0(:)
    real(8), allocatable :: x(:), y(:), z(:)
    real(8), allocatable :: rall(:,:,:), rmat(:,:)
    integer, allocatable :: ind(:)
    real(8), allocatable :: xx(:), yy(:), zz(:), wk(:), invrho(:)
    real(8), allocatable :: gx(:), gy(:), gz(:), g(:), lv(:), rowv(:)
    real(8), allocatable :: g11(:,:), g22(:,:), g33(:,:), g12(:,:), g13(:,:), g23(:,:)

    nu = p + 1
    nv = 2 * p
    np = nu * nv
    pi = 4.0d0 * atan(1.0d0)

    if (n /= np) stop 'stokesKernelMat: n mismatch'

    allocate(xg(nu), gwt(nu), pols(p+1), wt_theta(nu), ywt(np))
    allocate(utheta(nu), wsph(np), w0(np))
    allocate(x(np), y(np), z(np))
    allocate(rall(np,np,nu), rmat(np,np))
    allocate(ind(np))
    allocate(xx(np), yy(np), zz(np), wk(np), invrho(np))
    allocate(gx(np), gy(np), gz(np), g(np), lv(np), rowv(np))
    allocate(g11(np,np), g22(np,np), g33(np,np), g12(np,np), g13(np,np), g23(np,np))

    x = rxyz(1,:)
    y = rxyz(2,:)
    z = rxyz(3,:)

    call gl_nodes(nu, xg, gwt)
    xg = -xg

    do i = 1, nu
      call legepols(xg(i), p, pols)
      d = sum(pols)
      theta_i = acos(xg(i))
      wt_theta(i) = (pi / dble(p)) * gwt(i) * d / cos(theta_i / 2.0d0)
      utheta(i) = theta_i
    end do

    do k = 1, nv
      do i = 1, nu
        ywt((k-1)*nu + i) = wt_theta(i) / (8.0d0*pi)
        wsph((k-1)*nu + i) = sin(utheta(i))
      end do
    end do

    w0 = w / wsph

    do j = 1, nu
      call rot_mat(p, np, utheta(j), rmat)
      rall(:,:,j) = rmat
    end do

    g11 = 0.0d0; g22 = 0.0d0; g33 = 0.0d0
    g12 = 0.0d0; g13 = 0.0d0; g23 = 0.0d0

    do k = 1, nv
      t = 0
      do c = 1, nv
        src = mod(c - (k-1) - 1 + nv, nv) + 1
        do i = 1, nu
          t = t + 1
          ind(t) = i + (src-1)*nu
        end do
      end do

      do j = 1, nu
        indg = j + nu*(k-1)

        do ii = 1, np
          do jj = 1, np
            rmat(ii,jj) = rall(ind(ii), ind(jj), j)
          end do
        end do

        xx = matmul(rmat, x)
        yy = matmul(rmat, y)
        zz = matmul(rmat, z)
        wk = matmul(rmat, w0) * wsph

        do i = 1, np
          invrho(i) = 1.0d0 / sqrt((xx(i)-x(indg))**2 + (yy(i)-y(indg))**2 + (zz(i)-z(indg))**2)
          g(i) = ywt(i) * wk(i) * invrho(i)
          gx(i) = (xx(i)-x(indg)) * invrho(i)
          gy(i) = (yy(i)-y(indg)) * invrho(i)
          gz(i) = (zz(i)-z(indg)) * invrho(i)
        end do

        lv = g * (1.0d0 + gx*gx)
        rowv = matmul(lv, rmat)
        g11(indg,:) = rowv

        lv = g * (1.0d0 + gy*gy)
        rowv = matmul(lv, rmat)
        g22(indg,:) = rowv

        lv = g * (1.0d0 + gz*gz)
        rowv = matmul(lv, rmat)
        g33(indg,:) = rowv

        lv = g * gx * gy
        rowv = matmul(lv, rmat)
        g12(indg,:) = rowv

        lv = g * gx * gz
        rowv = matmul(lv, rmat)
        g13(indg,:) = rowv

        lv = g * gy * gz
        rowv = matmul(lv, rmat)
        g23(indg,:) = rowv
      end do
    end do

    do i = 1, np
      do j = 1, np
        a(3*i-2, 3*j-2) = g11(i,j)
        a(3*i-2, 3*j-1) = g12(i,j)
        a(3*i-2, 3*j  ) = g13(i,j)

        a(3*i-1, 3*j-2) = g12(i,j)
        a(3*i-1, 3*j-1) = g22(i,j)
        a(3*i-1, 3*j  ) = g23(i,j)

        a(3*i  , 3*j-2) = g13(i,j)
        a(3*i  , 3*j-1) = g23(i,j)
        a(3*i  , 3*j  ) = g33(i,j)
      end do
    end do

    deallocate(xg, gwt, pols, wt_theta, ywt)
    deallocate(utheta, wsph, w0, x, y, z)
    deallocate(rall, rmat, ind)
    deallocate(xx, yy, zz, wk, invrho, gx, gy, gz, g, lv, rowv)
    deallocate(g11, g22, g33, g12, g13, g23)
  end subroutine stokesKernelMat

  subroutine Sto3dSLPmat(nt, ns, xt, xs, ws, if_self, a)
    integer, intent(in) :: nt, ns, if_self
    real(8), intent(in) :: xt(3, nt), xs(3, ns), ws(ns)
    real(8), intent(out) :: a(3*nt, 3*ns)

    integer :: i, j
    integer :: i1, i2, i3, j1, j2, j3
    real(8) :: pi, scale
    real(8) :: dx, dy, dz, r2, ir, irw, ir3w
    real(8) :: g11, g22, g33, g12, g13, g23

    pi = 4.0d0 * atan(1.0d0)
    scale = 1.0d0 / (8.0d0 * pi)
    a = 0.0d0

    do i = 1, nt
      i1 = i
      i2 = nt + i
      i3 = 2*nt + i
      do j = 1, ns
        j1 = j
        j2 = ns + j
        j3 = 2*ns + j

        if (if_self /= 0 .and. nt == ns .and. i == j) cycle

        dx = xt(1, i) - xs(1, j)
        dy = xt(2, i) - xs(2, j)
        dz = xt(3, i) - xs(3, j)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2 <= 0.0d0) cycle

        ir = 1.0d0 / sqrt(r2)
        irw = ir * ws(j)
        ir3w = ir*ir*ir * ws(j)

        g11 = irw + dx*dx*ir3w
        g22 = irw + dy*dy*ir3w
        g33 = irw + dz*dz*ir3w
        g12 = dx*dy*ir3w
        g13 = dx*dz*ir3w
        g23 = dy*dz*ir3w

        a(i1, j1) = scale * g11
        a(i1, j2) = scale * g12
        a(i1, j3) = scale * g13
        a(i2, j1) = scale * g12
        a(i2, j2) = scale * g22
        a(i2, j3) = scale * g23
        a(i3, j1) = scale * g13
        a(i3, j2) = scale * g23
        a(i3, j3) = scale * g33
      end do
    end do
  end subroutine Sto3dSLPmat

  subroutine Sto3dSLPnmat(nt, ns, xt, nxt, xs, ws, if_self, t)
    integer, intent(in) :: nt, ns, if_self
    real(8), intent(in) :: xt(3, nt), nxt(3, nt), xs(3, ns), ws(ns)
    real(8), intent(out) :: t(3*nt, 3*ns)

    integer :: i, j
    integer :: i1, i2, i3, j1, j2, j3
    real(8) :: pi, scale
    real(8) :: dx, dy, dz, r2, ir3w, dnxir2, fac
    real(8) :: a11, a22, a33, a12, a13, a23

    pi = 4.0d0 * atan(1.0d0)
    scale = -3.0d0 / (4.0d0 * pi)
    t = 0.0d0

    do i = 1, nt
      i1 = i
      i2 = nt + i
      i3 = 2*nt + i
      do j = 1, ns
        j1 = j
        j2 = ns + j
        j3 = 2*ns + j

        if (if_self /= 0 .and. nt == ns .and. i == j) cycle

        dx = xt(1, i) - xs(1, j)
        dy = xt(2, i) - xs(2, j)
        dz = xt(3, i) - xs(3, j)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2 <= 0.0d0) cycle

        ir3w = ws(j) / (r2 * sqrt(r2))
        dnxir2 = (dx*nxt(1, i) + dy*nxt(2, i) + dz*nxt(3, i)) / r2
        fac = scale * dnxir2

        a11 = dx*dx*ir3w
        a22 = dy*dy*ir3w
        a33 = dz*dz*ir3w
        a12 = dx*dy*ir3w
        a13 = dx*dz*ir3w
        a23 = dy*dz*ir3w

        t(i1, j1) = fac * a11
        t(i1, j2) = fac * a12
        t(i1, j3) = fac * a13
        t(i2, j1) = fac * a12
        t(i2, j2) = fac * a22
        t(i2, j3) = fac * a23
        t(i3, j1) = fac * a13
        t(i3, j2) = fac * a23
        t(i3, j3) = fac * a33
      end do
    end do
  end subroutine Sto3dSLPnmat


  subroutine kerneldSMatrix(p, n, rxyz, nxyz, w, a)
    integer, intent(in) :: p, n
    real(8), intent(in) :: rxyz(3, n)
    real(8), intent(in) :: nxyz(3, n)
    real(8), intent(in) :: w(n)
    real(8), intent(out) :: a(3*n, 3*n)

    integer :: nu, nv, np
    integer :: i, j, k, c, src, t, ii, jj, indg
    real(8) :: pi, theta_i, d, ndotr
    real(8), allocatable :: xg(:), gwt(:), pols(:), wt_theta(:), ywt(:)
    real(8), allocatable :: utheta(:), wsph(:), w0(:)
    real(8), allocatable :: x(:), y(:), z(:), nx(:), ny(:), nz(:)
    real(8), allocatable :: rall(:,:,:), rmat(:,:)
    integer, allocatable :: ind(:)
    real(8), allocatable :: xx(:), yy(:), zz(:), wk(:), invrho(:)
    real(8), allocatable :: gx(:), gy(:), gz(:), g(:), lv(:), rowv(:)
    real(8), allocatable :: g11(:,:), g22(:,:), g33(:,:), g12(:,:), g13(:,:), g23(:,:)

    nu = p + 1
    nv = 2 * p
    np = nu * nv
    pi = 4.0d0 * atan(1.0d0)

    if (n /= np) stop 'kerneldSMatrix: n mismatch'

    allocate(xg(nu), gwt(nu), pols(p+1), wt_theta(nu), ywt(np))
    allocate(utheta(nu), wsph(np), w0(np))
    allocate(x(np), y(np), z(np), nx(np), ny(np), nz(np))
    allocate(rall(np,np,nu), rmat(np,np))
    allocate(ind(np))
    allocate(xx(np), yy(np), zz(np), wk(np), invrho(np))
    allocate(gx(np), gy(np), gz(np), g(np), lv(np), rowv(np))
    allocate(g11(np,np), g22(np,np), g33(np,np), g12(np,np), g13(np,np), g23(np,np))

    x = rxyz(1,:)
    y = rxyz(2,:)
    z = rxyz(3,:)
    nx = nxyz(1,:)
    ny = nxyz(2,:)
    nz = nxyz(3,:)

    call gl_nodes(nu, xg, gwt)
    xg = -xg

    do i = 1, nu
      call legepols(xg(i), p, pols)
      d = sum(pols)
      theta_i = acos(xg(i))
      wt_theta(i) = (pi / dble(p)) * gwt(i) * d / cos(theta_i / 2.0d0)
      utheta(i) = theta_i
    end do

    do k = 1, nv
      do i = 1, nu
        ywt((k-1)*nu + i) = wt_theta(i) / (8.0d0*pi)
        wsph((k-1)*nu + i) = sin(utheta(i))
      end do
    end do

    w0 = w / wsph

    do j = 1, nu
      call rot_mat(p, np, utheta(j), rmat)
      rall(:,:,j) = rmat
    end do

    g11 = 0.0d0; g22 = 0.0d0; g33 = 0.0d0
    g12 = 0.0d0; g13 = 0.0d0; g23 = 0.0d0

    do k = 1, nv
      t = 0
      do c = 1, nv
        src = mod(c - (k-1) - 1 + nv, nv) + 1
        do i = 1, nu
          t = t + 1
          ind(t) = i + (src-1)*nu
        end do
      end do

      do j = 1, nu
        indg = j + nu*(k-1)

        do ii = 1, np
          do jj = 1, np
            rmat(ii,jj) = rall(ind(ii), ind(jj), j)
          end do
        end do

        xx = matmul(rmat, x)
        yy = matmul(rmat, y)
        zz = matmul(rmat, z)
        wk = matmul(rmat, w0) * wsph

        do i = 1, np
          invrho(i) = 1.0d0 / sqrt((xx(i)-x(indg))**2 + (yy(i)-y(indg))**2 + (zz(i)-z(indg))**2)
          gx(i) = (xx(i)-x(indg)) * invrho(i)
          gy(i) = (yy(i)-y(indg)) * invrho(i)
          gz(i) = (zz(i)-z(indg)) * invrho(i)
          ndotr = nx(indg)*gx(i) + ny(indg)*gy(i) + nz(indg)*gz(i)
          g(i) = 6.0d0 * ywt(i) * wk(i) * ndotr * invrho(i) * invrho(i)
        end do

        lv = g * gx * gx
        rowv = matmul(lv, rmat)
        g11(indg,:) = rowv

        lv = g * gy * gy
        rowv = matmul(lv, rmat)
        g22(indg,:) = rowv

        lv = g * gz * gz
        rowv = matmul(lv, rmat)
        g33(indg,:) = rowv

        lv = g * gx * gy
        rowv = matmul(lv, rmat)
        g12(indg,:) = rowv

        lv = g * gx * gz
        rowv = matmul(lv, rmat)
        g13(indg,:) = rowv

        lv = g * gy * gz
        rowv = matmul(lv, rmat)
        g23(indg,:) = rowv
      end do
    end do

    do i = 1, np
      do j = 1, np
        a(3*i-2, 3*j-2) = g11(i,j)
        a(3*i-2, 3*j-1) = g12(i,j)
        a(3*i-2, 3*j  ) = g13(i,j)

        a(3*i-1, 3*j-2) = g12(i,j)
        a(3*i-1, 3*j-1) = g22(i,j)
        a(3*i-1, 3*j  ) = g23(i,j)

        a(3*i  , 3*j-2) = g13(i,j)
        a(3*i  , 3*j-1) = g23(i,j)
        a(3*i  , 3*j  ) = g33(i,j)
      end do
    end do

    deallocate(xg, gwt, pols, wt_theta, ywt)
    deallocate(utheta, wsph, w0, x, y, z, nx, ny, nz)
    deallocate(rall, rmat, ind)
    deallocate(xx, yy, zz, wk, invrho, gx, gy, gz, g, lv, rowv)
    deallocate(g11, g22, g33, g12, g13, g23)
  end subroutine kerneldSMatrix

end module stokes_mod
