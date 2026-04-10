module geometry_mod
  use sht_mod, only: gl_nodes
  implicit none
  private

  public :: build_geometry_bvp_f90

contains

  subroutine build_geometry_bvp_f90(p, np, shape_id, r, nx, xu, xv, w)
    integer, intent(in) :: p, np, shape_id
    real(8), intent(out) :: r(np, 3), nx(np, 3), xu(np, 3), xv(np, 3), w(np)

    integer :: nu, nv
    integer :: iu, iv, k
    integer, parameter :: n = 4, m = 3
    real(8), parameter :: amp = 0.5d0
    real(8) :: pi, twopi, theta, lambda
    real(8) :: rho, rho_t, rho_l
    real(8) :: erx, ery, erz, etx, ety, etz, elx, ely, elz
    real(8) :: xt1, xt2, xt3, xl1, xl2, xl3
    real(8) :: cx, cy, cz, jac
    real(8) :: pnm, yuline, cml, sml
    real(8), allocatable :: xg(:), gwt(:)

    nu = p + 1
    nv = 2 * p
    if (np /= nu*nv) stop 'build_geometry_bvp_f90: np mismatch'

    pi = 4.0d0 * atan(1.0d0)
    twopi = 2.0d0 * pi

    allocate(xg(nu), gwt(nu))
    call gl_nodes(nu, xg, gwt)
    xg = -xg

    k = 0
    do iv = 1, nv
      lambda = twopi * dble(iv-1) / dble(nv)

      do iu = 1, nu
        k = k + 1
        theta = acos(xg(iu))

        select case (shape_id)
        case (1) ! sphere
          rho = 1.0d0
          rho_t = 0.0d0
          rho_l = 0.0d0
        case (2) ! Y43
          call ynm_single_real_derivs(n, m, theta, lambda, pnm, yuline, rho_l)
          rho = 1.0d0 + amp * pnm
          rho_t = amp * yuline
          rho_l = amp * rho_l
        case default
          stop 'build_geometry_bvp_f90: unsupported shape_id'
        end select

        erx = sin(theta) * cos(lambda)
        ery = sin(theta) * sin(lambda)
        erz = cos(theta)

        etx = cos(theta) * cos(lambda)
        ety = cos(theta) * sin(lambda)
        etz = -sin(theta)

        elx = -sin(lambda)
        ely = cos(lambda)
        elz = 0.0d0

        r(k, 1) = rho * erx
        r(k, 2) = rho * ery
        r(k, 3) = rho * erz

        xt1 = rho_t * erx + rho * etx
        xt2 = rho_t * ery + rho * ety
        xt3 = rho_t * erz + rho * etz

        xl1 = rho_l * erx + rho * sin(theta) * elx
        xl2 = rho_l * ery + rho * sin(theta) * ely
        xl3 = rho_l * erz + rho * sin(theta) * elz

        xu(k, 1) = xt1
        xu(k, 2) = xt2
        xu(k, 3) = xt3

        xv(k, 1) = xl1
        xv(k, 2) = xl2
        xv(k, 3) = xl3

        cx = xt2*xl3 - xt3*xl2
        cy = xt3*xl1 - xt1*xl3
        cz = xt1*xl2 - xt2*xl1
        jac = sqrt(cx*cx + cy*cy + cz*cz)

        if (jac <= 0.0d0) stop 'build_geometry_bvp_f90: zero jacobian'

        nx(k, 1) = cx / jac
        nx(k, 2) = cy / jac
        nx(k, 3) = cz / jac
        w(k) = jac
      end do
    end do

    deallocate(xg, gwt)
  end subroutine build_geometry_bvp_f90


  subroutine ynm_single_real_derivs(n, m, theta, lambda, yreal, yureal, yvreal)
    integer, intent(in) :: n, m
    real(8), intent(in) :: theta, lambda
    real(8), intent(out) :: yreal, yureal, yvreal

    integer :: mm, l, q, mmin, mmax
    real(8) :: x, sx
    real(8) :: pmm, pmmp1, pll, pnm, pkm2, pkm1
    real(8) :: nc, pi, cml, sml
    real(8), allocatable :: pall(:), ppad(:), dpall(:), c(:)

    mmin = -n
    mmax = n + 1

    allocate(pall(2*n+1), ppad(2*n+3), dpall(2*n+1), c(2*n+2))
    pall = 0.0d0
    ppad = 0.0d0
    dpall = 0.0d0

    pi = 4.0d0 * atan(1.0d0)
    x = cos(theta)
    sx = sqrt(max(0.0d0, 1.0d0 - x*x))

    pmm = 1.0d0
    do mm = 0, n
      if (mm > 0) pmm = dble(2*mm-1) * sx * pmm

      if (n == mm) then
        pnm = pmm
      else
        pmmp1 = dble(2*mm+1) * x * pmm
        if (n == mm+1) then
          pnm = pmmp1
        else
          pkm2 = pmm
          pkm1 = pmmp1
          pll = pkm1
          do l = mm+2, n
            pll = (dble(2*l-1) * x * pkm1 - dble(l+mm-1) * pkm2) / dble(l-mm)
            pkm2 = pkm1
            pkm1 = pll
          end do
          pnm = pll
        end if
      end if

      nc = norm_const(n, mm, pi)
      pall(n+1+mm) = dble((-1)**mm) * nc * pnm
      if (mm > 0) then
        pall(n+1-mm) = dble((-1)**mm) * pall(n+1+mm)
      end if
    end do

    ppad(2:2*n+2) = pall(:)

    do q = 1, 2*n+2
      mm = mmin + (q-1)
      c(q) = sqrt((dble(n+mm)) * (dble(n+1-mm))) / 2.0d0
    end do

    do q = 1, 2*n+1
      dpall(q) = c(q) * ppad(q) - c(q+1) * ppad(q+2)
    end do

    cml = cos(dble(m) * lambda)
    sml = sin(dble(m) * lambda)

    yreal = pall(n+1+m) * cml
    yureal = -dpall(n+1+m) * cml
    yvreal = -dble(m) * pall(n+1+m) * sml

    deallocate(pall, ppad, dpall, c)
  end subroutine ynm_single_real_derivs


  real(8) function norm_const(n, mabs, pi)
    integer, intent(in) :: n, mabs
    real(8), intent(in) :: pi
    norm_const = sqrt((dble(2*n+1) / (4.0d0*pi)) * exp(log_gamma_ratio(n-mabs, n+mabs)))
  end function norm_const


  real(8) function log_gamma_ratio(a, b)
    integer, intent(in) :: a, b
    log_gamma_ratio = log_gamma(dble(a+1)) - log_gamma(dble(b+1))
  end function log_gamma_ratio

end module geometry_mod
