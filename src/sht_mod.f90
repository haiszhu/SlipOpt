module sht_mod
  implicit none
  private

  public :: gl_nodes
  public :: shc_expand
  public :: shc_shrink

  integer, save :: p_cached = -1
  integer, save, allocatable :: idx_compact_cache(:)

  interface
    subroutine legeexps(itype, n, x, u, v, whts)
      integer, intent(in) :: itype, n
      real(8), intent(out) :: x(n), u(n,n), v(n,n), whts(n)
    end subroutine legeexps
  end interface

contains

  subroutine gl_nodes(n, x, w)
    integer, intent(in) :: n
    real(8), intent(out) :: x(n), w(n)

    real(8), allocatable :: u(:,:), v(:,:), xr(:)
    integer :: itype

    if (n <= 0) return

    allocate(u(n,n), v(n,n), xr(n))
    itype = 1
    call legeexps(itype, n, xr, u, v, w)
    x = -xr
    deallocate(u, v, xr)
  end subroutine gl_nodes


  subroutine shc_expand(p, nc, shc, shc_full)
    integer, intent(in) :: p, nc
    complex(8), intent(in) :: shc((p+1)*(p+1), nc)
    complex(8), intent(out) :: shc_full((p+1)*(2*p+1), nc)

    integer :: k, ncoef

    call ensure_cache(p)
    ncoef = (p+1)*(p+1)

    shc_full = (0.0d0, 0.0d0)
    do k = 1, ncoef
      shc_full(idx_compact_cache(k), :) = shc(k, :)
    end do
  end subroutine shc_expand


  subroutine shc_shrink(p, nc, shc_full, shc)
    integer, intent(in) :: p, nc
    complex(8), intent(in) :: shc_full((p+1)*(2*p+1), nc)
    complex(8), intent(out) :: shc((p+1)*(p+1), nc)

    integer :: k, ncoef

    call ensure_cache(p)
    ncoef = (p+1)*(p+1)

    do k = 1, ncoef
      shc(k, :) = shc_full(idx_compact_cache(k), :)
    end do
  end subroutine shc_shrink


  subroutine ensure_cache(p)
    integer, intent(in) :: p

    integer :: nu, ncoef
    integer :: n, m, k

    if (p == p_cached) return

    if (allocated(idx_compact_cache)) deallocate(idx_compact_cache)

    nu = p + 1
    ncoef = nu * nu
    allocate(idx_compact_cache(ncoef))

    k = 0
    do n = 0, p
      do m = -n, n
        k = k + 1
        idx_compact_cache(k) = (m + p) * nu + (n + 1)
      end do
    end do

    p_cached = p
  end subroutine ensure_cache


end module sht_mod
