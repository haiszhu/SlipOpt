subroutine build_geometry_bvp_mex(p, np, shape_name_c, slen, r, nx, xu, xv, w)
  use geometry_mod, only: build_geometry_bvp_f90
  use iso_c_binding, only: c_char, c_null_char
  implicit none
  integer, intent(in) :: p, np, slen
  character(c_char), intent(in) :: shape_name_c(slen)
  real(8), intent(inout) :: r(3,np), nx(3,np), xu(3,np), xv(3,np), w(np)

  integer :: i, nlen, shape_id
  character(len=32) :: name

  name = ''
  nlen = 0
  do i = 1, slen
    if (shape_name_c(i) == c_null_char) exit
    if (nlen < len(name)) then
      nlen = nlen + 1
      name(nlen:nlen) = shape_name_c(i)
    end if
  end do

  select case (trim(name))
  case ('sphere')
    shape_id = 1
  case ('Y43')
    shape_id = 2
  case default
    stop 'build_geometry_bvp_mex: unsupported shape'
  end select

  call build_geometry_bvp_f90(p, np, shape_id, r, nx, xu, xv, w)
end subroutine build_geometry_bvp_mex
