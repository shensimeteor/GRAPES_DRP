!!for global
#include "inc_common.fpp"
program test_l2dm
use utility
use module_local
use module_local2dman
use module_types
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="test_l2dm"
    integer,parameter :: nx=360, ny=181, nx2=80, ny2=41
    real :: dx, dy, dx2, dy2
    real :: xcoord(nx), ycoord(ny), x2(nx2), y2(ny2)
    type(type_local2dman) :: l2dm, l2dm2
    integer :: i
    integer :: n_band, n_xmode=50, n_ymode=50
    integer, allocatable :: yidx_bgn(:), yidx_end(:)
    real,allocatable :: xl0(:), dy_ratio(:)
    character(len=100),allocatable :: strline(:)
    real :: yl0=7.5  !! in degree

    TIME_START(timer_bgn)
    dx=360./nx
    dy=180./(ny-1)

    call fspan(0., 360.-dx, nx, xcoord)
    call fspan(-90., 90.,  ny, ycoord)

    call local2dman_set_domain(l2dm, xcoord, ycoord)
    call allocread_line("xmodes_xl0.txt", n_band, strline)
    print*, "n_band=", n_band
    allocate(xl0(n_band), yidx_bgn(n_band), yidx_end(n_band))
    do i=1, n_band
        read(strline(i), *) yidx_bgn(i), yidx_end(i), xl0(i)
    end do
    call local2dman_set_xmodes(l2dm, yidx_bgn, yidx_end, xl0, n_xmode, .true.)
    allocate(dy_ratio(ny-1))
    dy_ratio=1.
    call local2dman_set_ymode(l2dm, dy_ratio, yl0, op_nmodes=n_ymode)
    print*, "after set x&y"
    TIME_CLICK(timer, timer_bgn, timer_last)
    call local2dman_rownormalize_xy(l2dm)
    print*, "after rownormalize"
    TIME_CLICK(timer, timer_bgn, timer_last)
    call local2dman_output(l2dm, "2dmodes")
    print*, "after output"
    TIME_CLICK(timer, timer_bgn, timer_last)

    !!
    call local2dman_dealloc(l2dm)
    call local2dman_input(l2dm, "2dmodes")
    print*, "after input"
    TIME_CLICK(timer, timer_bgn, timer_last)

    dx2=360./nx2
    dy2=180./ny2
    call fspan(0., 360.-dx2, nx2, x2)
    call fspan(-90., 90.,  ny2, y2)
    call local2dman_interp_togrid(l2dm, l2dm2, x2, y2, .true.)
    print*, "after interp_togrid"
    TIME_CLICK(timer, timer_bgn, timer_last)

    call local2dman_output(l2dm2, "2dmodes_2")
    print*, "after re-output"
    TIME_CLICK(timer, timer_bgn, timer_last)

    
end program 
