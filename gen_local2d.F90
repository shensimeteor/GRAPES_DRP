!!for global
#include "inc_common.fpp"
program gen_local2d
use ArgsManager
use utility
use module_local
use module_local2dman
use module_types
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="gen_local2d"
    integer,parameter :: nx=360, ny=181
    real :: dx, dy
    real :: xcoord(nx), ycoord(ny)
    type(type_local2dman) :: l2dm
    integer :: i
    integer :: n_band, n_xmode, n_ymode
    integer, allocatable :: yidx_bgn(:), yidx_end(:)
    real,allocatable :: xl0(:), dy_ratio(:)
    character(len=100),allocatable :: strline(:)
    character(len=150) :: dir_out, xmode_descfile
    real :: yl0  !! in degree

    call args_readargs()
    if(args_get_narg() /= 5) then
        RP_ERROR("need 5 argument: <descript file for xmodes_x0> <yl0: real in degree> <n_xmode> <n_ymode> <outdir>")
    end if
    call args_getarg(1, xmode_descfile)
    call args_getarg(2, yl0)
    call args_getarg(3, n_xmode)
    call args_getarg(4, n_ymode)
    call args_getarg(5, dir_out)

    dx=360./nx
    dy=180./(ny-1)

    call fspan(0., 360.-dx, nx, xcoord)
    call fspan(-90., 90.,  ny, ycoord)

    call local2dman_set_domain(l2dm, xcoord, ycoord)
    call allocread_line(trim(xmode_descfile), n_band, strline)
    print*, "n_band=", n_band
    allocate(xl0(n_band), yidx_bgn(n_band), yidx_end(n_band))
    do i=1, n_band
        read(strline(i), *) yidx_bgn(i), yidx_end(i), xl0(i)
    end do
    call local2dman_set_xmodes(l2dm, yidx_bgn, yidx_end, xl0, n_xmode, .true.)
    TRACEX(l2dm%n_xlm)
    allocate(dy_ratio(ny-1))
    dy_ratio=1.
    call local2dman_set_ymode(l2dm, dy_ratio, yl0, op_nmodes=n_ymode)

    
    call local2dman_rownormalize_xy(l2dm)
    call local2dman_output(l2dm, trim(dir_out))

    !!write auxiliary information
    open(11, file=trim(dir_out)//"/info.txt", status='replace')
    write(11,*) "Xmode info"
    write(11,*) "Decompose method: Fourier"
    write(11,*) "n_xmode=", n_xmode
    write(11,*) "band_xl0_info:"
    do i=1, n_band
        write(11, *) yidx_bgn(i), yidx_end(i), xl0(i)
    end do
    write(11,*) "n_ymode=", n_ymode
    write(11,*) "yl0=", yl0
    write(11,*) "max(dy_ratio)=", maxval(dy_ratio), "min(dy_ratio)=", minval(dy_ratio)
    write(11,*) "dy_ratio:"
    write(11,*) dy_ratio
    close(11)
end program 
