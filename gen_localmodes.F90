!!for global
#include "inc_common.fpp"
program gen_localmodes
use ArgsManager
use utility
use module_local
use module_types
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="gen_localmodes"
    integer,parameter :: nx=360, ny=181
    real :: dx, dy
    real :: xcoord(nx), Covx(nx, nx), ycoord(ny), Covy(ny, ny)
    type(type_localmodes) :: xlm, ylm
    integer :: i, narg
    real :: r0  !! in degree
    real :: eig_ratio
    
    call args_readargs()
    narg=args_get_narg()
    if( narg == 1) then
        call args_getarg(1, r0)
        ASSUREX((r0<100. .and. r0>0.))
        eig_ratio=0.98
    else if( narg == 2) then
        call args_getarg(1, r0)
        call args_getarg(2, eig_ratio)
        ASSUREX(r0<100. .and. r0>0.)
        ASSUREX((eig_ratio <= 1.0 .and. eig_ratio >= 0.))
    else
        RP_ERROR("<usage> ./gen_localmodes.exe <r0 length (in degree)> [<eig-decomp ratio]")
    end if

    dx=360./nx
    dy=180./(ny-1)

    call fspan(0., 360.-dx, nx, xcoord)
    call fspan(-90., 90.,  ny, ycoord)

    call localmodes_generate_decompC(xlm, xcoord, r0, 1, .true., 360.,op_eig2sumratio=eig_ratio )
    call localmodes_write_csv("global_xmodes.csv", xcoord, xlm)

    call localmodes_generate_decompC(ylm, ycoord, r0, 1, .false., 0., op_eig2sumratio=eig_ratio )
    call localmodes_write_csv("global_ymodes.csv", ycoord, ylm)
    
 !! test read/write_csv, passed
!    call localmodes_delete(ylm)
!    ycoord=0
!    call localmodes_read_csv("global_ymodes.csv", ycoord, ylm)
!    call localmodes_write_csv("test_ymodes.csv", ycoord, ylm)
    
end program 
