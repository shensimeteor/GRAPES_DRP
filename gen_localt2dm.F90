!!for global
#include "inc_common.fpp"
program gen_localt2dm
use ArgsManager
use utility
use module_local
use module_localtrue2dman
use module_types
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="gen_localt2dm"
    integer,parameter :: nx=90, ny=46
    real :: dx, dy
    real :: xcoord(nx), ycoord(ny)
    type(type_localtrue2dman) :: lt2dm
    integer :: i
    real :: L0
    character(len=100),allocatable :: strline(:)
    character(len=150) :: output_prefix
    real :: yl0  !! in degree
    real :: t, tbgn, tlast=-1

    call args_readargs()
    if(args_get_narg() /= 2) then
        RP_ERROR("need 2 argument: L0 output_prefix")
    end if
    call args_getarg(1, L0)
    call args_getarg(2, output_prefix)
    dx=360./nx
    dy=180./(ny-1)

    call fspan(0., 360.-dx, nx, xcoord)
    call fspan(-90., 90.,  ny, ycoord)
    
    TIME_START(tbgn)
    call localt2dm_set_domain(lt2dm, xcoord, ycoord)
    
    call localt2dm_generate_decompC2d(lt2dm, L0, op_nmodes=500, op_eig2sumratio=0.98)
    TIME_CLICK(t, tbgn, tlast)

    call localt2dm_write_binary(lt2dm, trim(output_prefix))
    TIME_CLICK(t, tbgn, tlast)

end program 
