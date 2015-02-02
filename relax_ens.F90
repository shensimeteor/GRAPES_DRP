! --------------------------------------------------------------------------
! Description:
!
!   relax ensemble: xa_i' <= coef_a * xa_i' + (1-coef_a) * xb_i
!   ./relax_ens.exe  <list_px>  <list_fgx>  <out_fgx_dir>
!
! History:
!
!   2014-09-15:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
program relax_ens
use utility
use ArgsManager
use module_types
use module_modelstate
use module_ensemble
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="relax_ens"
    type(type_modeldesc) :: md
    type(type_ensPx) :: ensPx, ensFGx
    character(len=100) :: list_px, list_fgx, outdir
    character(len=150), allocatable :: pxfiles(:), fgxfiles(:), newfiles(:)
    integer :: n, i, t, npx, nfgx
    real,allocatable :: xb(:), xa(:)
    real :: coef_a=0.85, coef_b=0.15

    call args_readargs()
    n=args_get_narg()
    if(n/=3) then
        RP_ERROR("3 arguments need: <list_px> <list_fgx> <outdir>")
    end if
    call args_getarg(1, list_px)
    call args_getarg(2, list_fgx)
    call args_getarg(3, outdir)

    call allocread_line(trim(list_px), npx, pxfiles)
    call allocread_line(trim(list_fgx), nfgx, fgxfiles)
    if(npx /= nfgx) then
        RP_ERROR("input Px/FGx file number conflict!")
    end if
    allocate(newfiles(npx))
    do i=1, npx
        newfiles(i)=trim(outdir)//'/'//trim(kickout_filepath(trim(fgxfiles(i))))
        call system("cp "//trim(fgxfiles(i))//" "//trim(newfiles(i)))
    end do

    call modeldesc_init(md)
    call modeldesc_getcoord(md, trim(pxfiles(1)))
    call Px_read(md, pxfiles, ensPx)
    call Px_read(md, fgxfiles, ensFGx)
    allocate(xb(md%n_total), xa(md%n_total))
    call Px_average(ensPx, xb)
    call Px_average(ensFGx, xa)
    call Px_full2pert(ensPx, xb)
    call Px_full2pert(ensFGx, xa)

    ensFGx%Px=coef_a*ensFGx%Px+coef_b*ensPx%Px

    call Px_pert2full(ensFGx, xa)

    call Px_write(md, newfiles, ensFGx)

end program relax_ens
