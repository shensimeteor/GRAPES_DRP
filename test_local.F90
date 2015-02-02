#include "inc_common.fpp"
program test_local
use utility
use module_local
use module_types
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="main"
    integer,parameter :: len_modes=101, n_coarse=41
    real :: xcoord(len_modes), Cov1(len_modes, len_modes), xcoarse(n_coarse), Cov2(n_coarse, n_coarse)
    real :: Cov12(len_modes, n_coarse)
    type(type_localmodes) :: lm, lm2
    integer :: i
    real :: r0=10.0
    call fspan(0., 100., len_modes, xcoord)
    call fspan(0., 100., n_coarse, xcoarse)

    call localmodes_generate_decompC(lm, xcoord, r0, 1, .false., 100.,op_eig2sumratio=0.9999 )
    call localmodes_write("Local_Modes.txt", xcoord, lm)
    call localmodes_read("Local_Modes.txt", op_lm=lm)

    call localmodes_interp_linear(xcoord, lm, xcoarse, lm2, .false.)
    print*,lm%lenn, lm%m
    print*,lm2%lenn, lm2%m
    Cov2=matmul(lm2%modes, transpose(lm2%modes))
    call write_float_matrix2(Cov2, fname="Cov2.txt",element_fmat="F")
    Cov1=matmul(lm%modes, transpose(lm%modes))
    call write_float_matrix2(Cov1, fname="Cov1.txt",element_fmat="F")
    Cov12=matmul(lm%modes, transpose(lm2%modes))
    call write_float_matrix2(Cov12, fname="Cov12.txt",element_fmat="F")
    
    call write_float_matrix2(lm%modes, fname="lm1.txt")
    call write_float_matrix2(lm2%modes, fname="lm2.txt")

    
end program test_local
