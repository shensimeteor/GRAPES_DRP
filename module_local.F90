! --------------------------------------------------------------------------
! Description:
!
!   localization, concentrate on modes (model independent)
!   localmodes_generate_* (different methods to generate modes: sin, decompose C) or read: in localmodes_read
!   localmodes_read_csv (read modes from file)
!   localmodes_interp_* (interp lm1 to lm2; linear interp or 3-spline interp)
!   localmodes_rownormalize (normalize every row of modes)
!   localmodes_write_csv (write modes to file)
!   localmodes_delete (delete lm)
!   localmodes_init ( most times it need not be called, but for some special cases )
!
! History:
!
!   2014-07-02:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
module module_local
use utility
use sort
use interp
use lapack95
use f95_precision
use module_types
implicit none
private
    ! variables
    integer,parameter :: CMETHOD_GAUSS=1
    public :: localmodes_generate_nolocal, localmodes_generate_decompC, &
              localmodes_interp_linear, localmodes_read_csv, localmodes_write_csv, &
              localmodes_rownormalize, localmodes_delete, localmodes_init, &
              localmodes_read, localmodes_write, localmodes_interp_linear_nmg, &
              localmodes_generate_Fourier_pole, localmodes_generate_Fourier,&
              localmodes_linear_combine !! this 2 is depressed
    real,parameter :: PI=3.1415926
contains

!!generate
    subroutine localmodes_generate_nolocal(lm, lenn)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_generate_nolocal"
        type(type_localmodes) :: lm
        integer :: lenn
        call localmodes_init(lm, lenn, 1)
        lm%modes=1
    end subroutine localmodes_generate_nolocal  !!}}}

    !!if <op_nmodes> is present, select top <op_nmodes> modes
    !!if <op_eig2sumratio> is present, select top modes contribute more than <op_eig2sumratio>
    !!if both present, is top X modes > op_eig2sumratio, then select X; else, select <op_nmodes> (X is upper limit)
    !!if no present, select all modes
    subroutine localmodes_generate_decompC(lm, x, r0, Cmethod, iscycle, cycle_period, op_nmodes, op_eig2sumratio)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_generate_decompC"
        type(type_localmodes) :: lm
        real :: x(:),r0
        integer :: Cmethod
        logical :: iscycle
        real :: cycle_period
        real :: dis, dis2
        integer :: nx,i,j, stat, final_nmodes
        real :: C(size(x),size(x))
        integer :: sort_idx(size(x))
        integer,optional :: op_nmodes
        real,optional :: op_eig2sumratio 
        real :: eigvalue(size(x)), sum_eig2, ratio, sum_head
        if(Cmethod/=CMETHOD_GAUSS) then
            RP_ERROR("Cmethod not supported")
        end if
        !!calculate C
        nx=size(x)
        C=1
        do i=1,nx
            do j=1,i-1
                dis=x(i)-x(j)
                if(iscycle) then
                    dis2=x(j)+cycle_period-x(i)
                    if(dis2<dis) dis=dis2
                end if
                C(i,j)=getCorr_Gauss(dis/r0)
                C(j,i)=C(i,j)
            end do
        end do
        !!decompose C & sort
        call syevd(C, eigvalue, jobz='V',uplo='L',info=stat)
        LPK_CHECK(stat)
        where(eigvalue<0)
            eigvalue=0
        end where
        call sorti(eigvalue, sort_idx, "descend")
        !!select modes
        sum_eig2=sum(eigvalue **2)
        if(present(op_nmodes)) then
            if(op_nmodes > nx ) then
                RP_ERROR("op_nmodes > nx!")
            end if
        end if
        if(present(op_eig2sumratio)) then
            if(op_eig2sumratio > 1) then
                RP_ERROR("op_eig2sumratio > 1!")
            end if
        end if
        if(present(op_nmodes) .and. .not. present(op_eig2sumratio)) then
            final_nmodes=op_nmodes
        else if(.not. present(op_nmodes) .and. present(op_eig2sumratio)) then
            sum_head=0
            do i=1, nx
                sum_head=sum_head+eigvalue(sort_idx(i))**2
                ratio=sum_head/sum_eig2
                TRACE4("top ",i," eigvalues contribute:",ratio)
                if(ratio >= op_eig2sumratio) exit
            end do
            final_nmodes=i
        else if(present(op_nmodes) .and. present(op_eig2sumratio)) then
            sum_head=0
            do i=1, op_nmodes
                sum_head=sum_head+eigvalue(sort_idx(i))**2
                ratio=sum_head/sum_eig2
                TRACE4("top ",i," eigvalues contribute:",ratio)
                if(ratio >= op_eig2sumratio) exit
            end do
            final_nmodes=i
        else !! both no present
            final_nmodes=nx
        end if
        ratio=sum(eigvalue(sort_idx(1:final_nmodes))**2 / sum_eig2)
        RP_LOG("select "//trim(tostring(ival=final_nmodes))//" modes, contribute "//trim(tostring(fval=ratio)))
        !!assign modes
        call localmodes_init(lm, nx, final_nmodes)
        do i=1, final_nmodes
            lm%modes(:,i)=C(:, sort_idx(i))*sqrt(eigvalue(sort_idx(i)))
        end do
    end subroutine localmodes_generate_decompC  !!}}}

    !! for only cycle coordinate, x must be uniformed-distributed, cycle_period=nx*dx
    !! modes={1, sin(2pi/L*x), cos(2pi/L*x), sin(4pi/L*x), sin(4pi/L*x), ... } then normalized
    subroutine localmodes_generate_Fourier(lm, nx, dx, r0, Cmethod, n_mode)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_generate_Fourier"
        type(type_localmodes) :: lm
        real :: dx, r0, cycle_period, dis, dis2
        integer :: Cmethod, n_mode, nx, i,j,t
        logical :: is_cycle=.true.
        real :: C(nx,nx), xcoord(nx), len2(n_mode), beta(n_mode), acc
        if(Cmethod/=CMETHOD_GAUSS) then
            RP_ERROR("Cmethod not supported")
        end if
        if(is_cycle) then
            cycle_period=nx*dx
        else
            cycle_period=(nx-1)*dx
        end if
        C=0
        do i=1,nx
            C(i,i)=1
            do j=1,i-1
                dis=(i-j)*dx
                if(is_cycle) then
                    dis2=(j+nx-i)*dx
                    if(dis2<dis) dis=dis2
                end if
                C(i,j)=getCorr_Gauss(dis/r0)
                C(j,i)=C(i,j)
            end do
        end do
       ! print*, r0
       ! print*, C(1:5,1:5)
        call localmodes_init(lm, nx,  n_mode)
        call fspan(0., (nx-1)*dx, nx, xcoord)
        do j=1,n_mode
            t=mod(j,2)
            i=j/2
            if(t==1) then
                lm%modes(:,j)=cos(i*2.*PI/cycle_period*xcoord)
            else
                lm%modes(:,j)=sin(i*2.*PI/cycle_period*xcoord)
            end if
        end do
        do j=1,n_mode
            len2(j)=sum(lm%modes(:,j)**2)
        end do
        do j=1,n_mode
            acc=0
            do i=1,nx
                do t=1,nx
                    acc=acc+C(i,t)*lm%modes(i,j)*lm%modes(t,j)
                end do
            end do
            beta(j)=acc/(len2(j)**2)
        end do
        do j=1,n_mode
            if(beta(j) <0) then
                RP_WARNING("for mode="//trim(tostring(ival=j))//"  beta(j)="// &
                           trim(tostring(fval=beta(j)))//"  WE set it 0")
                beta(j)=0
            end if
            lm%modes(:,j)=lm%modes(:,j)*sqrt(beta(j))
        end do
    end subroutine localmodes_generate_Fourier  !!}}}

    !!at pole, C(:,:)=1 , so localmodes is simple
    subroutine localmodes_generate_Fourier_pole(lm, nx, n_mode)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_generate_Fourier_pole"
        type(type_localmodes) :: lm
        integer :: n_mode, nx
        call localmodes_init(lm, nx, n_mode)
        lm%modes=0
        lm%modes(:,1)=1
    end subroutine localmodes_generate_Fourier_pole  !!}}}

!!DEPRESSED
 !   subroutine localmodes_generate_decompC_2dglobe(lm1, lm2, lon, lat, r0_km, Cmethod, op_nmodes, op_eig2sumratio)  !!{{{
 !   implicit none
 !       character(len=*),parameter :: PROCEDURE_NAME="localmodes_generate_decompC_2dglobe"
 !   end subroutine localmodes_generate_decompC_2dglobe  !!}}}

    !!DEPRESSED
    !!for region (not global) domain
!    subroutine localmodes_generate_sin(lm, x, r0, nmode, expand_deta_x )  !!{{{
!    implicit none
!        character(len=*),parameter :: PROCEDURE_NAME="localmodes_generate_sin"
!    end subroutine localmodes_generate_sin  !!}}}


!!delete
    subroutine localmodes_delete(lm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_delete"
        type(type_localmodes):: lm
        if(allocated(lm%modes)) then
            deallocate(lm%modes)
            lm%m=0
            lm%lenn=0
        end if
    end subroutine localmodes_delete  !!}}}


!! do interp 
    !!(assume x1 is regular(sx,nx,dx,iscycle determined))
    subroutine localmodes_interp_linear(x1, lm1, x2, lm2, iscycle)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_interp_linear"
        real :: x1(:), x2(:)
        type(type_localmodes) :: lm1, lm2
        integer :: nx1, nx2,i
        logical :: iscycle
        nx1=size(x1)
        nx2=size(x2)
        call localmodes_init(lm2, nx2, lm1%m)
        do i=1, lm1%m
            call interp_1d_linear_reg2sct(x1(1), nx1, x1(2)-x1(1), iscycle, lm1%modes(:,i), x2, lm2%modes(:,i))
        end do
    end subroutine localmodes_interp_linear  !!}}}

    !! for NorMal but not regular/cyclic grid (vertial )
    subroutine localmodes_interp_linear_nmg(x1, lm1, x2, lm2)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_interp_linear_nmg"
        real :: x1(:), x2(:)
        type(type_localmodes) :: lm1, lm2
        integer :: nx1, nx2,i
        nx1=size(x1)
        nx2=size(x2)
        call localmodes_init(lm2, nx2, lm1%m)
        do i=1, lm1%m
            call interp_1d_linear_nmg2sct(x1, lm1%modes(:,i), x2, lm2%modes(:,i))
        end do
    end subroutine localmodes_interp_linear_nmg  !!}}}


    !! lm=a1*lm1+a2*lm2
    subroutine localmodes_linear_combine(lm1, a1, lm2, a2, lm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_linear_combine"
        type(type_localmodes) :: lm1, lm2, lm
        real :: a1, a2
        integer :: i
        if( .not. all((/lm1%lenn, lm1%m/) == (/lm2%lenn, lm2%m/)) ) then
            RP_ERROR("input para. conflict")
        end if
        if(lm%m/=0) then
            call localmodes_delete(lm)
        end if
        call localmodes_init(lm, lm1%lenn, lm1%m)
        lm%modes=a1*lm1%modes+a2*lm2%modes
    end subroutine localmodes_linear_combine  !!}}}

!! rownormalize
    subroutine localmodes_rownormalize(lm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_rownormalize"
        type(type_localmodes) :: lm
        integer :: i
        real,allocatable ::  row_len(:) 
        allocate(row_len(lm%lenn))
        do i=1, lm%lenn
            row_len(i)=sum(lm%modes(i,:)**2)
        end do
        do i=1, lm%lenn
            lm%modes(i,:)=lm%modes(i,:)/sqrt(row_len(i))
        end do
    end subroutine localmodes_rownormalize  !!}}}


!! write
    !!DEPRESSED
    subroutine localmodes_write(lmfile, x, lm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_write"
        character(len=*) :: lmfile
        real :: x(:)
        type(type_localmodes) :: lm
        integer :: nx,i
        nx=size(x)
        if(nx /= lm%lenn) then
            RP_ERROR("nx/=lm%lenn")
        end if
        open(31, file=trim(lmfile))
        write(31,*) nx
        write(31,*) x
        write(31,*) lm%m
        do i=1,lm%m
            write(31,*) i
            write(31,*) lm%modes(:,i)
        end do
        close(31)
    end subroutine localmodes_write  !!}}}

    subroutine localmodes_write_csv(lmfile, x, lm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_write"
        character(len=*) :: lmfile
        real :: x(:)
        type(type_localmodes) :: lm
        integer :: nx,i,j
        character(len=40) :: strfmt
        nx=size(x)
        if(nx /= lm%lenn) then
            RP_ERROR("nx/=lm%lenn")
        end if
        open(31, file=trim(lmfile), status='replace')
        write(31, '(I,",",I)') nx, lm%m
        write(strfmt, * ) "(A, ", lm%m, "(A, I2.2, A))"
        write(31, trim(strfmt))  "coord, ", (( "mode",i,","),i=1,lm%m)
        write(strfmt, * ) "(F,A, ", lm%m, "(F, A))"
        do i=1,nx
            write(31, trim(strfmt)) x(i), ",", ((lm%modes(i,j),","), j=1,lm%m)
        end do
        close(31)
    end subroutine localmodes_write_csv  !!}}}


!! read
    !!DEPRESSED
    subroutine localmodes_read(lmfile, op_x, op_lm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_read"
        character(len=*) :: lmfile
        real, optional :: op_x(:)
        type(type_localmodes),optional :: op_lm
        real, allocatable :: temp(:)
        integer :: stat, i, nx, idx
        open(31,file=trim(lmfile),iostat=stat,status='old')
        if(stat /= 0) then
            RP_ERROR("open file "//trim(lmfile)//" error")
        end if
        read(31,*) nx
        if(present(op_x)) then
            if(size(op_x) < nx) then
                RP_ERROR("size(op_x) < nx")
            end if
            read(31, *) op_x(1:nx)
        else
            allocate(temp(nx))
            read(31, *) temp
            deallocate(temp)
        end if
        if(present(op_lm)) then
            call localmodes_delete(op_lm)
            read(31, *) idx
            call localmodes_init(op_lm, nx, idx) 
            do i=1, op_lm%m
                read(31,*) idx
                if(idx /= i) then
                    RP_ERROR("read localmodes "//trim(lmfile)//" error: expect "//tostring(ival=i))
                end if
                read(31,*) op_lm%modes(:,i)
            end do
        end if
        close(31)
    end subroutine localmodes_read  !!}}}

    !!  x should be allocatable , if x not allocated, allocate x
    subroutine localmodes_read_csv(lmfile, x, lm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_read_csv"
        character(len=*) :: lmfile
        real,allocatable :: x(:)
        type(type_localmodes) :: lm
        integer :: stat, i, j, nx
        open(31, file=trim(lmfile), status='old', iostat=stat)
        if(stat /= 0) then
            RP_ERROR("open file "//trim(lmfile)//" error")
        end if
        read(31, *) nx, i
        if(.not. allocated(x)) then
            allocate(x(nx))
        else if(size(x) < nx) then
            RP_ERROR("size(x) < nx")
        end if
        call localmodes_delete(lm)
        call localmodes_init(lm, nx, i)
        read(31, *) 
        do i=1,nx
            read(31, *) x(i), (lm%modes(i,j), j=1,lm%m)
        end do
        close(31)
    end subroutine localmodes_read_csv  !!}}}


    subroutine localmodes_init(lm, lenn, m)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localmodes_init"
        type(type_localmodes) :: lm
        integer :: lenn, m
        lm%lenn=lenn
        lm%m=m
        allocate(lm%modes(lenn, m))
    end subroutine localmodes_init  !!}}}

    function getCorr_Gauss(dis)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="getCorr_Gauss"
        real :: dis, getCorr_Gauss
        getCorr_Gauss=exp(-dis**2/2)
    end function getCorr_Gauss  !!}}}

end module module_local
