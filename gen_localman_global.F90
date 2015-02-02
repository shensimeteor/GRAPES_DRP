!!for global, gen xm/ym/zm
#include "inc_common.fpp"
program gen_localman_global
use ArgsManager
use utility
use module_local
use module_localman
use module_types
use lapack95
use f95_precision
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="gen_localman_global"
    integer,parameter :: nx=180, ny=91, nz=3 !!Localization Grid
    type(type_localman) :: lman
    real :: dx, dy, lat
    real :: xcoord(nx), ycoord(ny), zcoord(nz), Cz(nz,nz), eigs(nz), Czz(nz,nz)
    real :: zratio(nz)=(/1.0, 1.2, 1.5/) !! only for mz=3
    real :: L0km 
    real, parameter :: R0=6400, PI=atan(1.)*4, D2R=PI/180., R2D=180./PI
    integer :: my, mz, px, py, j, k, n_zmode
    integer, allocatable :: ybgn(:), yend(:), zbgn(:), zend(:)
    integer :: flag1, flag2
    real,allocatable :: xl0(:,:), yl0(:), dy_ratio(:)
    character(len=100),allocatable :: strline(:)
    character(len=150) :: dir_out, xmode_descfile
    character(len=200) :: description
    
    print*,PI

    call args_readargs()
    if(args_get_narg() /= 5) then
        print*, "need 5 argument: <L0_km> <my> <mz> <px> <py>"
        print*, "my/mz: control the sampling interval for y / z"
        print*, "px/py: the number of modes for xm/ym"
        stop 2
    end if
    call args_getarg(1, L0km)
    call args_getarg(2, my)
    call args_getarg(3, mz)
    call args_getarg(4, px)
    call args_getarg(5, py)
    if(mz /= 1 ) then
        print*," mz must ==1 NOW"
        stop 2
    end if

    dx=360./nx
    dy=180./(ny-1)
    call fspan(0., 360.-dx, nx, xcoord)
    call fspan(-90., 90.,  ny, ycoord)
    zcoord=(/800, 500, 200/)

    flag1=mod(ny-1, my-1) !!90: 46, 31, 19, 16, 11 .. are OK
    flag2=mod(my-1, 2) !! 90: 31, 19, 11 ..
    if(flag1 /= 0 .or. flag2 /= 0) then
        RP_ERROR(" for <my>: ny-1 should be divided by my; my should be even : e.g. my=")
    end if
    if(mz > nz) then
        RP_WARNING(" mz should no more than nz, mz is as nz ")
        mz=nz
    end if

    call localman_set_domain(lman, xcoord, ycoord, zcoord, .true.)
    allocate(ybgn(my), yend(my), zbgn(mz), zend(mz))
    allocate(xl0(my,mz), yl0(mz))
    !! todo: calc. xl0, yl0 (for now nz=1 case), set_xm/set_xm, output, test & plot output
    do j=1, my
        lat= -90 + (ny-1)/(my-1) * (j-1) * dy 
        if ( j == 1 .or. j == my ) then
            xl0(j,1) = -1
        else
            xl0(j,1) = L0km/(R0*cos(lat*D2R)) * R2D
        end if
        ybgn(j)=1+(ny-1)/(my-1) * (j-1)
        yend(j)=ybgn(j)
    end do 
    yl0=L0km/(R0) * R2D
    if(mz==3) then
        zbgn(1)=1
        zend(1)=1
        zbgn(2)=2
        zend(2)=2
        zbgn(3)=3
        zend(3)=3
        do j=1,my
            xl0(j,2)=xl0(j,1)*zratio(2)
            xl0(j,3)=xl0(j,1)*zratio(3)
            xl0(j,1)=xl0(j,1)*zratio(1)
        end do
        yl0(2)=yl0(1)*zratio(2)
        yl0(3)=yl0(1)*zratio(3)
        yl0(1)=yl0(1)*zratio(1)
    else !! mz==1
        zbgn(1)=1
        zend(1)=3
    end if
    print*, px
    print*, ybgn, yend
    print*, xl0
    print*, yl0
    print*, zbgn, zend
    print*, "TO SET XM"
    call localman_set_xm(lman, ybgn, yend, zbgn, zend, xl0, px)
    print*, "TO SET YM"
    call localman_set_ym_byEIG(lman, yl0(1), py)
    !! for zm
    print*, "TO SET ZM"
    call set_Cz_LiuCS(zcoord, Cz)
    Czz=Cz
    call decompose_C_syevd(Cz, eigs, 0.99, n_zmode)
    call localmodes_init(lman%zm, nz, n_zmode)
    do k=1, lman%zm%m
        lman%zm%modes(:,k)=Cz(:,k)*sqrt(eigs(k))
    end do
    lman%has_set_zm=.true.
!    call localman_set_zm_nolocal(lman)
    call localman_rownormalize(lman)
    call check_Cxyz_outfile(lman, xl0, yl0, Czz)
    write(description, *) "L0km=",L0km, "   ; my=",my,"  ;mz=",mz
    call localman_output_xm(lman, "global_xm.nc", trim(description))
    call localman_output_ym(lman, "global_ym.nc", trim(description))
    call localman_output_zm(lman, "global_zm.nc", trim(description))


contains
!! 1/(1+5(log(P2/P1))**2)
subroutine set_Cz_LiuCS(zcoord, Cz)  !!{{{
use utility
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="set_Cz_LiuCS"
    real :: zcoord(:), Cz(:,:), dis
    integer :: i,j, n
    n=size(zcoord)
    do i=1, n
        Cz(i,i)=1
        do j=1,i-1
           dis=log(zcoord(i)/zcoord(j))
           Cz(i,j)=1./(1+5*dis**2)
           Cz(j,i)=Cz(i,j)
        end do
    end do
end subroutine set_Cz_LiuCS  !!}}}

!!eigvalue & eigvector are order of size(C,1)
!!give :: eig2sumratio return nmodes
!!C: output the eigvector matrix
subroutine decompose_C_syevd(C, eigvalue, eig2sumratio, nmodes)  !!{{{
use utility
use sort
use lapack95
use f95_precision
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="decompose_C_syevd"
    real :: eig2sumratio
    integer :: nmodes
    real :: C(:,:), eigvalue(:)
    real,allocatable :: eval(:)
    integer :: n, stat, final_nmodes, i
    integer,allocatable :: sort_idx(:)
    real :: sum_eig2, sum_head, ratio
    n=size(C,1)
    allocate(eval(n),sort_idx(n))
    call syevd(C, eval, jobz='V',uplo='L',info=stat)
    LPK_CHECK(stat)
    where(eval<0)
        eval=0
    end where
    call sorti(eval, sort_idx, "descend")
    sum_eig2=sum(eval **2)
    sum_head=0
    do i=1, n
        sum_head=sum_head+eval(sort_idx(i))**2
        ratio=sum_head/sum_eig2
        TRACE4("top ",i," eigvalues contribute:",ratio)
        if(ratio >= eig2sumratio) exit
    end do
    final_nmodes=i
    nmodes=final_nmodes
    eigvalue(1:nmodes)=eval(sort_idx(1:nmodes))
    C(:,1:nmodes)=C(:,sort_idx(1:nmodes))
end subroutine decompose_C_syevd  !!}}}

    !!compare the Reconstructed Single-Point Correlation with True Single-Point Corr.
    !!write to lman_check.txt
    subroutine check_Cxyz_outfile(lman, xl0, yl0, Cz)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="check_Cxyz_outfile"
        type(type_localman) :: lman
        integer :: xidx=10, yidx=10, zidx=2
        integer :: xnumber=6, ynumber=6 !!one side: xidx, xidx+1, ..xidx+xnumber-1
        integer :: znumber=3  !!left side & right side
        integer,allocatable :: xindices(:), yindices(:), zindices(:)
        real :: xl0(:,:), yl0(:), Cz(:,:), dis
        real,allocatable :: xsingcorr(:,:,:), ysingcorr(:,:), zsingcorr(:)
        real,allocatable :: true_xsc(:,:,:), true_ysc(:,:), true_zsc(:)
        character(len=30) :: fmat
        integer :: i,j,k
        allocate(xindices(xnumber), yindices(ynumber), zindices(znumber))
        do i=1, xnumber
            xindices(i)=xidx+i-1
        end do
        do i=1, ynumber
            yindices(i)=yidx+i-1
        end do
        do i=1, znumber
            zindices(i)=i
        end do
        allocate(xsingcorr(xnumber,lman%my,lman%mz), &
                 ysingcorr(ynumber,lman%mz), zsingcorr(znumber))
        call check_Cx_singlepoint(lman, xidx, xindices, xsingcorr)
        call check_Cy_singlepoint(lman, yidx, yindices, ysingcorr)
        call check_Cz_singlepoint(lman, zidx, zindices, zsingcorr)
        !!call true SingCorr
        allocate(true_xsc(xnumber,lman%my,lman%mz), &
                 true_ysc(ynumber,lman%mz), true_zsc(znumber))
        do k=1,lman%mz
            do j=1,lman%my
                if(xl0(j,k)<0) then
                    true_xsc(:,j,k)=1.
                else
                    do i=1,xnumber
                        dis=abs(lman%xcoord(xindices(i))-lman%xcoord(xidx))
                        true_xsc(i,j,k)=getCorr_Gauss(dis/xl0(j,k))
                    end do
                end if
            end do
            do i=1,ynumber
                dis=abs(lman%ycoord(yindices(i))-lman%ycoord(yidx))
                true_ysc(i,k)=getCorr_Gauss(dis/yl0(k))
            end do
        end do
        true_zsc=Cz(zidx,zindices)
        !!printout
        open(11,file="lman_check.txt")
        write(11,*) "[XM_CHECK]"
        write(fmat,*) "(A,",xnumber,"F9.5)"
        do k=1, lman%mz
            write(11,"(A,I1.1)") "-iz=",k
            do j=1, lman%my
                write(11,"(A,I2.2)") "--iy=",j
                write(11,fmat) "T  ",true_xsc(:,j,k)
                write(11,fmat) "R  ",xsingcorr(:,j,k)
                write(11,fmat) "R-T",xsingcorr(:,j,k)-true_xsc(:,j,k)
            end do
        end do
        write(11,*) "[YM_CHECK]"
        write(fmat,*) "(A,",ynumber,"F9.5)"
        do k=1, lman%mz
            write(11,"(A,I1.1)") "-iz=",k
            write(11,fmat) "T  ",true_ysc(:,k)
            write(11,fmat) "R  ",ysingcorr(:,k)
            write(11,fmat) "R-T",ysingcorr(:,k)-true_ysc(:,k)
        end do
        write(11,*) "[ZM_CHECK]"
        write(fmat,*) "(A,",znumber,"F9.5)"
        write(11,fmat) "T  ", true_zsc(:)
        write(11,fmat) "R  ", zsingcorr(:)
        write(11,fmat) "R-T", zsingcorr(:)-true_zsc(:)
        close(11)
        deallocate(xsingcorr,ysingcorr,zsingcorr,true_xsc,true_ysc,true_zsc,&
                   xindices, yindices, zindices)
    end subroutine check_Cxyz_outfile  !!}}}

    !!give localman reconstructed single point corr (Lman*lman) , single point x is defined by idx
    !!other points is defined by indices, return singcorr
    subroutine check_Cx_singlepoint(lman, idx, indices, singcorr)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="check_Cx_singlepoint"
        type(type_localman) :: lman
        integer :: idx, indices(:), j, k
        real :: singcorr(:,:,:) !!(size(indices), my, mz)
        real,allocatable :: Cx(:,:)
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_xm)
        ASSUREX(all(shape(singcorr) == (/size(indices), lman%my, lman%mz/)))
        allocate(Cx(lman%nx,lman%nx))
        do k=1, lman%mz
            do j=1, lman%my
                Cx=matmul(lman%x3dm(j,k)%modes, transpose(lman%x3dm(j,k)%modes))
                singcorr(:,j,k)=Cx(idx,indices)
            end do
        end do
        deallocate(Cx)
    end subroutine check_Cx_singlepoint  !!}}}

    subroutine check_Cy_singlepoint(lman, idx, indices, singcorr)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="check_Cy_singlepoint"
        type(type_localman) :: lman
        integer :: idx, indices(:), k
        real :: singcorr(:,:) !!(size(indices), mz)
        real,allocatable :: Cy(:,:)
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_ym)
        ASSUREX(all(shape(singcorr) == (/size(indices), lman%mz/)))
        allocate(Cy(lman%ny,lman%ny))
        do k=1, lman%mz
            Cy=matmul(lman%y2dm(k)%modes, transpose(lman%y2dm(k)%modes))
            singcorr(:,k)=Cy(idx,indices)
        end do
        deallocate(Cy)
    end subroutine check_Cy_singlepoint  !!}}}

    subroutine check_Cz_singlepoint(lman, idx, indices, singcorr)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="check_Cz_singlepoint"
        type(type_localman) :: lman
        integer :: idx, indices(:)
        real :: singcorr(:) !!(size(indices))
        real,allocatable :: Cz(:,:)
        ASSUREX(lman%has_setdomain)
        ASSUREX(lman%has_set_zm)
        ASSUREX(all(shape(singcorr) == (/size(indices)/)))
        allocate(Cz(lman%ny,lman%ny))
        Cz=matmul(lman%zm%modes, transpose(lman%zm%modes))
        singcorr(:)=Cz(idx,indices)
        deallocate(Cz)
    end subroutine check_Cz_singlepoint  !!}}}

function getCorr_Gauss(dis)  !!{{{
    implicit none
    character(len=*),parameter :: PROCEDURE_NAME="getCorr_Gauss"
    real :: dis, getCorr_Gauss
    getCorr_Gauss=exp(-dis**2/2)
end function getCorr_Gauss  !!}}}

end program 


