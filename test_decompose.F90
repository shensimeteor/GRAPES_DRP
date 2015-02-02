#include 'inc_common.fpp'
program test_decompose
use utility
use lapack95
use f95_precision
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="test_decompose"
    integer :: nxc=90, nyc=46 !!coarse grid
   !! integer :: nxc=180, nyc=91
   !! integer :: nxc=60, nyc=31
    real :: dxc, dyc 
    integer :: nc, nmodes
    real :: L0=600
    real,allocatable :: lat(:), lon(:), C2d(:,:)
    real :: tbgn, t, tlast=-1
    real,allocatable :: eigvalue(:), eigvector(:,:)

    dxc =360./nxc
    dyc =180./(nyc-1)
    allocate(lon(nxc),lat(nyc))
    call fspan(0., 360.-dxc, nxc, lon)
    call fspan(-90., 90., nyc, lat)
    nc=nxc*nyc
    allocate(C2d(nc,nc), eigvalue(nc), eigvector(nc,nc))
    TRACEX(shape(C2d))
    TIME_START(tbgn)
    call getC2d(lon, lat, C2d, L0)
    TIME_CLICK(t, tbgn, tlast)
!!    call decompose_C2d(C2d, nc, eigvalue, eigvector)
    call decompose_C_syevd(C2d, eigvalue, 0.98, nmodes)
    print*,eigvalue(1:100)
    print*,nmodes
    TIME_CLICK(t, tbgn, tlast)
    deallocate(C2d)
contains

!!eigvalue & eigvector are order of size(C,1)
!!give :: eig2sumratio return nmodes
!!C: output the eigvector matrix
subroutine decompose_C_syevd(C, eigvalue, eig2sumratio, nmodes)  !!{{{
use sort
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

subroutine sphere_distance(lon1, lat1, lon2, lat2, R, dis)  !!{{{
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="sphere_distance"
    real,parameter :: PI=3.1415926, D2R=PI/180
    real :: lon1, lat1, lon2, lat2, R, dis
    real :: dlon,arc
    dlon=lon2-lon1
    arc=sin(lat1*D2R)*sin(lat2*D2R)+cos(lat1*D2R)*cos(lat2*D2R)*cos(dlon*D2R)
    if(arc>1 .and. arc<1.0001) then
        arc=1.
    else if(arc<-1 .and. arc>-1.0001) then
        arc=-1.
    end if
    arc=acos(arc)
    dis=R*arc
end subroutine sphere_distance  !!}}}

!!construct C2d (a homogeneous & isotropic correlation model), L0 in km
subroutine getC2d(lon, lat, C2d, L0)  !!{{{
use utility
implicit none
    character(len=*),parameter :: PROCEDURE_NAME="getC2d"
    real :: lon(:), lat(:), C2d(:,:)
    real :: L0, dis
  !!  real :: getCorr_Gauss
    real,parameter :: earth_r=6400
    integer :: nx, ny, n, ix, iy, jx, jy, idx, jdx
    nx=size(lon)
    ny=size(lat)
    n=size(C2d,1)
    ASSUREX(n==nx*ny)
    print*,C2d(1:5,1:5)
    C2d=-99.
    do iy=1, ny
        do ix=1, nx
            idx=(iy-1)*nx+ix
            do jy=1, ny
                do jx=1, nx
                    jdx=(jy-1)*nx+jx
                    if(idx < jdx) then
                        call sphere_distance(lon(ix),lat(iy),lon(jx),lat(jy),earth_r,dis)
                        C2d(idx,jdx)=getCorr_Gauss(dis/L0)
                    else if(idx == jdx) then
                        C2d(idx,jdx)=1
                    else
                        C2d(idx,jdx)=C2d(jdx,idx)
                    end if
                end do
            end do
        end do
    end do
!!test
    TRACEX(maxval(C2d))
    TRACEX(minval(C2d))
end subroutine getC2d  !!}}}

function getCorr_Gauss(dis)  !!{{{
    implicit none
    character(len=*),parameter :: PROCEDURE_NAME="getCorr_Gauss"
    real :: dis, getCorr_Gauss
    getCorr_Gauss=exp(-dis**2/2)
end function getCorr_Gauss  !!}}}


end program test_decompose
