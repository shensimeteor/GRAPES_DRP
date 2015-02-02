#include "inc_common.fpp"
module module_localtrue2dman !!localt2dm_, lt2dm
use utility
use module_types
use module_local
use lapack95
use f95_precision
use sort
implicit none

    ! variables
    real,parameter :: EARTH_R=6400 !!km
    
    public :: localt2dm_set_domain, localt2dm_generate_decompC2d, localt2dm_write_binary
contains

    ! methods
    subroutine localt2dm_set_domain(lt2dm, xcoord, ycoord)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localt2dm_set_domain"
        real :: xcoord(:), ycoord(:)
        type(type_localtrue2dman) :: lt2dm
        ASSUREX( .not. lt2dm%has_setdomain )
        lt2dm%nx=size(xcoord)
        lt2dm%ny=size(ycoord)
        allocate(lt2dm%xcoord(lt2dm%nx))
        allocate(lt2dm%ycoord(lt2dm%ny))
        lt2dm%xcoord=xcoord
        lt2dm%ycoord=ycoord
        lt2dm%has_setdomain=.true.
    end subroutine localt2dm_set_domain  !!}}}

    !!construct C2d, decompose C2d, select modes into lt2dm
    !!by Gaussian, L0 in km, decompose use syevd(may be too slow for big C)
    subroutine localt2dm_generate_decompC2d(lt2dm, L0, op_nmodes, op_eig2sumratio)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localt2dm_generate_decompC2d"
        type(type_localtrue2dman) :: lt2dm
        real :: L0, ratio
        integer,  optional :: op_nmodes
        real, optional :: op_eig2sumratio
        integer :: n, stat, final_nmodes,i 
        real,allocatable :: C(:,:), eigvalue(:)
        integer,allocatable :: sort_idx(:)
        ASSUREX(lt2dm%has_setdomain)
        n=lt2dm%nx * lt2dm%ny
        allocate(C(n,n),eigvalue(n),sort_idx(n))
        !!construct C2d
        call getC2d(lt2dm%xcoord, lt2dm%ycoord, C, L0)
        !!decompose C2d
        call syevd(C, eigvalue, jobz='V', uplo='L', info=stat)
        LPK_CHECK(stat)
        where(eigvalue <0) 
            eigvalue=0
        end where
        call sorti(eigvalue, sort_idx, "descend")
        print*,eigvalue(sort_idx(1:100))
        !!select modes
        if(present(op_nmodes) .and. .not. present(op_eig2sumratio)) then
            final_nmodes=op_nmodes
        else if(present(op_eig2sumratio) .and. .not. present(op_nmodes))then
            call ratio_2_nmodes(eigvalue(sort_idx), op_eig2sumratio, final_nmodes)
        else if(present(op_eig2sumratio) .and. present(op_nmodes)) then
            call ratio_2_nmodes(eigvalue(sort_idx), op_eig2sumratio, final_nmodes)
            if(final_nmodes>op_nmodes) final_nmodes=op_nmodes
        else
            final_nmodes=n
        end if
        ratio=sum(eigvalue(sort_idx(1:final_nmodes))**2)/sum(eigvalue**2)
        RP_LOG("select "//trim(tostring(ival=final_nmodes))//" modes, contribute "//trim(tostring(fval=ratio)))
        !!assign modes
        call localmodes_init(lt2dm%lm2d, n, final_nmodes)
        do i=1, final_nmodes
            lt2dm%lm2d%modes(:,i)=C(:,sort_idx(i))*sqrt(eigvalue(sort_idx(i)))
        end do
    end subroutine localt2dm_generate_decompC2d  !!}}}

    !! lt2dm is too larg to write in ascii: out: prefix.head in ascii; prefix.modes in binary
    subroutine localt2dm_write_binary(lt2dm, file_prefix)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="localt2dm_write_binary"
        type(type_localtrue2dman) :: lt2dm
        character(len=*) :: file_prefix
        character(len=100) :: strfmt
        integer :: i
        open(21, file=trim(file_prefix)//".head")
        write(21, "(3(I,','))") lt2dm%nx, lt2dm%ny, lt2dm%lm2d%m
        write(strfmt, *) '(',lt2dm%nx,"(F,','))"
        write(21,trim(strfmt)) lt2dm%xcoord
        write(strfmt, *) '(',lt2dm%ny,"(F,','))"
        write(21,trim(strfmt)) lt2dm%ycoord
        close(21)
        open(22, file=trim(file_prefix)//".modes", form="unformatted", access="direct", recl=lt2dm%nx*lt2dm%ny)
        do i=1, lt2dm%lm2d%m
            write(22, rec=i) lt2dm%lm2d%modes(:,i)
        end do
        close(22)
    end subroutine localt2dm_write_binary  !!}}}

!! private functions
    !!construct C2d (a homogeneous & isotropic correlation model), L0 in km
    subroutine getC2d(lon, lat, C2d, L0)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="getC2d"
        real :: lon(:), lat(:), C2d(:,:)
        real :: L0, dis
      !!  real :: getCorr_Gauss
        integer :: nx, ny, n, ix, iy, jx, jy, idx, jdx
        nx=size(lon)
        ny=size(lat)
        n=size(C2d,1)
        ASSUREX(n==nx*ny)
        C2d=-99.
        do iy=1, ny
            do ix=1, nx
                idx=(iy-1)*nx+ix
                do jy=1, ny
                    do jx=1, nx
                        jdx=(jy-1)*nx+jx
                        if(idx < jdx) then
                            call sphere_distance(lon(ix),lat(iy),lon(jx),lat(jy),EARTH_R,dis)
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

    !!pass in: eigvalue : sorted; eig2sumratio
    !!return: nmodes
    subroutine ratio_2_nmodes(eigvalue, eig2sumratio, nmodes)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="ratio_2_nmodes"
        real :: eigvalue(:), eig2sumratio
        integer :: nmodes
        integer :: i, n
        real :: sum_eig2, sum_head, ratio
        sum_eig2=sum(eigvalue**2)
        sum_head=0
        n=size(eigvalue)
        do i=1, n
            sum_head=sum_head+eigvalue(i)**2
            ratio=sum_head/sum_eig2
            print*,i,ratio
            if(ratio >= eig2sumratio) exit
        end do
        nmodes=i
    end subroutine ratio_2_nmodes  !!}}}

end module module_localtrue2dman

