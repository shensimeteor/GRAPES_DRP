! --------------------------------------------------------------------------
! Description:
!
!   something about the model state(modeldesc), model dependent (rewrite if model changes)
!   a interface to exchange data from file(structed) to memory(array)
!   interface convention:
!       1. modeldesc_init(md)
!       2. modeldesc_getcoord(md, ncfile)
!       3. modelstate_input(md, ncfile, x)
!       4. modelstate_output(md, ncfile, x)
!       5. modeldesc_final(md)
!
! History:
!
!   2014-07-02:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
module module_modelstate
use module_types
use utility
use netcdf
use module_grapesinput
implicit none

    
contains

    ! methods
    subroutine modeldesc_init(md)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="modeldesc_init"
        type(type_modeldesc) :: md
        integer :: i, idx_cnt
        md%n_xcoord=2
        md%n_ycoord=2
        md%n_zcoord=2
        allocate(md%xcoords(md%n_xcoord))
        allocate(md%ycoords(md%n_ycoord))
        allocate(md%zcoords(md%n_zcoord))
        md%nvar=4 !! u,v,pi,hum
        allocate(md%vars(md%nvar))
        idx_cnt=0
        md%vars(1)%varname="u"
        md%vars(1)%n_grids=360*180*36 !!u
        md%vars(1)%xcoord_id=2
        md%vars(1)%ycoord_id=1
        md%vars(1)%zcoord_id=3
        md%vars(1)%ibgn=1
        md%vars(1)%iend=md%vars(1)%n_grids
        md%vars(2)%varname="v"
        md%vars(2)%n_grids=360*179*36 !!v: pole will not analyze
        md%vars(2)%xcoord_id=1
        md%vars(2)%ycoord_id=2
        md%vars(2)%zcoord_id=3
        md%vars(2)%ibgn=md%vars(1)%iend+1
        md%vars(2)%iend=md%vars(1)%iend+md%vars(2)%n_grids
        md%vars(3)%varname="pi"
        md%vars(3)%n_grids=360*180*38  !!pi
        md%vars(3)%xcoord_id=1
        md%vars(3)%ycoord_id=1
        md%vars(3)%zcoord_id=1
        md%vars(3)%ibgn=md%vars(2)%iend+1
        md%vars(3)%iend=md%vars(2)%iend+md%vars(3)%n_grids
        md%vars(4)%varname="hum"
        md%vars(4)%n_grids=360*180*37  !!hum
        md%vars(4)%xcoord_id=1
        md%vars(4)%ycoord_id=1
        md%vars(4)%zcoord_id=2
        md%vars(4)%ibgn=md%vars(2)%iend+1
        md%vars(4)%iend=md%vars(2)%iend+md%vars(3)%n_grids
        md%n_total=md%vars(4)%iend
        TRACEX(md%n_total)
        TRACEX(md%vars(:)%ibgn)
        TRACEX(md%vars(:)%iend)
        md%has_init=.true.
    end subroutine modeldesc_init  !!}}}

    subroutine modeldesc_final(md)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="modeldesc_final"
        type(type_modeldesc) :: md
        integer :: i
        md%has_init=.false.
        deallocate(md%vars)
        deallocate(md%xcoords, md%ycoords, md%zcoords)
        md%nvar=0
        md%n_xcoord=0
        md%n_ycoord=0
        md%n_zcoord=0
        md%n_total=0
    end subroutine modeldesc_final  !!}}}

    subroutine modeldesc_getcoord(md)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="modeldesc_getcoord"
        type(type_modeldesc) :: md
        integer :: ncid, varid
        real :: zh(37) !! full level zcoord(zzet) values
        real :: dz1, dz37
        integer :: i
        !!x
        md%xcoords(1)%n=360
        allocate(md%xcoords(1)%val(md%xcoords(1)%n))
        call fspan(0.,359.,360,md%xcoords(1)%val)
        md%xcoords(1)%iscycle = .true. 
        md%xcoords(2)%n=360
        allocate(md%xcoords(2)%val(md%xcoords(2)%n))
        call fspan(0.5,359.5,360,md%xcoords(2)%val)
        md%xcoords(2)%iscycle = .true. 
        !!y
        md%ycoords(1)%n=180
        allocate(md%ycoords(1)%val(md%ycoords(1)%n))
        call fspan(-89.5,89.5,180,md%ycoords(1))
        md%ycoords(1)%iscycle = .false.
        md%ycoords(2)%n=179
        call fspan(-89,89,179,md%ycoords(2))
        allocate(md%ycoords(2)%val(md%ycoords(2)%n))
        md%ycoords(2)%iscycle = .false.
        !!todo :zcoords :
        zh=(/0, 65.24981, 230.1493, 474.9776, 784.4373, 1147.267, 1554.654, 1999.950, 2478.427, 2986.432, &
           3521.796, 4083.160, 4670.249, 5283.126, 5922.162, 6588.492, 7283.020, 8006.647, 8759.948, &
           9543.588, 10357.67, 11202.07, 12078.03, 12986.87, 13935.47, 14995.40, 16194.72, 17627.36, &
           19000.00, 20500.00, 22000.00, 23500.00, 25300.00, 27000.00, 29000.00, 30700.00, 32500.00/)
        md%zcoords(1)%n=37 !!full level
        allocate(md%zcoords(1)%val(md%zcoords(1)%n))
        md%zcoords(1)%val=zh
        md%zcoords(1)%iscycle=.false.
        md%zcoords(2)%n=38 !!0-37 half level
        allocate(md%zcoords(2)%val(md%zcoords(2)%n))
        dz1=zh(2)-zh(1)
        dz37=zh(37)-zh(36)
        md%zcoords(2)%val(1)=zh(1)-0.5*dz1
        md%zcoords(2)%val(38)=zh(37)+0.5*dz37
        do i=2,37
            md%zcoords(2)%val(i)=0.5*(zh(i-1)+zh(i))
        end do
        md%zcoords(2)%iscycle=.false.
        md%zcoords(3)%n=36 !!1-36 half level
        allocate(md%zcoords(3)%val(md%zcoords(3)%n))
        do i=1,36
            md%zcoords(3)%val(i)=0.5*(zh(i)+zh(i+1))
        end do
        md%zcoords(3)%iscycle=.false.
    end subroutine modeldesc_getcoord  !!}}}

    subroutine modelstate_input(md, gifile, x)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="modelstate_input"
        type(type_modeldesc) :: md
        character(len=*) :: gifile
        type(type_grapesinput) :: gidata
        real :: x(:)
        integer :: ncid, varid, i
        if(size(x) /= md%n_total) then
            RP_ERROR("size of x /= md%n_total")
        end if
        call alloc_grapesinput(gidata, 360, 180, 38, 5)
        call read_grapesinput(gidata, trim(gifile))
        x(md%vars(1)%ibgn: md%vars(1)%iend) = reshape(gidata%u(:,:,2:37), (/360*180*36/))
        x(md%vars(2)%ibgn: md%vars(2)%iend) = reshape(gidata%v(:,2:180,2:37), (/360*179*36/))
        x(md%vars(3)%ibgn: md%vars(3)%iend) = reshape(gidata%pi, (/360*180*38/))
        x(md%vars(4)%ibgn: md%vars(4)%iend) = reshape(gidata%hum, (/360*180*37/))
        call dealloc_grapesinput(gidata)
    end subroutine modelstate_input  !!}}}

    !!gifile (a copy of xb) must exist to provide other variables of x
    subroutine modelstate_output(md, gifile, x)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="modelstate_output"
        type(type_modeldesc) :: md
        character(len=*) :: gifile
        type(type_grapesinput) :: gidata
        real :: x(:)
        real,allocatable :: th(:,:,:)
        integer :: ncid, varid, i, lexist, nx, ny, nz
        if(size(x) /= md%n_total) then
            RP_ERROR("size of x /= md%n_total")
        end if
        inquire(file=trim(gifile), exist=lexist)
        if(.not. lexist) then
            RP_ERROR("template file "//trim(gifile)//" NOT exist")
        end if
        call alloc_grapesinput(gidata, 360, 180, 38, 5)
        call read_grapesinput(gidata, trim(gifile))
        gidata%u(:,:,2:37)=reshape(x(md%vars(1)%ibgn:md%vars(1)%iend), (/360,180,36/))
        gidata%v(:,2:180,2:37)=reshape(x(md%vars(2)%ibgn:md%vars(2)%iend), (/360,179,36/))
        gidata%v(:,1,:)=0. !!south pole
        gidata%pi(:,:,:)=reshape(x(md%vars(3)%ibgn:md%vars(3)%iend), (/360,180,38/))
        gidata%moist_2(:,:,2:38)=reshape(x(md%vars(4)%ibgn:md%vars(3)%iend), (/360,180,37/))
        allocate(th(360,180,37))
        call CalTH_fromHydroStatic(gidata, gidata%pi, gidata%moist_2(:,:,2:38), th) 
        gidata%th(:,:,2:38)=th
        gidata%pip=gidata%pi-gidata%piref
        gidata%thp=gidata%th-gidata%thref
        call write_grapesinput(gidata, trim(gifile))
        call dealloc_grapesinput(gidata)
        deallocate(th)
    end subroutine modelstate_output  !!}}}

    function modelstate_get_idx(md, ivar, ix, iy, iz)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="modelstate_get_idx"
        type(type_modeldesc) :: md
        integer :: ivar, ix, iy, iz, modelstate_get_idx
        integer :: nx, ny, nz, nvar
        nvar=md%nvar
        ASSUREX( ivar<=nvar .and. ivar >=1 )
        nx=md%xcoords(md%vars(ivar)%xcoord_id)%n
        ny=md%ycoords(md%vars(ivar)%ycoord_id)%n
        nz=md%zcoords(md%vars(ivar)%zcoord_id)%n
        ASSUREX( ix<=nx .and. ix >=1)
        ASSUREX( iy<=ny .and. iy >=1)
        ASSUREX( iz<=nz .and. iz >=1)
        modelstate_get_idx=md%vars(ivar)%ibgn-1+(iz-1)*nx*ny+(iy-1)*nx+ix
    end function modelstate_get_idx  !!}}}

    !! if offset > 0, move towards east
    subroutine modelstate_xshift(md, x, newx, offset)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="modelstate_xshift"
        type(type_modeldesc) :: md
        real :: x(:), newx(:)
        integer :: offset, i, nx, j
        real,allocatable :: field(:,:,:)
        integer,allocatable :: idx(:)
        do i=1, md%nvar
            call modelstate_alloc3dfield(md, i, field)
            nx=size(field,1)
            allocate(idx(nx))
            idx=(/(j-offset, j=1,nx)/)
            where(idx <= 0) 
                idx=idx+nx
            end where
            where(idx > nx)
                idx=idx-nx
            end where
            field=field(idx,:,:)
            newx(md%vars(i)%ibgn:md%vars(i)%iend)=reshape(field, (/size(field)/))
            deallocate(idx, field)
        end do
    end subroutine modelstate_xshift  !!}}}

    subroutine modelstate_alloc3dfield(md, ivar, field)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="modelstate_alloc3dfield"
        integer :: ivar
        type(type_modeldesc) :: md
        real,allocatable :: field(:,:,:)
        integer :: nx, ny, nz, ix, iy, iz
        ASSUREX( ivar<=md%nvar .and. ivar >=1 )
        ix=md%vars(ivar)%xcoord_id
        iy=md%vars(ivar)%ycoord_id
        iz=md%vars(ivar)%zcoord_id
        nx=md%xcoords(ix)%n
        ny=md%ycoords(iy)%n
        nz=md%zcoords(iz)%n
        allocate(field(nx,ny,nz))
    end subroutine modelstate_alloc3dfield  !!}}}

    subroutine modelstate_get3dfield(md, ivar, x, f3d)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="modelstate_get3dfield"
        integer :: ivar
        type(type_modeldesc) :: md
        real :: f3d(:,:,:), x(:)
        integer :: nx, ny, nz, ix, iy, iz
        ASSUREX( ivar<=md%nvar .and. ivar >=1 )
        ix=md%vars(ivar)%xcoord_id
        iy=md%vars(ivar)%ycoord_id
        iz=md%vars(ivar)%zcoord_id
        nx=md%xcoords(ix)%n
        ny=md%ycoords(iy)%n
        nz=md%zcoords(iz)%n
        ASSUREX(all((/nx,ny,nz/)==shape(f3d)))
        f3d=reshape(x(md%vars(ivar)%ibgn:md%vars(ivar)%iend), (/nx,ny,nz/))
    end subroutine modelstate_get3dfield  !!}}}

end module module_modelstate
