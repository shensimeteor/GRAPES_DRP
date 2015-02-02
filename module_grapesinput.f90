#include "inc_common.fpp"
module module_grapesinput
use utility
  implicit none

  ! variables
  type :: type_grapesinput
    integer :: NX, NY, NZ, ND !! the configuration
    !! the data
    integer :: mdate(5)
    real,allocatable :: dthdz(:,:,:), thref(:,:,:), piref(:,:,:), zsx(:,:,:), zsy(:,:,:), &
                        zz(:,:,:), pip(:,:,:), pi(:,:,:), u(:,:,:), v(:,:,:), &
                        wzet(:,:,:), thp(:,:,:), th(:,:,:), moist_2(:,:,:), resv(:,:,:)
    real,allocatable :: tslb(:,:,:) ,smois(:,:,:) 
    real,allocatable :: zs(:,:), zst(:,:), fisx(:,:), fisy(:,:), tsk(:,:), &
                        ht(:,:), tmn(:,:), xland(:,:), lu_index(:,:), snowc(:,:), &
                        xice(:,:), xlat(:,:), xlong(:,:), soil_type(:,:), vegfrac(:,:)
    real,allocatable :: dx(:), dy(:), fu(:), fv(:), cosv(:), &
                        cosfi(:), tanu(:), tanv(:) !! of NY
    real,allocatable :: zh(:), dk(:), d2k(:) !! of NZ
    real,allocatable :: zzs(:), dzs(:)  !! of ND
    real :: ztop
  end type 

contains

  ! methods
 

  subroutine alloc_grapesinput(gidata, nx, ny, nz, nd) !!{{{
    type(type_grapesinput) :: gidata
    integer :: nx, ny, nz, nd
    gidata%NX=nx
    gidata%NY=ny
    gidata%NZ=nz
    gidata%ND=nd
    allocate(  gidata%dthdz(nx,ny,nz), gidata%thref(nx,ny,nz), gidata%piref(nx,ny,nz), gidata%zsx(nx,ny,nz), gidata%zsy(nx,ny,nz), &
               gidata%zz(nx,ny,nz), gidata%pip(nx,ny,nz), gidata%pi(nx,ny,nz), gidata%u(nx,ny,nz), gidata%v(nx,ny,nz), &
               gidata%wzet(nx,ny,nz), gidata%thp(nx,ny,nz), gidata%th(nx,ny,nz), gidata%moist_2(nx,ny,nz), gidata%resv(nx,ny,nz) )
    allocate(  gidata%tslb(nx,ny,nd), gidata%smois(nx,ny,nd))
    allocate(  gidata%zs(nx,ny), gidata%zst(nx,ny), gidata%fisx(nx,ny), gidata%fisy(nx,ny), gidata%tsk(nx,ny), &
               gidata%ht(nx,ny), gidata%tmn(nx,ny), gidata%xland(nx,ny), gidata%lu_index(nx,ny), gidata%snowc(nx,ny), &
               gidata%xice(nx,ny), gidata%xlat(nx,ny), gidata%xlong(nx,ny), gidata%soil_type(nx,ny), gidata%vegfrac(nx,ny) )
    allocate(  gidata%dx(ny), gidata%dy(ny), gidata%fu(ny), gidata%fv(ny), gidata%cosv(ny), &
               gidata%cosfi(ny), gidata%tanu(ny), gidata%tanv(ny) )
    allocate(  gidata%zh(nz), gidata%dk(nz), gidata%d2k(nz), gidata%zzs(nd), gidata%dzs(nd) )
  end subroutine alloc_grapesinput  !!}}}

  subroutine dealloc_grapesinput(gidata)  !!{{{
  implicit none
      character(len=*),parameter :: PROCEDURE_NAME="dealloc_grapesinput"
      deallocate(  gidata%dthdz, gidata%thref, gidata%piref, gidata%zsx, gidata%zsy, &
                 gidata%zz, gidata%pip, gidata%pi, gidata%u, gidata%v, &
                 gidata%wzet, gidata%thp, gidata%th, gidata%moist_2, gidata%resv )
      deallocate(  gidata%tslb, gidata%smois)
      deallocate(  gidata%zs, gidata%zst, gidata%fisx, gidata%fisy, gidata%tsk, &
                 gidata%ht, gidata%tmn, gidata%xland, gidata%lu_index, gidata%snowc, &
                 gidata%xice, gidata%xlat, gidata%xlong, gidata%soil_type, gidata%vegfrac )
      deallocate(  gidata%dx, gidata%dy, gidata%fu, gidata%fv, gidata%cosv, &
                 gidata%cosfi, gidata%tanu, gidata%tanv )
      deallocate(  gidata%zh, gidata%dk, gidata%d2k, gidata%zzs(nd), gidata%dzs(nd) )
      gidata%NX=0
      gidata%NY=0
      gidata%NZ=0
      gidata%ND=0
  end subroutine dealloc_grapesinput  !!}}}

  subroutine read_grapesinput(gidata, filepath) !!{{{
    type(type_grapesinput) :: gidata
    character(len=*) :: filepath
    integer :: nx, ny, nz, nd, k
    nx=gidata%NX
    ny=gidata%NY
    nz=gidata%NZ
    nd=gidata%ND
    open(10,file=trim(filepath),form='unformatted',convert="big_endian", status='old')
    read(10) gidata%mdate
    do k=1, nz
      read(10) gidata%dthdz(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%thref(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%piref(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%zsx(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%zsy(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%zz(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%pip(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%pi(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%u(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%v(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%wzet(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%thp(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%th(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%moist_2(:,:,k)
    enddo
    do k=1, nz
      read(10) gidata%resv(:,:,k)
    enddo
    !!
    do k=1, nd
      read(10) gidata%tslb(:,:,k)
    end do
    do k=1, nd
      read(10) gidata%smois(:,:,k)
    end do
    !!
    read(10) gidata%zs 
    read(10) gidata%zst
    read(10) gidata%fisx
    read(10) gidata%fisy
    read(10) gidata%tsk
    read(10) gidata%ht
    read(10) gidata%tmn
    read(10) gidata%xland
    read(10) gidata%lu_index
    read(10) gidata%snowc
    read(10) gidata%xice
    read(10) gidata%xlat
    read(10) gidata%xlong
    read(10) gidata%soil_type
    read(10) gidata%vegfrac
    !! 
    read(10) gidata%dx
    read(10) gidata%dy
    read(10) gidata%fu
    read(10) gidata%fv
    read(10) gidata%cosv
    read(10) gidata%cosfi
    read(10) gidata%tanu
    read(10) gidata%tanv
    !!
    read(10) gidata%zh
    read(10) gidata%dk
    read(10) gidata%d2k
    read(10) gidata%zzs
    read(10) gidata%dzs
    !!
    read(10) gidata%ztop 
   
    close(10)
  end subroutine read_grapesinput  !!}}}

  subroutine write_grapesinput(gidata, filepath) !!{{{
    type(type_grapesinput) :: gidata
    character(len=*) :: filepath
    integer :: nx, ny, nz, nd, k
    nx=gidata%NX
    ny=gidata%NY
    nz=gidata%NZ
    nd=gidata%ND
    open(10,file=trim(filepath),form='unformatted',convert="big_endian")
    write(10) gidata%mdate
    do k=1, nz
      write(10) gidata%dthdz(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%thref(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%piref(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%zsx(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%zsy(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%zz(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%pip(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%pi(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%u(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%v(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%wzet(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%thp(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%th(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%moist_2(:,:,k)
    enddo
    do k=1, nz
      write(10) gidata%resv(:,:,k)
    enddo
    !!
    do k=1, nd
      write(10) gidata%tslb(:,:,k)
    end do
    do k=1, nd
      write(10) gidata%smois(:,:,k)
    end do
    !!
    write(10) gidata%zs 
    write(10) gidata%zst
    write(10) gidata%fisx
    write(10) gidata%fisy
    write(10) gidata%tsk
    write(10) gidata%ht
    write(10) gidata%tmn
    write(10) gidata%xland
    write(10) gidata%lu_index
    write(10) gidata%snowc
    write(10) gidata%xice
    write(10) gidata%xlat
    write(10) gidata%xlong
    write(10) gidata%soil_type
    write(10) gidata%vegfrac
    !! 
    write(10) gidata%dx
    write(10) gidata%dy
    write(10) gidata%fu
    write(10) gidata%fv
    write(10) gidata%cosv
    write(10) gidata%cosfi
    write(10) gidata%tanu
    write(10) gidata%tanv
    !!
    write(10) gidata%zh
    write(10) gidata%dk
    write(10) gidata%d2k
    write(10) gidata%zzs
    write(10) gidata%dzs
    !!
    write(10) gidata%ztop 
   
    close(10)
  end subroutine write_grapesinput  !!}}}

    !!bin, directly output: pi/piref/pip, th/thref/thp, moist_2 u/v/zz
    subroutine reformat_grapesinput(gidata, filepath)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="reformat_grapesinput"
        type(type_grapesinput) :: gidata
        character(len=*) :: filepath
        integer :: i,cnt
        open(10,file=trim(filepath),form='unformatted',access="direct",recl=gidata%NX*gidata%NY)
        cnt=0
        do i=1,gidata%NZ
            cnt=cnt+1
            write(10,rec=cnt) gidata%pi(:,:,i)
        end do
        do i=1,gidata%NZ
            cnt=cnt+1
            write(10,rec=cnt) gidata%piref(:,:,i)
        end do
        do i=1,gidata%NZ
            cnt=cnt+1
            write(10,rec=cnt) gidata%pip(:,:,i)
        end do
        do i=1,gidata%NZ
            cnt=cnt+1
            write(10,rec=cnt) gidata%th(:,:,i)
        end do
        do i=1,gidata%NZ
            cnt=cnt+1
            write(10,rec=cnt) gidata%thref(:,:,i)
        end do
        do i=1,gidata%NZ
            cnt=cnt+1
            write(10,rec=cnt) gidata%thp(:,:,i)
        end do
        do i=1,gidata%NZ
            cnt=cnt+1
            write(10,rec=cnt) gidata%moist_2(:,:,i)
        end do
        do i=1,gidata%NZ
            cnt=cnt+1
            write(10,rec=cnt) gidata%u(:,:,i)
        end do
        do i=1,gidata%NZ
            cnt=cnt+1
            write(10,rec=cnt) gidata%v(:,:,i)
        end do
        do i=1,gidata%NZ
            cnt=cnt+1
            write(10,rec=cnt) gidata%zz(:,:,i)
        end do
        close(10)
    end subroutine reformat_grapesinput  !!}}}

    !!use HydronStatic Balance to calculate th
    subroutine CalTH_fromHydroStatic(gidata, pi, q, th)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="CalTH_fromHydroStatic"
        type(type_grapesinput) :: gidata
        real :: th(:,:,:), pi(:,:,:), q(:,:,:)
        real,allocatable :: zhalf(:,:,:)
        real,parameter :: cp = 1004.5, gravity  = 9.80665
        integer :: i,j,k
        if(.not. all(gidata%NX == (/size(th,1), size(q,1), size(pi,1)/))) then
            RP_ERROR("input parameter 1st dimension size conflicts")
        end if
        if(.not. all(gidata%NY == (/size(th,2), size(q,2), size(pi,2)/))) then
            RP_ERROR("input parameter 2st dimension size conflicts")
        end if
        if(.not. all(gidata%NZ == (/size(th,3)+1, size(q,3)+1, size(pi,3)/))) then
            RP_ERROR("input parameter 1st dimension size conflicts")
        end if
        !!get zhalf
        allocate(zhalf(gidata%NX, gidata%NY, gidata%NZ))
        zhalf(:,:,1)=gidata%zz(:,:,2)-0.5*(gidata%zz(:,:,3)-gidata%zz(:,:,2))
        zhalf(:,:,gidata%NZ)=gidata%zz(:,:,gidata%NZ)+0.5*(gidata%zz(:,:,gidata%NZ)- &
                             gidata%zz(:,:,gidata%NZ-1))
        do k=2, gidata%NZ-1 !!2-37
            zhalf(:,:,k)=(gidata%zz(:,:,k)+gidata%zz(:,:,k+1))*0.5
        end do
        !!get th
        do k=1, gidata%NZ-1 !! full level 1-37
            th(:,:,k)=-gravity * (zhalf(:,:,k+1)-zhalf(:,:,k)) / &
                     (cp*(1+0.608*q(:,:,k))*(pi(:,:,k+1)-pi(:,:,k)))
        end do
        th(:,:,1)=th(:,:,2)
    end subroutine CalTH_fromHydroStatic  !!}}}

end module module_grapesinput
