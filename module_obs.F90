! --------------------------------------------------------------------------
! Description:
!
!   something on obs ( so model dependent): read/write obs; normalize
!   for real assimilations, this module can be treated as a obs interface, 
!   hinding obs details
!   todo(for WRF/GRAPES where filtered the QC flaged obs is needed
!   1. obs_read_odobs : read od & obs
!   2. obs_read_obs_byod : read obs by od given (only read the obs set by od, ie. ignore the filtered)
!   3. obs_filter_od : filter od ( whose qcflag is flaged) to a new_od
!   4. obs_filter_obs : filter y
!
! History:
!
!   2014-07-08:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
module module_obs
use utility
use module_types
implicit none
private
    ! variables
    
    public :: obs_read_txt, obsdesc_final, obs_normalize, obs_write_txt, obs_group, &
              obs_getqcflag, obs_filter_od, obs_filter_obs

    real :: MISSING_VALUE=-888888.00

contains

    ! methods
    !! read y & allocate/read od (if not read yet)
    !! <filenames> array corsponding the <tindices>, read od & y
    !! if any in <filenames> not exist, just ignore and continue reading others
    subroutine obs_read_txt(od, y, filenames, tindices)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="obs_read_txt"
        type(type_obsdesc) :: od
        real :: y(:)
        character(len=*) :: filenames(:)
        integer :: tindices(:)
        character(len=150) :: strline
        integer :: n_input, i, nline, n_obs, cnt, stat, temp
        logical :: lexist
        if(size(filenames)  /= size(tindices)) then
            RP_ERROR("INPUT filenames & tindices not the same")
        end if
        n_input=size(filenames)
        n_obs=0
        do i=1,n_input
            inquire(file=trim(filenames(i)), exist=lexist)
            if(lexist) then
                call query_lines(filenames(i), nline)
                n_obs=n_obs+nline-1
            end if
        end do
        if(size(y) < n_obs) then
            RP_ERROR("size(y) < real n_obs")
        end if
        if(.not. od%has_read) then
            call private_alloc_obsdesc(od, n_obs)
            od%qcflag(:)=0
        end if
        cnt=0
        do i=1,n_input
            inquire(file=trim(filenames(i)), exist=lexist)
            if( .not. lexist) cycle
            open(11, file=trim(filenames(i)), iostat=stat)
            if(stat/=0) then
                RP_ERROR("FAIL to open "//trim(filenames(i)))
            end if
            read(11, *, iostat=stat)
            if(stat/=0) then
                RP_ERROR("FAIL to read "//trim(filenames(i)))
            end if
            do 
                read(11, "(A)",iostat=stat) strline
                if(stat/=0) exit
                cnt=cnt+1
                od%tidx(cnt)=tindices(i)
                od%zpos(cnt)=1
                read(strline, *) od%varname(cnt), temp, &
                                          od%xpos(cnt), od%ypos(cnt), od%zpos(cnt), &
                                          y(cnt), od%error(cnt)
            end do
            close(11)
        end do
        if(cnt /= n_obs) then
            RP_ERROR("cnt /= n_obs")
        end if
        od%has_read=.true.
    end subroutine obs_read_txt  !!}}}

    subroutine obs_write_txt(od, y, file_prefix, nt_window)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="obs_write_txt"
        type(type_obsdesc) :: od
        real :: y(:)
        character(len=*) :: file_prefix
        integer :: nt_window
        character(len=50) :: filename
        integer :: i, cnt(nt_window)
        logical :: lcreate(nt_window)
        if(.not. od%has_read) then
            RP_ERROR("obsdesc not read yet!")
        end if
        lcreate=.false.
        cnt=0
        do i=1, od%n_obs
            if(.not. lcreate(od%tidx(i)+1)) then
                write(filename,"(A, I2.2)") trim(file_prefix), od%tidx(i)
                open(10+od%tidx(i)+1, file=trim(filename)) 
                write(10+od%tidx(i)+1, *) "var      posidx         lon         lat     pressure    value       error"
                lcreate(od%tidx(i)+1)=.true.
            end if
            cnt(od%tidx(i)+1)=cnt(od%tidx(i)+1)+1
!!            write(10+od%tidx(i)+1, *), od%varname(i), cnt(od%tidx(i)+1), od%xpos(i), od%ypos(i), od%zpos(i), y(i), od%error(i)
            write(10+od%tidx(i)+1,"(A,I,F,F,F,G,G)") trim(od%varname(i)), cnt(od%tidx(i)+1), od%xpos(i), od%ypos(i), od%zpos(i), y(i), od%error(i)
        end do
        do i=1, nt_window
            if(lcreate(i)) then
                close(10+i)
            end if
        end do
    end subroutine obs_write_txt  !!}}}

    !!give ys(:,:), for i if any of y(i,:) or od%xpos/ypos/zpos is missing, set qcflag=1
    subroutine obs_getqcflag(ys, od)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="obs_getqcflag"
        type(type_obsdesc) :: od
        real :: ys(:,:), diff
        integer :: nys, nobs, iobs, i
        diff=1.0
        nobs=size(ys,1)
        nys=size(ys,2)
        ASSUREX(nobs == od%n_obs)
        do iobs=1,nobs
            do i=1, nys
                if(abs(ys(i,iobs)-MISSING_VALUE)<diff) then
                    od%qcflag(iobs)=1
                    exit
                end if
            end do
            if(od%qcflag(iobs) == 0) then
                if(abs(od%xpos(iobs)-MISSING_VALUE)<diff) then
                    od%qcflag(iobs)=1
                else if(abs(od%ypos(iobs)-MISSING_VALUE)<diff) then
                    od%qcflag(iobs)=1
                else if(abs(od%zpos(iobs)-MISSING_VALUE)<diff) then
                    od%qcflag(iobs)=1
                end if
            end if
        end do
    end subroutine obs_getqcflag  !!}}}
    
    subroutine obs_filter_od(od, newod)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="obs_filter_od"
        type(type_obsdesc) :: od, newod
        integer :: newod_nobs, i, j
        newod_nobs=count(od%qcflag == 0)
        call private_alloc_obsdesc(newod, newod_nobs)
        j=0
        do i=1, od%n_obs
            if(od%qcflag(i)==0) then
                j=j+1
                newod%xpos(j)=od%xpos(i)
                newod%ypos(j)=od%ypos(i)
                newod%zpos(j)=od%zpos(i)
                newod%error(j)=od%error(i)
                newod%qcflag(j)=0
                newod%tidx(j)=od%error(i)
                newod%varname(j)=od%varname(i)
                newod%has_read=.true.
            end if
        end do
    end subroutine obs_filter_od  !!}}}

    subroutine obs_filter_obs(od, newod, y, newy)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="obs_filter_obs"
        type(type_obsdesc) :: od, newod
        real :: y(:), newy(:)
        integer :: i, j
        ASSUREX(size(y)==od%n_obs)
        ASSUREX(size(newy)==newod%n_obs)
        ASSUREX(od%has_read)
        ASSUREX(newod%has_read)
        j=0
        do i=1, od%n_obs
            if(od%qcflag(i)==0) then
                j=j+1
                newy(j)=y(i)
            end if
        end do
    end subroutine obs_filter_obs  !!}}}

!!!! MODEL DEPENDENT!!
    subroutine obs_group(od)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="obs_group"
        type(type_obsdesc) :: od
        integer :: tmp, i
        if(od%has_group) then
            if(od%has_group_set) then
                call obs_delete_group(od)
            end if
            od%group(OG_L2)%group_name="lev200"
            od%group(OG_L2)%n=count(od%zpos == 200.)
            allocate(od%group(OG_L2)%idx(od%group(OG_L2)%n))
            call getind_flt(od%zpos, 200., od%group(OG_L2)%idx, tmp, 0.1)
            od%group(OG_L5)%group_name="lev500"
            od%group(OG_L5)%n=count(od%zpos == 500.)
            allocate(od%group(OG_L5)%idx(od%group(OG_L5)%n))
            call getind_flt(od%zpos, 500., od%group(OG_L5)%idx, tmp, 0.1)
            od%group(OG_L8)%group_name="lev800"
            od%group(OG_L8)%n=count(od%zpos == 800.)
            allocate(od%group(OG_L8)%idx(od%group(OG_L8)%n))
            call getind_flt(od%zpos, 800., od%group(OG_L8)%idx, tmp, 0.1)
            od%has_group_set=.true.
        else
            RP_LOG("In obs_group, as you set od%has_group .false., no group will be set")
        end if
    end subroutine obs_group  !!}}}

    subroutine obs_delete_group(od)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="obs_delete_group"
        type(type_obsdesc) :: od
        integer :: i
        do i=1, od%n_group
            deallocate(od%group(i)%idx)
            od%group(i)%n=0
            od%group(i)%group_name=" "
        end do
        od%has_group_set=.false.
    end subroutine obs_delete_group  !!}}}

    subroutine obsdesc_final(od)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="obsdesc_final"
        type(type_obsdesc) :: od
        call private_dealloc_obsdesc(od)
    end subroutine obsdesc_final  !!}}}

    !! y -> normalized y
    subroutine obs_normalize(od, y)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="obs_normalize"
        real :: y(:)
        type(type_obsdesc) :: od
        y=y/od%error
    end subroutine obs_normalize  !!}}}


!! private
    subroutine private_alloc_obsdesc(od, n_obs)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="private_alloc_obsdesc"
        type(type_obsdesc) :: od
        integer :: n_obs
        od%n_obs=n_obs
        allocate(od%xpos(n_obs))
        allocate(od%ypos(n_obs))
        allocate(od%zpos(n_obs))
        allocate(od%error(n_obs))
        allocate(od%qcflag(n_obs))
        allocate(od%tidx(n_obs))
        allocate(od%varname(n_obs))
    end subroutine private_alloc_obsdesc  !!}}}

    subroutine private_dealloc_obsdesc(od)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="private_dealloc_obsdesc"
        type(type_obsdesc) :: od
        deallocate(od%xpos)
        deallocate(od%ypos)
        deallocate(od%zpos)
        deallocate(od%error)
        deallocate(od%qcflag)
        deallocate(od%tidx)
        deallocate(od%varname)
        od%n_obs=0
    end subroutine private_dealloc_obsdesc  !!}}}

    subroutine query_lines(filename, nline)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="query_lines"
        character(len=*) :: filename
        integer :: nline, stat, i
        open(199,file=trim(filename),status="OLD",iostat=stat)
        if(stat /= 0) then
            RP_ERROR("NO "//trim(filename)//" FIND!")
        end if
        i=0 
        do
            read(199, *, iostat=stat)
            if(stat /= 0) exit
            i = i+1
        end do
        nline=i
    end subroutine query_lines  !!}}}


end module module_obs
