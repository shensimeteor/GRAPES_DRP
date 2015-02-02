! --------------------------------------------------------------------------
! Description:
!   NOW ASSUME :: xcoord/ycoord must be in uniform-grid, (dx/dy constant)
!                 xcoord/ycoord must both in increasing order
!                 global model (xcoord is cyclic)
!                 only support 1 grid(xcoord/ycoord)
!
!   local2dman_set_domain
!   local2dman_set_xmodes  :: will NOT do row_normalize
!   local2dman_set_ymode   :: ...
!   local2dman_interp_togrid  :: ..
!   local2dman_interp_toscatter :: ..
!   local2dman_rownormalize_xy :: do row_normalize for xmodes/ymode
!   local2dman_input
!   local2dman_output
!   local2dman_dealloc
!
! History:
!
!   2014-08-29:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
module module_local2dman
use utility
use module_types
use module_local
implicit none

    ! variables

contains

    ! methods
    subroutine local2dman_set_domain(l2dm, xcoord, ycoord)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="local2dman_set_domain"
        real :: xcoord(:), ycoord(:)
        type(type_local2dman) :: l2dm
        if(l2dm%has_setdomain) then
            RP_ERROR("l2dm has set domained already!")
        end if
        l2dm%nx=size(xcoord)
        l2dm%ny=size(ycoord)
        allocate(l2dm%xcoord(l2dm%nx))
        l2dm%xcoord=xcoord
        allocate(l2dm%ycoord(l2dm%ny))
        l2dm%ycoord=ycoord
        allocate(l2dm%xlm_idx(l2dm%ny))
        l2dm%has_setdomain=.true.
    end subroutine local2dman_set_domain  !!}}}

    !!set x modes ; if xl0(i)<0, means the xl0(i)=infinite, i.e.: the pole; x muse be cycle
    !!xlm_idx not surrounded by yidx_bgn(:):yidx_end(:) will be interpolated by xlms
    !!is_xlm_ysym : if true, only to give half yidx_bgn/yidx_end/xl0, the left will be recovered by symmetric 
    subroutine local2dman_set_xmodes(l2dm, yidx_bgn, yidx_end, xl0, n_mode, isxcycle)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="local2dman_set_xmodes"
        type(type_local2dman) :: l2dm
        real :: xl0(:), dx, w1, w2, w
        integer :: yidx_bgn(:), yidx_end(:)
        integer :: nxlm,i, nx, n_mode, nband, j, cnt
        logical :: isxcycle
        if(.not. l2dm%has_setdomain) then
            RP_ERROR("l2dm has NOT set domain yet!")
        end if
        if( .not. all((/size(yidx_bgn),size(yidx_end)/)==size(xl0))) then
            RP_ERROR("input parameter size conflicts")
        end if
        nxlm=l2dm%ny-sum(yidx_end-yidx_bgn+1)
        if(nxlm < 0 ) then
            RP_ERROR("yidx_bgn/yidx_end have OVERLAY numbers")
        end if
        if(.not. isxcycle) then
            RP_ERROR("NOW only support cyclic xcoord")
        end if
        nband=size(yidx_bgn)
        nxlm=nxlm+nband
        allocate(l2dm%xlms(nxlm))
        nx=l2dm%nx
        dx=l2dm%xcoord(2)-l2dm%xcoord(1)
        l2dm%xlm_idx=-1
        do i=1, nband
            l2dm%xlm_idx(yidx_bgn(i):yidx_end(i))=i
            if(xl0(i)<0) then
                call localmodes_generate_Fourier_pole(l2dm%xlms(i), nx, n_mode)
            else
                call localmodes_generate_Fourier(l2dm%xlms(i), nx, dx, xl0(i), 1, n_mode)
            end if
        end do
        cnt=nband
        do i=1, nband-1
            do j=yidx_end(i)+1, yidx_bgn(i+1)-1
                cnt=cnt+1
                l2dm%xlm_idx(j)=cnt
                w1=l2dm%ycoord(yidx_bgn(i+1))-l2dm%ycoord(j)
                w2=l2dm%ycoord(j)-l2dm%ycoord(yidx_end(i))
                w=l2dm%ycoord(yidx_bgn(i+1))-l2dm%ycoord(yidx_end(i))
                w1=w1/w
                w2=w2/w
                call localmodes_linear_combine(l2dm%xlms(i), w1, l2dm%xlms(i+1), w2, l2dm%xlms(cnt))
            end do
        end do
        if(any(l2dm%xlm_idx <= 0)) then
            RP_ERROR("(any(l2dm%xlm_idx <= 0)")
        END IF
        l2dm%n_xlm=cnt
        l2dm%has_setx=.true.
    end subroutine local2dman_set_xmodes  !!}}}

    !! dy_ratio: size(ny)-1
    subroutine local2dman_set_ymode(l2dm, dy_ratio, yl0, op_nmodes, op_eig2sumratio)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="local2dman_set_ymode"
        type(type_local2dman) :: l2dm
        real :: dy_ratio(:), yl0, ycoord(size(dy_ratio)+1)
        integer :: i,j, ny
        integer,optional :: op_nmodes
        real,optional :: op_eig2sumratio
        if(.not. l2dm%has_setdomain) then
            RP_ERROR("l2dm has NOT set domain yet!")
        end if
        if( size(dy_ratio) /= l2dm%ny-1) then
            RP_ERROR("size(dy_ratio) /= l2dm%ny-1")
        end if
        ny=l2dm%ny
        ycoord(1)=0
        do i=2, ny
            ycoord(i)=ycoord(i-1)+dy_ratio(i-1)*(l2dm%ycoord(i)-l2dm%ycoord(i-1))
        end do
        if(present(op_nmodes) .and. present(op_eig2sumratio)) then
            call localmodes_generate_decompC(l2dm%ylm, ycoord, yl0, 1, .false., 0.,&
                      op_nmodes=op_nmodes, op_eig2sumratio=op_eig2sumratio)
        else if(present(op_nmodes)) then
            call localmodes_generate_decompC(l2dm%ylm, ycoord, yl0, 1, .false., 0.,&
                      op_nmodes=op_nmodes)
        else if(present(op_eig2sumratio)) then
            call localmodes_generate_decompC(l2dm%ylm, ycoord, yl0, 1, .false., 0.,&
                      op_eig2sumratio=op_eig2sumratio)
        else
            call localmodes_generate_decompC(l2dm%ylm, ycoord, yl0, 1, .false., 0.)
        end if
        l2dm%has_sety=.true.
    end subroutine local2dman_set_ymode  !!}}}

    subroutine local2dman_dealloc(l2dm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="local2dman_dealloc"
        type(type_local2dman) :: l2dm
        deallocate(l2dm%xcoord, l2dm%ycoord, l2dm%xlm_idx)
        if(allocated(l2dm%xlms)) then
            deallocate(l2dm%xlms)
        end if
        if(l2dm%ylm%m <= 0) then
            call localmodes_delete(l2dm%ylm)
        end if
        l2dm%n_xlm=0
        l2dm%nx=0
        l2dm%ny=0
        l2dm%has_setx=.false.
        l2dm%has_sety=.false.
        l2dm%has_setdomain=.false.
    end subroutine local2dman_dealloc  !!}}}

    !!use l2dm1 to generate l2dm2(if already set, dealloc)
    subroutine local2dman_interp_togrid(l2dm1, l2dm2, xcoord, ycoord, isxcycle)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="local2dman_interp_togrid"
        type(type_local2dman) :: l2dm1, l2dm2
        type(type_localmodes) :: tmp_lm
        real :: xcoord(:), ycoord(:), weight(size(ycoord),2)
        integer :: neibor(size(ycoord),2)
        logical :: isxcycle
        logical,allocatable :: has_xlm(:)
        integer :: i,j, cnt_xlms
        if(.not. (l2dm1%has_setx .and. l2dm1%has_sety)) then
            RP_ERROR("l2dm NOT setx .or. sety")
        end if
        if(l2dm2%has_setdomain) then
            call local2dman_dealloc(l2dm2)
        end if
        call local2dman_set_domain(l2dm2, xcoord, ycoord)
        call interp1d_get_neibor_weight(l2dm1%ycoord, l2dm2%ycoord, neibor, weight)
        cnt_xlms=0
        allocate(has_xlm(l2dm1%n_xlm))
        has_xlm=.false.
        do i=1, l2dm2%ny
           !! write(*,'(5I)') i,neibor(i,1),neibor(i,2),l2dm1%xlm_idx(neibor(i,1)),l2dm1%xlm_idx(neibor(i,2))
            if(l2dm1%xlm_idx(neibor(i,1))==l2dm1%xlm_idx(neibor(i,2))) then
                if(.not. has_xlm(l2dm1%xlm_idx(neibor(i,1)))) then
                    cnt_xlms=cnt_xlms+1
                   !! print*, "interp:", l2dm2%ycoord(i)
                    has_xlm(l2dm1%xlm_idx(neibor(i,1)))=.true.
                end if
                l2dm2%xlm_idx(i)=cnt_xlms
            else
                cnt_xlms=cnt_xlms+1
               !! print*, "interp:", l2dm2%ycoord(i)
                l2dm2%xlm_idx(i)=cnt_xlms
            end if
        end do
        l2dm2%n_xlm=cnt_xlms
        allocate(l2dm2%xlms(cnt_xlms))
        cnt_xlms=0
        has_xlm=.false.
        do i=1, l2dm2%ny
            if(l2dm1%xlm_idx(neibor(i,1))==l2dm1%xlm_idx(neibor(i,2))) then
                if(.not. has_xlm(l2dm1%xlm_idx(neibor(i,1)))) then
                    cnt_xlms=cnt_xlms+1
                    has_xlm(l2dm1%xlm_idx(neibor(i,1)))=.true.
                    call localmodes_interp_linear(l2dm1%xcoord, l2dm1%xlms(l2dm1%xlm_idx(neibor(i,1))), &
                                                  l2dm2%xcoord, l2dm2%xlms(cnt_xlms), isxcycle)
                end if
            else
                cnt_xlms=cnt_xlms+1
                call localmodes_linear_combine(l2dm1%xlms(l2dm1%xlm_idx(neibor(i,1))),  weight(i,1), &
                                               l2dm1%xlms(l2dm1%xlm_idx(neibor(i,2))),  weight(i,2),  tmp_lm)
                call localmodes_interp_linear(l2dm1%xcoord, tmp_lm, l2dm2%xcoord, l2dm2%xlms(cnt_xlms), isxcycle)
            end if
        end do
        call localmodes_interp_linear(l2dm1%ycoord, l2dm1%ylm, ycoord, l2dm2%ylm, .false.)
        l2dm2%has_setdomain=.true.
        l2dm2%has_setx=.true.
        l2dm2%has_sety=.true.
    end subroutine local2dman_interp_togrid  !!}}}

    !!use l2dm to generate xlm/ylm
    subroutine local2dman_interp_toscatter(l2dm, xpos, ypos, xlm, ylm, isxcycle)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="local2dman_interp_toscatter"
        type(type_local2dman) :: l2dm
        real :: xpos(:), ypos(:), dy
        type(type_localmodes) :: xlm, ylm, tmp_lm
        integer :: nsct, i, iy
        logical :: isxcycle
        if(.not. (l2dm%has_setx .and. l2dm%has_sety)) then
            RP_ERROR("l2dm NOT setx .or. sety")
        end if
        nsct=size(xpos)
        if(nsct /= size(ypos)) then
            RP_ERROR("xpos & ypos size conflicts!")
        end if
        dy=l2dm%ycoord(2)-l2dm%ycoord(1)
        call localmodes_init(xlm, nsct, l2dm%xlms(1)%m)
        do i=1, nsct
            iy=nint((ypos(i)-l2dm%ycoord(1))/dy+1) 
            call localmodes_interp_linear(l2dm%xcoord, l2dm%xlms(l2dm%xlm_idx(iy)), &
                                          (/xpos(i)/), tmp_lm, isxcycle)
            xlm%modes(i,:)=tmp_lm%modes(1,:)
            call localmodes_delete(tmp_lm)
        end do
        call localmodes_interp_linear(l2dm%ycoord, l2dm%ylm, ypos, ylm, .false.)
    end subroutine local2dman_interp_toscatter  !!}}}

    subroutine local2dman_rownormalize_xy(l2dm)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="local2dman_rownormalize_xy"
        type(type_local2dman) :: l2dm
        integer :: i
        if(.not. (l2dm%has_setx .and. l2dm%has_sety)) then
            RP_ERROR("l2dm NOT setx .or. sety")
        end if
        do i=1, l2dm%n_xlm
            call localmodes_rownormalize(l2dm%xlms(i))
        end do
        call localmodes_rownormalize(l2dm%ylm)
    end subroutine local2dman_rownormalize_xy  !!}}}

    !! interp x1(n1)->x2(n2), neibor(n2,2), weight(n2,2)
    !! x1/x2 must increasing order; must NOT cyclic
    subroutine interp1d_get_neibor_weight(x1, x2, neibor, weight)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="interp1d_get_neibor_weight"
        real :: x1(:), x2(:), weight(:,:)
        integer :: neibor(:,:)
        integer :: n1, n2, i, j, j1, j2
        n1=size(x1)
        n2=size(x2)
        j1=1
        j2=n1-1
        do i=1,n2
            do j=j1,j2
                if(x2(i) > x1(j) .and. x2(i) < x1(j+1)) then
                    neibor(i,1)=j
                    neibor(i,2)=j+1
                    weight(i,1)=(x1(j+1)-x2(i))/(x1(j+1)-x1(j))
                    weight(i,2)=(x2(i)-x1(j))/(x1(j+1)-x1(j))
                    j1=j
                    exit
                 else if(x2(i) == x1(j)) then
                     neibor(i,1)=j
                     neibor(i,2)=j
                     weight(i,1)=1
                     weight(i,2)=0
                     j1=j
                     exit
                 else if(x2(i)==x1(n1)) then
                     neibor(i,1)=n1
                     neibor(i,2)=n1
                     weight(i,1)=1
                     weight(i,2)=0
                     exit
                 end if
             end do
         end do
    end subroutine interp1d_get_neibor_weight  !!}}}

    subroutine local2dman_output(l2dm, outdir)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="local2dman_output"
        type(type_local2dman) :: l2dm
        integer :: ifile
        character(len=*) :: outdir
        character(len=100) :: prefix, strfmt, filename
        if(.not. (l2dm%has_setx .and. l2dm%has_sety)) then
            RP_ERROR("l2dm NOT setx .or. sety")
        end if
        call system("test -d "//trim(outdir)//" || mkdir "//trim(outdir))
        open(11, file=trim(outdir)//"/local2dman.head", status='replace')
        write(11,"(3(I,','))") l2dm%nx, l2dm%ny, l2dm%n_xlm
        write(strfmt, *) '(',l2dm%nx,"(F,','))"
        write(11,trim(strfmt)) l2dm%xcoord
        write(strfmt, *) '(',l2dm%ny,"(F,','))"
        write(11,trim(strfmt)) l2dm%ycoord
        write(strfmt, *) '(',l2dm%ny,"(I,','))"
        write(11,trim(strfmt)) l2dm%xlm_idx
        close(11)
        prefix=trim(outdir)//"/xmodes"
        do ifile=1, l2dm%n_xlm
            write(filename, '(A, I3.3, A)') trim(prefix), ifile, '.csv'
           !! print*,ifile, l2dm%xlms(ifile)%lenn, l2dm%xlms(ifile)%m
            call localmodes_write_csv(trim(filename), l2dm%xcoord, l2dm%xlms(ifile))
        end do
        filename=trim(outdir)//"/ymode.csv"
        call localmodes_write_csv(trim(filename), l2dm%ycoord, l2dm%ylm)
    end subroutine local2dman_output  !!}}}

    subroutine local2dman_input(l2dm, inputdir)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="local2dman_input"
        type(type_local2dman) :: l2dm
        character(len=*) :: inputdir
        character(len=100) :: filename
        character(len=3) :: str
        integer :: nx ,ny, n_xlm, i
        real,allocatable :: xcoord(:), ycoord(:)
        real,allocatable :: tmp(:)
        if(.not. all((/l2dm%has_setdomain, l2dm%has_setx, l2dm%has_sety/) == .false.)) then
            call local2dman_dealloc(l2dm)
        end if
        open(11, file=trim(inputdir)//"/local2dman.head", status="old")
        read(11, *) nx, ny, n_xlm
        allocate(xcoord(nx))
        allocate(ycoord(ny))
        read(11, *) xcoord
        read(11, *) ycoord
        call local2dman_set_domain(l2dm, xcoord, ycoord)
        l2dm%n_xlm=n_xlm
        allocate(l2dm%xlms(n_xlm))
        read(11, *) l2dm%xlm_idx
        do i=1, n_xlm
            write(str, "(I3.3)") i
            filename=trim(inputdir)//"/xmodes"//str//".csv"
            call localmodes_read_csv(trim(filename), tmp, l2dm%xlms(i))
            deallocate(tmp)
        end do
        filename=trim(inputdir)//"/ymode.csv"
        call localmodes_read_csv(trim(filename), tmp, l2dm%ylm)
        deallocate(tmp, xcoord, ycoord)
        l2dm%has_setdomain=.true.
        l2dm%has_setx=.true.
        l2dm%has_sety=.true.
    end subroutine local2dman_input  !!}}}

    !!output prefix.head, prefix.xmode2d, prefix.ymode, the prefix.*mode* is binary
    subroutine local2dman_output_binary(l2dm, out_prefix)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="local2dman_output_binary"
        type(type_local2dman) :: l2dm 
        character(len=*) :: out_prefix
        character(len=100) :: strfmt
        integer :: i,j 
        open(11, file=trim(out_prefix)//".l2dmhead", status='replace')
        write(11,"(5(I,','))") l2dm%nx, l2dm%ny, l2dm%n_xlm, l2dm%xlms(1)%m, l2dm%ylm%m
        write(strfmt, *) '(',l2dm%nx,"(F,','))"
        write(11,trim(strfmt)) l2dm%xcoord
        write(strfmt, *) '(',l2dm%ny,"(F,','))"
        write(11,trim(strfmt)) l2dm%ycoord
        write(strfmt, *) '(',l2dm%ny,"(I,','))"
        write(11,trim(strfmt)) l2dm%xlm_idx
        close(11)
        open(21, file=trim(out_prefix)//".xmode2d", form="unformatted", access="direct", recl=l2dm%xlms(1)%lenn)
        do i=1, l2dm%ny
            do j=1, l2dm%xlms(i)%m
                write(21, rec=(i-1)*l2dm%xlms(1)%m+j) l2dm%xlms(l2dm%xlm_idx(i))%modes(:,j)
            end do
        end do
        close(21)
        open(31, file=trim(out_prefix)//".ymode", form="unformatted", access="direct", recl=l2dm%ylm%lenn)
        do i=1, l2dm%ylm%m
            write(31, rec=i) l2dm%ylm%modes(:,i)
        end do
        close(31)
    end subroutine local2dman_output_binary  !!}}}

end module module_local2dman
