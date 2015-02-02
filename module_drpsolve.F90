! --------------------------------------------------------------------------
! Description:
!
!   drpsolve, should be independent of the real model/obs
!   but can be situation dependent, i.e. what does pxfile/pyfile represent may differs the code
!
! History:
!
!   2014-08-07:
!
!       [shensi]: First creation.
! --------------------------------------------------------------------------
#include "inc_common.fpp"
module module_drpsolve
use utility
use lapack95
use f95_precision
use module_types
use module_ensemble
use module_modelstate
use module_obs
implicit none

contains
    !! depressed, use config_get_namelist
    subroutine drp_setup(dc)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_setup"
        type(type_drpconfig) :: dc
        integer :: i,j
        character(len=4) :: suffix
        character(len=2) :: suffix2
        dc%n_obs_slot=7 
        dc%with_Bterm=.true.
        dc%with_Cterm=.false.
        dc%Balpha_option=3
        dc%do_Py_EOF=.false.
        dc%pertPy_option=2
        dc%inflate_option=0
        dc%inflator=1.0
        dc%do_ETKF=.true.
        dc%max_nt=50
        dc%eps=0.01
        dc%analysis_time=0
        dc%n_mem=20 !!
        dc%is_outerloop=.true.
        dc%is_outerloop_1st=.true. !!
        dc%max_n_obs=2000
        allocate(dc%pxfiles(dc%n_mem))
        allocate(dc%px00files(dc%n_mem))
        allocate(dc%pyfiles(dc%n_obs_slot, dc%n_mem))
        allocate(dc%obfile(dc%n_obs_slot))
        allocate(dc%ybfile(dc%n_obs_slot))
        allocate(dc%gxfiles(dc%n_mem))
        dc%xbfile="xb"
        dc%xgfile="xg"
        dc%xb00file="xb00"
        dc%xg00file="xg00"
        dc%betagfile="betag"
        dc%xafile="xa"
        dc%xa00file="xa00"
        !dc%betafile="beta"
        do i=1, dc%n_mem
            write(suffix, '(A1,I3.3)') 'e',i
            dc%pxfiles(i)="ENS/px."//trim(suffix)
            dc%px00files(i)="ENS/px00."//trim(suffix)
            do j=0, dc%n_obs_slot-1
                write(suffix2, '(I2.2)') j
                dc%pyfiles(j,i)="ENS/py"//trim(suffix2)//"."//trim(suffix)
            end do
            dc%gxfiles(i)="ENS/gx."//trim(suffix)
        end do
        do j=0, dc%n_obs_slot-1
            write(suffix2, '(I2.2)') j
            dc%ybfile(j)="yb"//trim(suffix2)
            dc%obfile(j)="ob"//trim(suffix2)
        end do
        !!diagnosis output settings
        dc%to_output_Bbalpha=.true.
        dc%to_output_jval=.true.
        dc%to_output_jgrad=.true.
    end subroutine drp_setup  !!}}}

    !!here situation: Pxfile is Px full field, Pyfile is Py full field
    !!not fulfill: PyEOF case
    subroutine drp_run(dc, ld)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_run"
        type(type_drpconfig) :: dc
        type(type_localdesc) :: ld
        type(type_modeldesc) :: md
        type(type_obsdesc) :: od
        type(type_ensPx) :: ensPx 
        type(type_ensPy) :: ensPy
        type(type_beta) :: beta, betag
        real,allocatable :: xb(:), xg(:), xb00(:), xg00(:), &
                            balpha(:,:), xa(:), xa00(:), pydcorr(:), xtemp(:)
        real,allocatable :: d_max(:), d(:), yb(:), ybase(:), ob(:), dy_sim(:)
        logical :: to_io_00 !! to read/write **00 file
        logical :: to_read_yb !! if use the averge of Py, no need to read yb
        logical :: to_read_beta, to_write_beta
        real :: lamda, test1, test2
        real, allocatable :: Gx(:,:) 
        integer :: i
       
        !!1. set md, read xb/xg/xb00/xg00
        RP_LOG("DRP read xb/xg/Px ---")
        call modeldesc_init(md)
        call modeldesc_getcoord(md, trim(dc%xbfile))
        allocate(xb(md%n_total))
        allocate(xg(md%n_total))
        call modelstate_input(md, trim(dc%xbfile), xb)
        call modelstate_input(md, trim(dc%xgfile), xg)
        to_io_00=(dc%analysis_time /= 0 .and. dc%is_outerloop)
        if(to_io_00) then
            allocate(xb00(md%n_total))
            allocate(xg00(md%n_total))
            call modelstate_input(md, trim(dc%xb00file), xb00)
            call modelstate_input(md, trim(dc%xg00file), xg00)
        end if
        !!2. read Px, calculate the perturbed Px
        call Px_read(md, dc%pxfiles, ensPx)
        call Px_full2pert(ensPx, xb)
        TRACEX(sum(ensPx%Px))
        TRACEX(sum(abs(ensPx%Px)))
        TIME_CLICK(timer, timer_bgn, timer_last)

        !!3. read od & ob & yb
        RP_LOG("DRP read&normalize d/Py[/yb] ---")
        allocate(d_max(dc%max_n_obs))
        call obs_read_txt(od, d_max, dc%obfile, (/(i-1,i=1,dc%n_obs_slot)/))
        allocate(ob(od%n_obs),d(od%n_obs))
        ob=d_max(1:od%n_obs)
        deallocate(d_max)
        to_read_yb=dc%do_read_yb
        allocate(yb(od%n_obs))
        if(to_read_yb) then  !!yb read from file
            call obs_read_txt(od, yb, dc%ybfile, (/(i-1,i=1,dc%n_obs_slot)/))
        end if
        call obs_group(od)
!!        call test_print_obsgroup(od)
        
        !!4. filter od/d
        !not in Baro Model

        !!5. read Py, calculate the perturbed Py, calculate d
        call Py_read(od, dc%pyfiles, (/(i-1,i=1,dc%n_obs_slot)/), ensPy)
        allocate(ybase(od%n_obs)) !! ptb_Py=Py-ybase
        if( (.not. dc%do_read_yb) .or. (dc%pertPy_option == 2)) then
            !!if pertPy_option==1 and do_read_yb=false, not happens
            ybase=0
            do i=1, ensPy%n_mem
                ybase=ybase+ensPy%Py(:,i)
            end do
            ybase=ybase/ensPy%n_mem
        else
            ybase=yb
        end if
        if(dc%do_read_yb) then !! use yb to cal. d=yo-yb
            d=ob-yb
        else !!do_read_yb=false
            yb=ybase
            d=ob-yb
        end if
        call Py_full2pert(ensPy, ybase)
        TRACEX(sum(ensPy%Py))
        TRACEX(sum(abs(ensPy%Py)))
        allocate(pydcorr(ensPy%n_mem))
        call drp_get_Pydcorr(dc, ensPy, d, pydcorr, .false.)

        !!6. normalize d/Py
        call obs_normalize(od, d)
        call obs_normalize(od, ob)
        call obs_normalize(od, yb)
        do i=1, ensPy%n_mem
            call obs_normalize(od, ensPy%Py(:,i))
        end do
        TIME_CLICK(timer, timer_bgn, timer_last)

  
        !!x. PyEOF ? , Px->Px, Py->Py
        
        !!7. prepare Balpha & Px->Px*balpha, Py->Py*balpha
        RP_LOG("DRP prepare Balpaha & transform Px/Py")
        allocate(balpha(ensPx%n_mem, ensPx%n_mem))
        call drp_get_balpha(dc, balpha, ensPx%n_mem)
        ensPx%Px=matmul(ensPx%Px, balpha)
        ensPy%Py=matmul(ensPy%Py, balpha)
        TIME_CLICK(timer, timer_bgn, timer_last)


        !!8. inflate
        if(dc%inflate_option == 1) then
            RP_LOG("DRP inflate : use inflator = "//tostring(fval=dc%inflator))
            call drp_update_useinflator(balpha, ensPx, ensPy, dc%inflator)
            TIME_CLICK(timer, timer_bgn, timer_last)
        else if(dc%inflate_option == 2) then
            TRACEX(dc%low_limit_inflator)
            TRACEX(dc%up_limit_inflator)
            call drp_adaptive_inflator(dc, od, d, ensPy, lamda)
            RP_LOG("DRP inflate : adaptive inflator is ="//tostring(fval=lamda))
            call drp_update_useinflator(balpha, ensPx, ensPy, lamda)
            if(dc%to_output_inflator) then
                open(63,file="diagnosis_inflator.csv")
                write(63,*) lamda
                close(63)
            end if
            TIME_CLICK(timer, timer_bgn, timer_last)
        end if
        !!9. localize & define Control Variable beta
        RP_LOG("DRP localize Px/Py")
        call Px_localize(ensPx, ld)
        call Py_localize(ensPx, ensPy, od, ld)

        TRACEX(ensPx%n_total_mem)
        call drp_init_beta(ensPx, beta)
        call drp_init_beta(ensPx, betag)
        TIME_CLICK(timer, timer_bgn, timer_last)

        !!10. minimise
        RP_LOG("DRP begin minimize")
        to_read_beta=(dc%is_outerloop .and. (.not. dc%is_outerloop_1st))
        if(to_read_beta) then
            call drp_read_beta(trim(dc%betagfile), betag)
            RP_LOG("DRP read betag")
        end if
        if(dc%minisolver == 0) then
            call drp_solve_direct(dc, od, ensPy, betag, beta, d)
        else
            call drp_minimize_cg(dc, od, ensPy, betag, beta, d) 
        end if
        TIME_CLICK(timer, timer_bgn, timer_last)

        !!11. update betag & output betag
        betag%array=beta%array+betag%array
        to_write_beta=dc%is_outerloop
        if(to_write_beta) then
            RP_LOG("DRP output betag")
            call drp_write_beta(trim(dc%betagfile), betag)
            TIME_CLICK(timer, timer_bgn, timer_last)
        end if

        !!11.x output ob's
        if(dc%to_output_ob) then
            allocate(dy_sim(od%n_obs))
            call Py_multiply_vector(od, ensPy, beta%array, dy_sim) 
            call write_float_matrix2( reshape((/ob, yb, d, dy_sim, dy_sim-d/), (/od%n_obs,5/)), &
                 fname="diagnosis_nmobs.csv",element_fmat="F14.4,','")
            deallocate(dy_sim)
        end if

        !!12. get xa & output xa
        RP_LOG("DRP output xa")
        allocate(xa(md%n_total))
        call Px_multiply_vector(md, ensPx, beta%array, xa)
        TRACEX(maxval(xa))
        xa=xg+xa
        call system("cp -f "//trim(dc%xbfile)//" "//trim(dc%xafile))
        call modelstate_output(md, trim(dc%xafile), xa)
        TIME_CLICK(timer, timer_bgn, timer_last)
     !   !!test
     !   xtemp=reshape(matmul(ensPx%Px, reshape(ensPy%Py(1,:), (/ensPy%n_mem,1/))), (/ensPx%n_state/))/(ensPx%n_mem-1)
     !   TRACEX(maxval(xtemp))
     !   TRACEX(maxloc(xtemp))
     !   call Py_get_row(od, ensPy, 1, beta%array) 
     !   call Px_multiply_vector(md, ensPx, beta%array, xtemp)
     !   TRACEX(maxval(xtemp)/(ensPx%n_mem-1))
     !   TRACEX(maxloc(xtemp))
    

        !!13. get & output Gx
        if(dc%do_ETKF) then
            RP_LOG("DRP ETKF : calculate Gx")
            allocate(Gx(ensPx%n_state, ensPx%n_mem))
            call drp_ETKF(ensPx, ensPy, Gx)
            ensPx%Px=Gx
            RP_LOG("DRP ETKF : output Gx")
            do i=1, ensPx%n_mem
                call system("cp -f "//trim(dc%pxfiles(i))//" "//trim(dc%gxfiles(i)))
            end do 
            call Px_write(md, dc%gxfiles, ensPx)
            TIME_CLICK(timer, timer_bgn, timer_last)
        end if

        !!14. for xa00 output 
        if(to_io_00) then
            RP_LOG("DRP xa00: calculate xa00")
            allocate(xa00(md%n_total))
            deallocate(ensPx%Px)
            call Px_read(md, dc%px00files, ensPx)
            ensPx%n_total_mem=ensPy%n_total_mem
            call Px_full2pert(ensPx, xb00)
            ensPx%Px=matmul(ensPx%Px, balpha)
            call Px_multiply_vector(md, ensPx, beta%array, xa00)
            xa00=xg00+xa00
            RP_LOG("DRP xa00: output xa00")
            call system("cp -f "//trim(dc%xb00file)//" "//trim(dc%xa00file))
            call modelstate_output(md, trim(dc%xa00file), xa00)
            TIME_CLICK(timer, timer_bgn, timer_last)
        end if

        !!15. deallocate
        deallocate(xb,xg,xa,balpha,d,yb)
        call Px_delete(ensPx)
        call Py_delete(ensPy)
        if(to_io_00) then
            deallocate(xb00,xg00,xa00)
        end if
    end subroutine drp_run  !!}}}

    subroutine drp_minimize_cg(dc, od, ensPy, beta_g, beta, d)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_minimize_cg"
        type(type_drpconfig) :: dc
        type(type_obsdesc) :: od
        type(type_ensPy) :: ensPy
        type(type_beta) :: beta_g, beta
        real :: d(:)
        real,allocatable :: grad(:), grads(:,:), dk(:), adk(:)
        integer :: ny, m, i, t
        real :: jval, jvals(3), grad_nm, grads_nm(3), grad_nm2
        real :: stop_grad_nm, alpha, gammar, grad_nm2_last
        real,allocatable :: jo_val(:)
        logical :: stop_by_grad
        character(len=100) :: string, strfmt
        m=ensPy%n_total_mem
        allocate(grad(m))
        allocate(grads(m,3))
        ny=ensPy%n_state
        if(dc%max_nt < 1) then
            return
        end if
        if(dc%to_output_jval) then
            open(58, file="diagnosis_jval.csv", status='replace')
            if(od%has_group .and. od%has_group_set) then
                allocate(jo_val(od%n_group))
                string="; JoGroups: "
                do i=1, od%n_group
                    string=trim(string)//trim(od%group(i)%group_name)//","
                end do
                write(58, *) 'Iterations, J, Jb, Jo, Jc'//trim(string)
                call drp_get_J(dc, od, beta_g, ensPy, d, beta, jval, jvals, jo_val)
                write(strfmt, "(A, I3, A)") '(I3, ",",', 4+od%n_group, '(F,","))'
                write(58, trim(strfmt)) 0, jval, jvals, jo_val
            else
                write(58, *) 'Iterations, J, Jb, Jo, Jc'
                call drp_get_J(dc, od, beta_g, ensPy, d, beta, jval, jvals)
                write(58, '(I3, ",", 4(F,","))') 0,jval,jvals, jo_val
            end if
        end if
        call drp_get_GradJ(dc, od, beta_g, ensPy, d, beta, grad, grads)
        grad_nm2=sum(grad**2)
        grad_nm=sqrt(grad_nm2)
        RP_LOG("Minimize_CG begin, ||Grad(J@xg)||="//trim(tostring(fval=grad_nm)))
        TIME_CLICK(timer, timer_bgn, timer_last)
        stop_grad_nm=grad_nm*dc%eps
        if(dc%to_output_jgrad) then
            do i=1, 3
                grads_nm(i)=sqrt(sum(grads(:,i)**2))
            end do
            open(59,file="diagnosis_jgrad.csv", status="replace")
            write(59, *) 'Iterations, ||Grad(J)||, ||Grad(Jb)||, ||Grad(Jo)||, ||Grad(Jc)||'
            write(59, '(I3,",", 4(F,","))') 0,grad_nm, grads_nm
        end if
        allocate(dk(m))
        allocate(adk(m))
        dk=-grad
        stop_by_grad=.false.
        do i=1, dc%max_nt
            !!1. get alpha (the length-ratio of dk) to go to next beta
            call drp_get_Ax(dc, od, ensPy, dk, adk)
!            TRACEX(dk(maxloc(abs(dk))))
!            TRACEX(adk(maxloc(abs(adk))))
            alpha=sum(-grad*dk)/sum(dk*adk)
            !!2. update beta
            beta%array=beta%array+alpha*dk
            !!3. calculate jval & grad & check stop
            if(dc%to_output_jval) then
                if(od%has_group .and. od%has_group_set) then
                    call drp_get_J(dc, od, beta_g, ensPy, d, beta, jval, jvals, op_jogroups=jo_val)
                    write(strfmt, "(A, I3, A)") '(I3, ",",', 4+od%n_group, '(F,","))'
                    write(58, trim(strfmt)) i,jval,jvals, jo_val
                else
                    call drp_get_J(dc, od, beta_g, ensPy, d, beta, jval, jvals)
                    write(58, '(I3,",", 4(F,","))') i,jval,jvals
                end if
            end if
            call drp_get_GradJ(dc, od, beta_g, ensPy, d, beta, grad, grads)
            grad_nm2_last=grad_nm2
            grad_nm2=sum(grad**2)
            grad_nm=sqrt(grad_nm2)
            if(dc%to_output_jgrad) then
                do t=1, 3
                    grads_nm(t)=sqrt(sum(grads(:,t)**2))
                end do
                write(59, '(I3,",", 4(F,","))') i,grad_nm, grads_nm
            end if
            if(grad_nm <= stop_grad_nm) then
                RP_LOG("Minimize_CG finish after "//trim(tostring(ival=i,fmat="(I3)"))// &
                       "iterations with ||Grad(J)||="//trim(tostring(fval=grad_nm)))
                stop_by_grad=.true.
                exit
            end if
            !!4. update dk
            gammar=grad_nm2/grad_nm2_last
            dk=-grad+gammar*dk
            RP_LOG("Minimize_CG after iteration:"//trim(tostring(ival=i, fmat="(I3)")))
            TIME_CLICK(timer, timer_bgn, timer_last)
        end do
        if(.not. stop_by_grad) then
            RP_LOG("Minimize_CG ended after all "//trim(tostring(ival=dc%max_nt,fmat="(I3)")) &
                   //"iterations with ||Grad(J)||="//trim(tostring(fval=grad_nm)))
        end if
        if(dc%to_output_jval) close(58)
        if(dc%to_output_jgrad) close(59)
        deallocate(grad, grads, dk, adk)
    end subroutine drp_minimize_cg  !!}}}

    !!unlike drp_minimize_cg, it solves the linear equations directly(GradJ=0)
    !!only for non-localized solution
    subroutine drp_solve_direct(dc, od, ensPy, beta_g, beta, d)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_solve_direct"
        type(type_drpconfig) :: dc
        type(type_ensPy) :: ensPy
        type(type_beta) :: beta_g, beta
        type(type_obsdesc) :: od
        real :: d(:)
        real,allocatable :: A(:,:), B(:,:)
        integer :: n,i,stat
        real :: grad_nm2, grad_nm, grads_nm(3)
        real,allocatable :: grad(:), grads(:,:), jo_val(:)
        real :: jval, jvals(3)
        character(len=200) :: string, strfmt
        if(dc%with_Cterm) then
            RP_WARNING("with_Cterm = true not supported now!")
        end if
        RP_LOG("DRP solve directly")
        n=ensPy%n_total_mem
        if(dc%to_output_jval) then
            open(58, file="diagnosis_jval.csv", status='replace')
            if(od%has_group .and. od%has_group_set) then
                allocate(jo_val(od%n_group))
                string="; JoGroups: "
                do i=1, od%n_group
                    string=trim(string)//trim(od%group(i)%group_name)//","
                end do
                write(58, *) 'state, J, Jb, Jo, Jc'//trim(string)
                call drp_get_J(dc, od, beta_g, ensPy, d, beta, jval, jvals, op_jogroups=jo_val)
                write(strfmt, "(A, I3, A)") '(A, ",",', 4+od%n_group, '(F,","))'
                write(58, trim(strfmt)) "guess", jval, jvals, jo_val
            else
                write(58, *) 'state, J, Jb, Jo, Jc'
                call drp_get_J(dc, od, beta_g, ensPy, d, beta, jval, jvals)
                write(58, '(A, ",", 4(F,","))') "guess",jval,jvals
            end if
        end if
        if(dc%to_output_jgrad) then
            allocate(grad(n))
            allocate(grads(n,3))
            call drp_get_GradJ(dc, od, beta_g, ensPy, d, beta, grad, grads)
            grad_nm2=sum(grad**2)
            grad_nm=sqrt(grad_nm2)
            do i=1, 3
                grads_nm(i)=sqrt(sum(grads(:,i)**2))
            end do
            open(59,file="diagnosis_jgrad.csv", status="replace")
            write(59, *) 'State, ||Grad(J)||, ||Grad(Jb)||, ||Grad(Jo)||, ||Grad(Jc)||'
            write(59, '(A,",", 4(F,","))') "guess",grad_nm, grads_nm
        end if
        allocate(A(n,n), B(n,1))
        A=matmul(transpose(ensPy%Py), ensPy%Py)
        do i=1,n
            A(i,i)=A(i,i)+1
        end do
        B=reshape(matmul(transpose(ensPy%Py),d)-beta_g%array, (/n,1/))
        call posv(A,B,'L',stat)
        LPK_CHECK(stat)
        beta%array=B(:,1)
        if(dc%to_output_jval) then
            if(od%has_group .and. od%has_group_set) then
                call drp_get_J(dc, od, beta_g, ensPy, d, beta, jval, jvals, jo_val)
                write(strfmt, "(A, I3, A)") '(A, ",",', 4+od%n_group, '(F,","))'
                write(58, trim(strfmt)) "analys", jval, jvals, jo_val
            else
                call drp_get_J(dc, od, beta_g, ensPy, d, beta, jval, jvals)
                write(58, '(A, ",", 4(F,","))') "analys",jval,jvals
            end if
            close(58)
        end if
        if(dc%to_output_jgrad) then
            call drp_get_GradJ(dc, od, beta_g, ensPy, d, beta, grad, grads)
            grad_nm2=sum(grad**2)
            grad_nm=sqrt(grad_nm2)
            do i=1, 3
                grads_nm(i)=sqrt(sum(grads(:,i)**2))
            end do
            write(59, '(A,",", 4(F,","))') "analys",grad_nm, grads_nm
            close(59)
            deallocate(grad,grads)
        end if
        deallocate(A,B)
    end subroutine drp_solve_direct  !!}}}

    !! J(beta)=0.5*(beta_g+beta)'*(beta_g+beta)+0.5*[ensPy*beta-d]'*[ensPy*beta-d]+0.5*Jc
    !! jval is total j cost function value; jvals(3)
    subroutine drp_get_J(dc, od, beta_g, ensPy, d, beta, jval, jvals, op_jogroups)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_get_J"
        type(type_obsdesc) :: od
        type(type_drpconfig) :: dc
        type(type_beta) :: beta_g, beta
        type(type_ensPy) :: ensPy
        real :: d(:)
        real :: jo, jb, jc, jval, jvals(3)
        real, optional :: op_jogroups(:) 
        real :: y_sim(size(d))
        integer :: i
        if(dc%with_Bterm) then
            jb=0.5*sum((beta_g%array+beta%array)**2)
        else
            jb=0
        end if
        print*,"test: before call Py_multiply_vector"
        call Py_multiply_vector(od, ensPy, beta%array, y_sim) 
        TIME_CLICK(timer, timer_bgn, timer_last)
        jo=0.5*sum((y_sim-d)**2)
        if(dc%with_Cterm) then
            RP_WARNING("with_Cterm = true not supported now!")
            jc=0
        else
            jc=0
        end if
        jval=jb+jo+jc
        jvals=(/jb,jo,jc/)
        if(present(op_jogroups)) then
            if(size(op_jogroups) < od%n_group) then
                RP_ERROR("para op_jogroups size too small")
            end if
            do i=1, od%n_group
                op_jogroups(i)=0.5*sum((y_sim(od%group(i)%idx)-d(od%group(i)%idx))**2)
            end do
        end if
    end subroutine drp_get_J  !!}}}

    !!GradJ(beta)=(beta_g+beta) + ensPy'*[ensPy*beta-d] + GradJc(beta)
    !!jgrad is gradient(J), jgrads is gradient(Jb,Jo,Jc)
    subroutine drp_get_GradJ(dc, od, beta_g, ensPy, d, beta, jgrad, jgrads)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_get_GradJ"
        type(type_obsdesc) :: od
        type(type_drpconfig) :: dc
        type(type_beta) :: beta_g, beta
        type(type_ensPy) :: ensPy
        real :: d(:)
        real :: jgrad(:), jgrads(:,:)
        real :: jbgrad(size(jgrad)), jograd(size(jgrad)), jcgrad(size(jgrad))
        real :: y_sim(size(d))
        integer :: i
        if(dc%with_Bterm) then 
            jbgrad=beta_g%array+beta%array
        else
            jbgrad=0
        end if
        jgrads(:,1)=jbgrad
        print*,"test: before call Py_multiply_vector"
        call Py_multiply_vector(od, ensPy, beta%array, y_sim) 
        TIME_CLICK(timer, timer_bgn, timer_last)
        y_sim=y_sim-d
        print*,"test: before call Py_transpose_multiply_vector"
        call Py_transpose_multiply_vector(od, ensPy, y_sim, jograd)
        TIME_CLICK(timer, timer_bgn, timer_last)
        jgrads(:,2)=jograd
        if(dc%with_Cterm) then
            RP_WARNING("with_Cterm = true not supported now!")
            jcgrad=0
        else
            jcgrad=0
        end if
        jgrads(:,3)=jcgrad
        jgrad=jbgrad+jograd+jcgrad
    end subroutine drp_get_GradJ  !!}}}

    !!Ax(vct)=vct+ensPy'*ensPy*vct + Jc_AX ?
    subroutine drp_get_Ax(dc, od, ensPy, vct, ax)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_get_Ax"
        type(type_drpconfig) :: dc
        type(type_obsdesc) :: od
        type(type_ensPy) :: ensPy
        real :: vct(:), ax(:)
        real :: axb(size(ax)), axo(size(ax)), axc(size(ax))
        real,allocatable :: y_sim(:)
        if(dc%with_Bterm) then 
            axb=vct
        else
            axb=0
        end if
        allocate(y_sim(od%n_obs))
        print*,"test: before call Py_multiply_vector"
        call Py_multiply_vector(od, ensPy, vct, y_sim) 
        TIME_CLICK(timer, timer_bgn, timer_last)
!        TRACEX(y_sim(maxloc(abs(y_sim))))
        print*,"test: before call Py_transpose_multiply_vector"
        call Py_transpose_multiply_vector(od, ensPy, y_sim, axo)
        TIME_CLICK(timer, timer_bgn, timer_last)
!        TRACEX(axo(maxloc(abs(axo))))
        if(dc%with_Cterm) then
            RP_WARNING("with_Cterm = true not supported now!")
            axc=0
        else
            axc=0
        end if
        ax=axb+axo+axc
        deallocate(y_sim)
    end subroutine drp_get_Ax  !!}}}

    !! not sure right or wrong if use PyEOF (check Lorenz96-DRP-4DVar's code)
    !! now assume balpha is square matrix ( though no need to be a square)
    subroutine drp_get_balpha(dc, balpha, m)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_get_balpha"
        type(type_drpconfig) :: dc
        real :: balpha(:,:)
        real,allocatable :: balpha_2(:,:)
        integer :: m, n, i
        n=size(balpha,1)
        if(dc%Balpha_option == 1) then
            balpha=-1./m
            do i=1, n
                balpha(i,i)=1+balpha(i,i)
            end do
            balpha=balpha/sqrt(m-1.)
        else if(dc%Balpha_option == 2) then
            balpha=-1./(m+1)
            do i=1, n
                balpha(i,i)=1+balpha(i,i)
            end do
            balpha=balpha/sqrt(m*1.0)
        else if(dc%Balpha_option == 3) then
            balpha=0
            do i=1, n
                balpha(i,i)=1.0/sqrt(m-1.0)
            end do
        else
            RP_ERROR("Balpha_option  Not Recognized")
        end if
        if(dc%to_output_Bbalpha) then
            call write_float_matrix2(balpha,fname="diagnosis_balpha.csv",element_fmat="F14.4,','")
            allocate(balpha_2(size(balpha,1),size(balpha,1)))
            balpha_2=matmul(balpha, transpose(balpha))
            call write_float_matrix2(balpha_2,fname="diagnosis_Balpha.csv",element_fmat="F14.4,','")
            deallocate(balpha_2)
        end if
    end subroutine drp_get_balpha  !!}}}

    !!get perturbed analysis ensemble Gx
    subroutine drp_ETKF(ensPx, ensPy, Gx )  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_ETKF"
        type(type_ensPx) :: ensPx
        type(type_ensPy) :: ensPy
        real :: Gx(:,:)
        real,allocatable :: Py2(:,:)!!Py'*Py=Py2
        real,allocatable :: eigvalue(:), TransMtrx(:,:)
        integer :: nx, n_mem, stat, i
        nx=size(Gx,1)
        n_mem=size(Gx,2)
        if(nx /= ensPx%n_state) then
            RP_ERROR("size(Gx,1) error!")
        end if
        if(n_mem /= ensPx%n_mem) then
            RP_ERROR("size(Gx,2) error!")
        end if
        allocate(Py2(n_mem, n_mem), eigvalue(n_mem), TransMtrx(n_mem,n_mem))
        Py2=matmul(transpose(ensPy%Py),ensPy%Py) !!Py2=Py'*Py 
        call syevd(Py2, eigvalue, jobz='V',uplo='L',info=stat)  !!Py2=EV of (Py'*Py)
        where(eigvalue < 0) 
            eigvalue=0
        end where
        eigvalue=1/sqrt(eigvalue+1)
        TRACEX(eigvalue)
        !!TransMtrx = EV * sqrt(eigvalue + I) * EV'
        TransMtrx=0
        do i=1, n_mem
            TransMtrx=TransMtrx+eigvalue(i)*matmul(Py2(:,(/i/)),transpose(Py2(:,(/i/))))
        end do
        !!Gx=sqrt(n_mem-1)*Px*TransMtrx
        Gx=matmul(ensPx%Px, TransMtrx)
        Gx=Gx*sqrt(n_mem-1.0)
    end subroutine drp_ETKF  !!}}}

    !!get inflator: lamda, ensPy must have multiply balpha already 
    subroutine drp_adaptive_inflator(dc, od, d, ensPy, lamda)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_adaptive_inflator"
        type(type_obsdesc) :: od
        type(type_drpconfig) :: dc
        real :: d(:), lamda, d2, tr_py2
        type(type_ensPy) :: ensPy
        integer ::  na, my, i
        my=od%n_obs
        d2=sum(d*d) !!d2=d'*d
        tr_py2=0
        do i=1,my
            tr_py2=tr_py2+sum(ensPy%Py(i,:)**2)
        end do  !!tr_py2=tr(Py*Py')
        lamda=(d2-my)/tr_py2
        TRACE(lamda)
        if(lamda<dc%low_limit_inflator) then
            lamda=dc%low_limit_inflator
        else if(dc%up_limit_inflator > 0 .and. lamda > dc%up_limit_inflator) then
            lamda=dc%up_limit_inflator
        end if
    end subroutine drp_adaptive_inflator  !!}}}
    
    !!Px<=Px*sqrt(lamda), Py<=Py*sqrt(lamda), balpha=balpha*sqrt(lamda)
    subroutine drp_update_useinflator(balpha, ensPx, ensPy, lamda )  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_update_useinflator"
        type(type_ensPx) :: ensPx
        type(type_ensPy) :: ensPy
        real :: balpha(:,:)
        real :: lamda
        ensPx%Px=ensPx%Px*sqrt(lamda)
        ensPy%Py=ensPy%Py*sqrt(lamda)
        balpha=balpha*sqrt(lamda)
    end subroutine drp_update_useinflator  !!}}}

    subroutine drp_get_Pydcorr(dc, ensPy, d, corr, minus_bias)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_get_Pydcorr"
        type(type_drpconfig) :: dc
        type(type_ensPy) :: ensPy
        real :: d(:), corr(:)
        logical :: minus_bias
        integer :: m,p,i
        real :: avg_py, avg_d, stdv_py, stdv_d
        m=size(d)
        p=ensPy%n_mem
        if(p>size(corr)) then
            RP_ERROR("p>size(corr), the input array-corr has too small size")
        end if
        if(minus_bias) then
            avg_d=sum(d)/m
        else
            avg_d=0
        end if
        stdv_d=sqrt(sum((d-avg_d)**2))
        do i=1, p
            if(minus_bias) then
                avg_py=sum(ensPy%Py(:,i))/m
            else
                avg_py=0
            end if
            stdv_py=sqrt(sum((ensPy%Py(:,i)-avg_py)**2))
            corr(i)=sum((ensPy%Py(:,i)-avg_py)*(d-avg_d))
            if(stdv_py <0.01) then
                corr(i)=-9.99
            else
                corr(i)=corr(i)/(stdv_d*stdv_py)
            end if
        end do
        if(dc%to_output_pydcorr) then
            call write_float_matrix2(reshape(corr,(/p,1/)),fname="diagnosis_Pydcorr.csv",element_fmat="F10.6,','")
        end if
    end subroutine drp_get_Pydcorr  !!}}}

!!about type_beta
    !! beta must init before read
    subroutine drp_read_beta(beta_file, beta)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_read_beta"
        character(len=*) :: beta_file
        character(len=20) :: strfmt
        type(type_beta) :: beta
        integer :: i,j,k,t
        real,allocatable :: mtrx(:,:,:,:)
        open(41, file=trim(beta_file), status="old")
        read(41, *) i,j,k,t
        if( .not. all((/i,j,k,t/) == (/beta%m1, beta%m2, beta%m3, beta%m4/))) then
            RP_ERROR("size of beta not conform with the beta_file")
        end if
        allocate(mtrx(beta%m4, beta%m3, beta%m2, beta%m1))
        do i=1, beta%m4
            do j=1, beta%m3
                do k=1, beta%m2
                    read(41, *) (mtrx(i,j,k,t),t=1,beta%m1)
                end do
            end do
        end do
        close(41)
        beta%array=reshape(mtrx, (/size(beta%array)/))
        deallocate(mtrx) 
    end subroutine drp_read_beta  !!}}}

    subroutine drp_write_beta(beta_file, beta)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_write_beta"
        character(len=*) :: beta_file
        character(len=20) :: strfmt
        type(type_beta) :: beta
        integer :: i,j,k,t
        real,allocatable :: mtrx(:,:,:,:)
        allocate(mtrx(beta%m4, beta%m3, beta%m2, beta%m1))
        mtrx=reshape(beta%array,(/beta%m4, beta%m3, beta%m2, beta%m1/))
        open(41, file=trim(beta_file), status="replace")
        write(41, '(4(I,","))') beta%m1, beta%m2, beta%m3, beta%m4
        strfmt=" "
        write(strfmt, "('(',I3,'(G, "",""))')" ) beta%m1
        do i=1, beta%m4
            do j=1, beta%m3
                do k=1, beta%m2
                    write(41, trim(strfmt)) (mtrx(i,j,k,t),t=1,beta%m1)
                end do
            end do
        end do
        close(41)
        deallocate(mtrx)
    end subroutine drp_write_beta  !!}}}

    subroutine drp_init_beta(ensPx, beta)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="drp_init_beta"
        type(type_ensPx) :: ensPx
        type(type_beta) :: beta
        beta%m1=ensPx%n_mem
        beta%m2=ensPx%lman%x3dm(1,1)%m
        beta%m3=ensPx%lman%y2dm(1)%m
        beta%m4=ensPx%lman%zm%m
        allocate(beta%array(beta%m4 * beta%m3 * beta%m2 * beta%m1))
        beta%array=0
    end subroutine drp_init_beta  !!}}}

    subroutine test_print_obsgroup(od)  !!{{{
    implicit none
        character(len=*),parameter :: PROCEDURE_NAME="test_print_obsgroup"
        type(type_obsdesc) :: od
        integer :: i, ngroup
        ngroup=size(od%group)
        do i=1, ngroup
            TRACEX(od%group(i)%group_name)
            TRACEX(od%group(i)%n)
            TRACEX(od%group(i)%idx)
        end do
    end subroutine test_print_obsgroup  !!}}}

end module module_drpsolve
