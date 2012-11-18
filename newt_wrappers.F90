MODULE newt_wrappers

implicit none

CONTAINS

    subroutine run_newt(funcv, N, tolf, input, output, error)
        USE nrtype
        USE nr
        IMPLICIT NONE
        INTERFACE
            FUNCTION funcv(x)
            USE nrtype
            implicit none
            real, DIMENSION(:), INTENT(IN) :: x
            real, dimension(size(x)) :: funcv
            END FUNCTION funcv
        END INTERFACE
        integer, intent(IN) :: N
        double precision, intent(IN) :: tolf
        integer :: i
        double precision, dimension(N) :: x,f
        double precision, dimension(N), intent(IN) :: input
        double precision, dimension(N), intent(OUT) :: output
        double precision, optional, intent(OUT) :: error
        logical :: check
        double precision :: dummy
        x = input
        if (present(error)) then
            call newt(funcv,x,tolf,check,error)
        else
            call newt(funcv,x,tolf,check,dummy)
        endif
        f=funcv(x)
        output = x
        if (.not. present(error) .and. check) then
            write(*,*) 'Convergence problems.'
            call Driver_abortFlash('Convergence problems in run_newt')
        endif
        !write(*,'(1x,a5,t10,a1,t22,a1)') 'Index','x','f'
        !do i=1,N
        !    write(*,'(1x,i2,2x,2f24.10)') i,x(i),f(i)
        !end do

        return
    END subroutine run_newt

    !subroutine run_annewt(funcv, N, tolf, input, output, failed)
    !    USE nrtype
    !    USE nrutil, ONLY : nrerror
    !    USE nr
    !    IMPLICIT NONE
    !    INTERFACE
    !        FUNCTION funcv(x)
    !        USE nrtype
    !        implicit none
    !        real, DIMENSION(:), INTENT(IN) :: x
    !        real, dimension(size(x)) :: funcv
    !        END FUNCTION funcv
    !    END INTERFACE
    !    integer, parameter :: maxcalls = 10
    !    integer, intent(IN) :: N
    !    double precision, intent(IN) :: tolf
    !    logical, intent(OUT) :: failed
    !    integer :: i, ncalls
    !    double precision, dimension(N) :: x,f,xsav
    !    double precision, dimension(N), intent(IN) :: input
    !    double precision, dimension(N), intent(OUT) :: output
    !    double precision :: new_tolf
    !    logical :: check
    !    x = input
    !    xsav = x
    !    call annewt(funcv,x,tolf,check,failed)
    !    !failed = .true.
    !    !check = .true.
    !    !ncalls = 0
    !    !new_tolf = tolf
    !    !do while (failed .eq. .true. .or. check .eq. .true.)
    !    !    if (mod(ncalls,2) .eq. 0) then
    !    !        x = (1.0 + ncalls/maxcalls*1.d-5)*x !Just slightly perturb
    !    !    else
    !    !        x = (1.0 - ncalls/maxcalls*1.d-5)*x !Just slightly perturb
    !    !    endif
    !    !    call annewt(funcv,x,new_tolf,check,failed)
    !    !    ncalls = ncalls + 1
    !    !    if (ncalls .gt. maxcalls) exit
    !    !    new_tolf = 2.d0*new_tolf
    !    !enddo
    !    !if (ncalls .gt. maxcalls) then
    !    !    if (check .eq. .true.) write(*,*) 'Convergence problems.'
    !    !    if (failed .eq. .true.) write(*,*) 'Maxits exceeded.'
    !    !    call nrerror('exceeded maxcalls in run_annewt')
    !    !endif
    !    !f=funcv(x)
    !    output = x
    !    !if (check) then
    !    !    write(*,*) 'Convergence problems.'
    !    !    call Driver_abortFlash('Convergence problems in run_annewt')
    !    !endif
    !    !if (failed) then
    !    !    write(*,*) 'Annewt failed'
    !    !    call Driver_abortFlash('Failed in run_annewt')
    !    !endif
    !    !write(*,'(1x,a5,t10,a1,t22,a1)') 'Index','x','f'
    !    !do i=1,N
    !    !    write(*,'(1x,i2,2x,2f24.10)') i,x(i),f(i)
    !    !end do

    !    return
    !END subroutine run_annewt

END MODULE newt_wrappers
