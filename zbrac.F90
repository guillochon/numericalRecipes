SUBROUTINE zbrac(func,x1,x2,succes)
    USE nrtype; USE nrutil, ONLY : nrerror
    IMPLICIT NONE
    REAL(SP), INTENT(INOUT) :: x1,x2
    LOGICAL(LGT), INTENT(OUT) :: succes
    INTERFACE
        FUNCTION func(x)
        USE nrtype
        IMPLICIT NONE
        REAL(SP), dimension(:), INTENT(IN) :: x
        REAL(SP), dimension(size(x)) :: func
        END FUNCTION func
    END INTERFACE
    INTEGER(I4B), PARAMETER :: NTRY=500
    REAL(SP), PARAMETER :: FACTOR=1.6_sp
    INTEGER(I4B) :: j
    REAL(SP) :: f1,f2
    REAL(SP), dimension(1) :: f1a, f2a, x1a, x2a
    if (x1 == x2) call nrerror('zbrac: you have to guess an initial range')
    x1a(1) = x1
    x2a(1) = x2
    f1a=func(x1a)
    f2a=func(x2a)
    succes=.true.
    do j=1,NTRY
        print *, f1a, f2a
        if ((f1a(1) > 0.0 .and. f2a(1) < 0.0) .or. &
            (f1a(1) < 0.0 .and. f2a(1) > 0.0)) RETURN
        if (abs(f1a(1)) < abs(f2a(1))) then
            x1a(1)=x1a(1)+FACTOR*(x1a(1)-x2a(1))
            x1 = x1a(1)
            f1a=func(x1a)
        else
            x2a(1)=x2a(1)+FACTOR*(x2a(1)-x1a(1))
            x2 = x2a(1)
            f2a=func(x2a)
        end if
    end do
    x1 = x1a(1)
    x2 = x2a(1)
    succes=.false.
END SUBROUTINE zbrac
