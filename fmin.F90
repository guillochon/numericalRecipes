MODULE fminln
    USE nrtype
    USE nrutil, ONLY : nrerror
CONTAINS
!BL
    FUNCTION fmin(funcv,fvec,x)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: x
    REAL(SP), DIMENSION(size(x)), INTENT(OUT) :: fvec
    REAL(SP) :: fmin
    INTERFACE
        FUNCTION funcv(x)
        USE nrtype
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: x
        REAL(SP), DIMENSION(size(x)) :: funcv
        END FUNCTION funcv
    END INTERFACE
    fvec=funcv(x)
    fmin=0.5_sp*dot_product(fvec,fvec)
    END FUNCTION fmin
END MODULE fminln
