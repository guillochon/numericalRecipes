SUBROUTINE eosjac(funcv,x,df,fmin,fvec)
    USE nrtype; USE nrutil, ONLY : assert_eq, nrerror
    USE eos_helmData, ONLY: nr_eos_deriv
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: x
    REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df
    REAL(SP), OPTIONAL, INTENT(OUT) :: fmin
    REAL(SP), DIMENSION(size(x)), OPTIONAL, INTENT(OUT) :: fvec
    INTERFACE
        FUNCTION funcv(x)
        USE nrtype
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: x
        REAL(SP), DIMENSION(size(x)) :: funcv
        END FUNCTION funcv
    END INTERFACE
    INTEGER(I4B) :: j,n
    REAL(DP), DIMENSION(size(x)) :: xsav,xph,h,heps,val
    n=assert_eq(size(x),size(df,1),size(df,2),'eosjac')
    val=funcv(x)
    if (present(fvec)) fvec = val
    if (present(fmin)) fmin = 0.5_sp*dot_product(val,val) !Rolled this into eosjac to reduce eos calls
    df = nr_eos_deriv
END SUBROUTINE eosjac
