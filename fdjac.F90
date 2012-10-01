SUBROUTINE fdjac(funcv,x,fvec,df)
    USE nrtype; USE nrutil, ONLY : assert_eq, nrerror
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: fvec
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
    REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df
    INTERFACE
        FUNCTION funcv(x)
        USE nrtype
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: x
        REAL(SP), DIMENSION(size(x)) :: funcv
        END FUNCTION funcv
    END INTERFACE
    REAL(DP), PARAMETER :: EPS=1.0e-04_dp
    INTEGER(I4B) :: j,n
    REAL(DP), DIMENSION(size(x)) :: xsav,xph,h
    n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
    xsav=x
    h=EPS*abs(xsav)
    do j=1,n
        if (h(j) .eq. 0.0) then
            print *, h, xsav
            call nrerror('h is zero!')
        endif
    enddo
    where (h <= 0.0) h=EPS
    xph=xsav+h
    h=xph-xsav
    do j=1,n
        x(j)=xph(j)
        df(:,j)=(funcv(x)-fvec(:))/h(j)
        x(j)=xsav(j)
    end do
    if (minval(abs(df)) .eq. 0.0d0) then
        !print *, 'jac', df(:,j), h(j), x, xsav, fvec 
        !call nrerror('0 in jacobian')
    endif
END SUBROUTINE fdjac
