SUBROUTINE newt(funcv,x,tolf,check,error)
    USE nrtype; USE nrutil, ONLY : nrerror,vabs
    USE nr, ONLY : fdjac,lnsrch,lubksb,ludcmp
    USE fminln
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
    REAL(SP), INTENT(IN) :: tolf
    REAL(SP), OPTIONAL, INTENT(OUT) :: error
    LOGICAL(LGT), INTENT(OUT) :: check
    INTEGER(I4B), PARAMETER :: MAXITS=1000
    REAL(SP), PARAMETER :: TOLX=epsilon(x),STPMX=100.0
    REAL(SP) :: tolmin
    INTEGER(I4B) :: its
    INTEGER(I4B), DIMENSION(size(x)) :: indx
    REAL(SP) :: d,f,fold,stpmax
    REAL(SP), DIMENSION(size(x)) :: g,p,xold
    REAL(SP), DIMENSION(size(x)), TARGET :: fvec
    REAL(SP), DIMENSION(size(x),size(x)) :: fjac
    INTERFACE
        FUNCTION funcv(x)
        USE nrtype
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: x
        REAL(SP), DIMENSION(size(x)) :: funcv
        END FUNCTION funcv
    END INTERFACE
    tolmin = tolf*1.e-2_dp
    f=fmin(funcv,fvec,x)
    !print *, 'fvec/tolf', fvec, tolf, x
    !if (maxval(abs(fvec(:))) < 0.01_dp*tolf) then
    error = maxval(abs(fvec(:)))
    if (error < tolf) then
        check=.false.
        RETURN
    end if
    stpmax=STPMX*max(vabs(x(:)),real(size(x),sp))
    do its=1,MAXITS
        call fdjac(funcv,x,fvec,fjac)
        g(:)=matmul(fvec(:),fjac(:,:))
        xold(:)=x(:)
        fold=f
        p(:)=-fvec(:)
        call ludcmp(fjac,indx,d)
        call lubksb(fjac,indx,p)
        call lnsrch(funcv,xold,fold,g,p,x,f,stpmax,check,fmin,fvec)
        error = maxval(abs(fvec(:)))
        if (error < tolf) then
            check=.false.
            RETURN
        end if
        if (check) then
            check=(maxval(abs(g(:))*max(abs(x(:)),1.0_dp) / &
                max(f,0.5_dp*size(x))) < tolmin)
            RETURN
        end if
        if (maxval(abs(x(:)-xold(:))/max(abs(x(:)),1.0_dp)) < TOLX) &
            RETURN
    end do
    print *, 'newt values'
    print *, x, xold, fvec, f, fold, tolf, fjac, g, stpmax
    call nrerror('MAXITS exceeded in newt')
END SUBROUTINE newt
