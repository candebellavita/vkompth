MODULE sco_arrays
USE iso_fortran_env, ONLY : WP => REAL64
IMPLICIT NONE
    INTEGER, PARAMETER :: N_kn= 10000
    REAL(WP), DIMENSION(:), ALLOCATABLE :: X_above , KN_corrected
    REAL(WP) break_lim
    REAL(WP) X(N_kn) , KN_corr_tot(N_kn), T(N_kn+4), S(N_kn+4)
    INTEGER, DIMENSION(:), ALLOCATABLE :: X_below_pos
    INTEGER I, dim_below , dim_above, N, Nest

CONTAINS
SUBROUTINE sco_arrdef(N, T, S)
    REAL(WP), DIMENSION(:), ALLOCATABLE :: X_above , KN_corrected
    REAL(WP) break_lim, Xmin, Xmax
    REAL(WP) X(N_kn) , KN_corr_tot(N_kn), T(N_kn+4), S(N_kn+4)
    INTEGER, DIMENSION(:), ALLOCATABLE :: X_below_pos
    INTEGER I, dim_below , dim_above, N



    Xmin = log(1e-5)
    Xmax = log(1000.)
    CALL sco_linspace(Xmin, Xmax, N_kn, X)    ! the output is the X array, whose components are evenly spaced numbers between 1e-5 and 1000. (arbitraty limits)
    X = exp(X)
    break_lim = 1e-4
    ! X_below_pos is an array whose components are the components index of the X array for which X is ≤ break_lim
    ! X_above is the array whose components are the X components that exceed the break_lim
    ! In the next DO block we define the size of both arrays: dim_below and dim_above

    dim_below = 0
    DO I=1, N_kn
        IF(X(I) .LE. break_lim) then
        dim_below = dim_below + 1
        ENDIF
    ENDDO

    dim_above = N_kn - dim_below
    ! Now that we know the size of the arrays, we assigne the memory storage necessary
    IF (ALLOCATED(X_below_pos)) DEALLOCATE(X_below_pos)
    ALLOCATE(X_below_pos (dim_below))
    IF (ALLOCATED(X_above)) DEALLOCATE(X_above)
    ALLOCATE(X_above (dim_above))
    IF (ALLOCATED(KN_corrected)) DEALLOCATE(KN_corrected)
    ALLOCATE (KN_corrected (dim_above))

    dim_below=0
    DO I=1,N_kn
        IF (X(I).LE. break_lim) then
        dim_below= dim_below + 1
        X_below_pos(dim_below) = I
        ENDIF
    ENDDO

    dim_above = 0

    DO I = 1, N_kn
        IF (X(I) .GT. break_lim) THEN
        dim_above = dim_above + 1
        X_above(dim_above) = X(I)
        END IF
    ENDDO

    KN_corrected = ( (1. + X_above) / (X_above ** 3) ) * ( (2. * X_above + 2. * X_above **2) / &
    (1. + 2. * X_above) - log(1. + 2. * X_above)) +log(1. + 2. * X_above) / (2. * X_above) - &
    (1. + 3. * X_above) / ((1. + 2. * X_above) ** 2)

    dim_above = 0
    DO I=1, N_kn
        IF (I.LE.dim_below) THEN
        KN_corr_tot(I) = KN_corrected(1)
        ELSE
        dim_above = dim_above + 1
        KN_corr_tot(I) = KN_corrected(dim_above)
        ENDIF
    ENDDO
!    write(*,*) 'x', x
!    write(*,*) 'KNCORRTOT', KN_CORR_TOT
    Nest= N_kn + 4      ! overestimation of total number of knots
    CALL sco_InterpolatedUnivariateSpline(N_kn, X, kn_corr_tot, Nest, N, T, S)

END SUBROUTINE sco_arrdef


SUBROUTINE sco_linspace(from, to, n,  array)
USE iso_fortran_env, ONLY : WP => REAL64
! Generates evenly spaced numbers from 'from' to 'to' (inclusive)
! from, to : the lower and upper boundaries of the numbers to generate
! n dimension of the array
! array: Array of evenly spaced numbers
    REAL(WP), INTENT(in) :: from, to
    integer :: n, i
    REAL(WP), INTENT(out) :: array(n)
    REAL(WP) :: interval


    interval = to - from

    IF (n.EQ.0) return

    IF (n.EQ.1) THEN

        array(1) = from
        return

    END IF

    DO i=1, n

        array(i) = from + interval * (i - 1) / (n - 1)

    END DO

ENDSUBROUTINE sco_linspace


SUBROUTINE sco_InterpolatedUnivariateSpline(M,X,Y,nest,N,T,C)
USE iso_fortran_env, ONLY : WP => REAL64
    INTEGER, INTENT(IN) :: M, nest
    INTEGER, INTENT(OUT) :: N
    INTEGER :: k, iopt, lwrk, ier
    REAL(WP), INTENT(IN) ::  X(M), Y(M)
    REAL(WP), INTENT(OUT) :: T(Nest), c(nest)
    REAL(WP) :: W(M), s, xb, xe, fp
    REAL(WP), DIMENSION(:), ALLOCATABLE ::  wrk
    INTEGER, ALLOCATABLE :: iwrk(:)

    iopt = 0
    k = 3    ! degree of the spline
    s = 0.   ! smoothing factor
    W = 1.0  ! weights
!   nest = M + k + 1     over estimation of total number of dots
!    wrk_dim = nest * 3 * (k+2) + m * (k+1)
    lwrk = (m*(k+1)+nest*(7+3*k))+1
    IF (ALLOCATED(wrk)) DEALLOCATE(wrk)
    ALLOCATE (wrk (lwrk))

    IF (ALLOCATED(iwrk)) DEALLOCATE(iwrk)
    ALLOCATE (iwrk(nest))

    xb = X(1)*0.99
    xe = X(M)*1.01
!    write(*,*) 'iopt', iopt
!    write(*,*) 'm',m
!    write(*,*) 'x',x
!    write(*,*) 'y',y
!    write(*,*) 'w',w
!    write(*,*) 'xb',xb, 'xe', xe
!    write(*,*) 'k',k, 's', s, 'nest', nest

    CALL sco_curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
!    write(*,*) 't',t
!    write(*,*) 'c',c
!    write(*,*) 'fp',fp
!    write(*,*) 'wrk',wrk
!    write(*,*) 'lwrk', lwrk
!    write(*,*) 'iwrk', iwrk
!    write(*,*) 'ier', ier
END SUBROUTINE sco_InterpolatedUnivariateSpline

ENDMODULE sco_arrays
