SUBROUTINE vkdualbb(ear,ne,param,IFL,photar,photer)
    IMPLICIT NONE
    INTEGER ifl,ne,mesh_size,i,j, ENEMAX, MAXNE, mode, neRef, NX, nestsol
    parameter(mesh_size=2999, ENEMAX=10000, MAXNE=10000, NX=299)

    REAL ear(0:ne), param(17), photar(ne), photer(ne), reflag
    INTEGER, save :: pne
    REAL, save :: photrms(MAXNE), photlags(MAXNE), photsss(MAXNE), photreal(MAXNE), photimag(MAXNE)
    REAL, save :: pkTs1, pkTe1, pgam1, pLsize1, peta1, pqpo_freq, paf, pDHext1, pear(0:MAXNE)
    REAL, save :: pkTs2, pkTe2, pgam2, pLsize2, peta2, pDHext2, pphi

    DOUBLE PRECISION Lsize1, kTe1, kTs1, gam1, eta1, DHext1
    DOUBLE PRECISION Lsize2, kTe2, kTs2, gam2, eta2, DHext2
    DOUBLE PRECISION :: tau, qpo_freq, phi, af
    DOUBLE PRECISION :: Tsss(mesh_size+4), Ssss(mesh_size+4), Treal(mesh_size+4), Sreal(mesh_size+4)
    DOUBLE PRECISION :: Timag(mesh_size+4), Simag(mesh_size+4), X(NX)
    DOUBLE PRECISION :: Xmin, Xmax

    DOUBLE PRECISION :: dTe1_mod, dTs1_mod, dTe1_arg, dTs1_arg, Hexo01_out
    DOUBLE PRECISION :: dTe2_mod, dTs2_mod, dTe2_arg, dTs2_arg, Hexo02_out
    character*(128) dTe1mod,dTs1mod,pdTe1mod,pdTs1mod
    character*(128) dTe2mod,dTs2mod,pdTe2mod,pdTs2mod

    INTEGER Nsss, Nreal, Nimag, dim_int
    DOUBLE PRECISION bwRef(ne+1, 2), fracrms(ne+1), plag_scaled(ne+1), plag(ne+1)
    DOUBLE PRECISION bwX(NX, 2), fracrmsT(NX), plag_scaledT(NX)
    DOUBLE PRECISION SSS_band(ne+1), Re_band(ne+1), Im_band(ne+1)
    DOUBLE PRECISION SSS_band1(NX), Re_band1(NX), Im_band1(NX)
    DOUBLE PRECISION SSS_band2(NX), Re_band2(NX), Im_band2(NX)
    DOUBLE PRECISION SSS_bandT(NX), Re_bandT(NX), Im_bandT(NX)
    DOUBLE COMPLEX new_complex1(NX), new_complex2(NX), cj

    LOGICAL samecall

    DATA pdTe1mod,pdTs1mod/'dTe1_mod','dTs1_mod'/
    DATA pdTe2mod,pdTs2mod/'dTe2_mod','dTs2_mod'/

    cj = (0.0,1.0)

    !This model does not return model variances.
    photer = 0

    IF (ne.gt.MAXNE) THEN
        write(*,*) 'energies must be shorter than ', MAXNE
        write(*,*) 'please recompile the module with a larger MAXNE'
        photar = 0
        return
    END IF

    !mode parameter (1:rms; 2:plag; 3:sss; 4:real; 5:imag)
    mode = int(param(16))
    IF (mode.gt.5) THEN
        WRITE(*,*) 'Warning: parameter MODE should be between 1 and 5.'
        mode = 5
    ELSE IF (mode.lt.1) THEN
        WRITE(*,*) 'Warning: parameter MODE should be between 1 and 5.'
        mode = 1
    ENDIF
    reflag = param(17)
    samecall = .FALSE.

    IF (pkTs1.eq.param(1).and.pkTs2.eq.param(2).and.pkTe1.eq.param(3).and.pkTe2.eq.param(4).and. &
        pgam1.eq.param(5).and.pgam2.eq.param(6).and.pLsize1.eq.param(7).and.pLsize2.eq.param(8).and. &
        peta1.eq.param(9).and.peta2.eq.param(10).and.pqpo_freq.eq.param(11).and. &
        paf.eq.param(12).and.pDHext1.eq.param(13).and.pDHext2.eq.param(14).and.pphi.eq.param(15).and.&
        pne.eq.ne) then
        i=0
        DO WHILE (ear(i).eq.pear(i).and.i.le.ne)
            i=i+1
        END DO
        if(i.gt.ne) THEN
            samecall = .TRUE.
        END IF
    ENDIF

    IF (samecall) THEN
            if (mode.eq.1) then
                !write(*,*) 'rms...'
                photar = photrms(1:ne)
                return
            else if (mode.eq.2) then
                !write(*,*) 'Lags...'
                photar = photlags(1:ne)+reflag*(ear(1:ne)-ear(0:ne-1))
                return
            else if (ifl.eq.3) then
                !write(*,*) 'SSS...'
                photar = photsss(1:ne)
                return
            else if (mode.eq.4) then
                !write(*,*) 'Real...'
                photar = photreal(1:ne)
                return
            else if (mode.eq.5) then
                !write(*,*) 'Imag...'
                photar = photimag(1:ne)
                return
            end if
    end if

    !write(*,*) 'New call... calculating ', ifl
    pear(0:ne) = ear(0:ne)
    pne = ne
    pkTs1 = param(1)
    pkTs2 = param(2)
    pkTe1 = param(3)
    pkTe2 = param(4)
    pgam1 = param(5)
    pgam2 = param(6)
    pLsize1 = param(7)
    pLsize2 = param(8)
    peta1 = param(9)
    peta2 = param(10)
    pqpo_freq = param(11)
    paf = param(12)
    pDHext1 = param(13)
    pDHext2 = param(14)
    pphi = param(15)

    kTs1 = param(1)
    kTs2 = param(2)
    kTe1 = param(3)
    kTe2 = param(4)
    gam1 = param(5)
    gam2 = param(6)
    Lsize1 = param(7)
    Lsize2 = param(8)
    eta1 = param(9)
    eta2 = param(10)
    qpo_freq = param(11)
    af = param(12)
    DHext1 = param(13)
    DHext2 = param(14)
    phi = param(15)

    ! If gamma<0 then, it is -tau. Otherwise, we get tau from kTe and gamma.
    if (gam1.lt.0) then
        tau = -gam1
    else
        tau = sqrt(2.25+3./((kTe1 / 511.)*((gam1+0.5)**2-2.25))) - 1.5
    endif

    CALL sco_MODEL_LOGbb(af, Lsize1, kTe1, kTs1, tau, qpo_freq, DHext1, eta1, Nsss, Ssss,Tsss, &
      Nreal, Sreal, Treal, Nimag, Simag, Timag, dTe1_mod, dTs1_mod, dTe1_arg, dTs1_arg, Hexo01_out)

    write(dTe1mod,*) dTe1_mod
    call fpmstr(pdTe1mod,dTe1mod)
    write(dTs1mod,*) dTs1_mod
    call fpmstr(pdTs1mod,dTs1mod)

    dim_int = mesh_size+4

    Xmin = log(1e-2)
    Xmax = log(1000.)
    CALL sco_linspace(Xmin, Xmax, NX, X)
    X = exp(X)

    IF (ne .LT. ENEMAX) THEN
        DO I= 1, NX-1
          bwX(i,1) = X(i)
          bwX(i,2) = X(i+1)
        END DO

        CALL sco_band_integrated_amplitude(NX,bwX,dim_int, Nsss, Tsss, Ssss, Nreal, Treal, &
        Sreal, Nimag, Timag, Simag, fracrmsT, plag_scaledT, SSS_band1,Re_band1, Im_band1)
        new_complex1 = (Re_band1+cj*Im_band1)
    ELSE
        write(*,*) 'energies is longer than ', ENEMAX
        write(*,*) 'please recompile using a larger ENEMAX'
    ENDIF

    kTs1 = param(1)
    kTs2 = param(2)
    kTe1 = param(3)
    kTe2 = param(4)
    gam1 = param(5)
    gam2 = param(6)
    Lsize1 = param(7)
    Lsize2 = param(8)
    eta1 = param(9)
    eta2 = param(10)
    qpo_freq = param(11)
    af = param(12)
    DHext1 = param(13)
    DHext2 = param(14)
    phi = param(15)

    ! If gamma<0 then, it is -tau. Otherwise, we get tau from kTe and gamma.
    if (gam2.lt.0) then
        tau = -gam2
    else
        tau = sqrt(2.25+3./((kTe2 / 511.)*((gam2+0.5)**2-2.25))) - 1.5
    endif

    CALL sco_MODEL_LOGbb(af, Lsize2, kTe2, kTs2, tau, qpo_freq, DHext2, eta2, Nsss, Ssss, &
    Tsss, Nreal, Sreal, Treal, Nimag, Simag, Timag, dTe2_mod, dTs2_mod, dTe2_arg, dTs2_arg, Hexo02_out)

    write(dTe2mod,*) dTe2_mod
    call fpmstr(pdTe2mod,dTe2mod)
    write(dTs2mod,*) dTs2_mod
    call fpmstr(pdTs2mod,dTs2mod)

    IF (ne .LT. ENEMAX) THEN
        DO I= 1, NX-1
          bwX(i,1) = X(i)
          bwX(i,2) = X(i+1)
        END DO

        CALL sco_band_integrated_amplitude(NX,bwX,dim_int, Nsss, Tsss, Ssss, Nreal, Treal, &
        Sreal, Nimag, Timag, Simag, fracrmsT, plag_scaledT, SSS_band2,Re_band2, Im_band2)
        new_complex2 = (Re_band2+cj*Im_band2)*exp(cj*phi)
    ELSE
        write(*,*) 'energies is longer than ', ENEMAX
        write(*,*) 'please recompile using a larger ENEMAX'
    ENDIF

    ! We now add the two components
    SSS_bandT = SSS_band1 + SSS_band2
    SSS_bandT(NX) = 1.0
    Re_bandT = REALPART(new_complex1+new_complex2)/SSS_bandT
    Im_bandT = IMAGPART(new_complex1+new_complex2)/SSS_bandT

    nestsol = NX + 4
    ! We now fit Splines to the solution
    CALL sco_InterpolatedUnivariateSpline(NX,X,SSS_bandT,nestsol,Nsss,Tsss,Ssss)
    CALL sco_InterpolatedUnivariateSpline(NX,X,Re_bandT,nestsol,Nreal,Treal,Sreal)
    CALL sco_InterpolatedUnivariateSpline(NX,X,Im_bandT,nestsol,Nimag,Timag,Simag)

    ! We finally obtain the solution on the bwRef vector based on ear
    IF (ne .LT. ENEMAX) THEN
        neRef = ne+1
        bwRef(1,1) = 3.0
        bwRef(1,2) = 4.0
        DO I= 1, ne
          bwRef(i+1,1) = ear(i-1)
          bwRef(i+1,2) = ear(i)
        END DO
        dim_int = mesh_size+4

        CALL sco_band_integrated_amplitude(neRef,bwRef,dim_int, Nsss, Tsss, Ssss, Nreal, Treal, &
                    Sreal, Nimag, Timag, Simag, fracrms, plag_scaled, SSS_band,Re_band, Im_band)

    ELSE
      write(*,*) 'energies is longer than ', ENEMAX
      write(*,*) 'please recompile using a larger ENEMAX'
    ENDIF

    photrms = real(fracrms(2:ne+1)*(ear(1:ne)-ear(0:ne-1)))
    photlags = real(plag_scaled(2:ne+1)*(ear(1:ne)-ear(0:ne-1)))
    photsss = real(SSS_band(2:ne+1))
    photreal = real(Re_band(2:ne+1))
    photimag = real(Im_band(2:ne+1))


    if (mode.eq.1) then
        photar = photrms(1:ne)
        return
    else if (mode.eq.2) then
        photar = photlags(1:ne)+reflag*(ear(1:ne)-ear(0:ne-1))
        return
    else if (mode.eq.3) then
        photar = photsss(1:ne)
        return
    else if (mode.eq.4) then
        photar = photreal(1:ne)
        return
    else if (mode.eq.5) then
        photar = photimag(1:ne)
        return
    end if

END SUBROUTINE





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
