SUBROUTINE vkddka(ear,ne,param,IFL,photar,photer)
    IMPLICIT NONE
    INTEGER ifl,ne,mesh_size,i,ENEMAX, MAXNE, mode, neRef, NX, nestsol, ier
    parameter(mesh_size=2999, ENEMAX=10000, MAXNE=10000, NX=299)

    DOUBLE PRECISION ear(0:ne), param(15), photar(ne), photer(ne)
    LOGICAL, save :: firstcall
    INTEGER, save :: pne
    DOUBLE PRECISION, save :: photrms(MAXNE), photlags(MAXNE), photsss(MAXNE), photreal(MAXNE), photimag(MAXNE), photpow(MAXNE)
    DOUBLE PRECISION, save :: pkTs1, pkTe1, pgam1, pLsize1, peta1, pqpo_freq, paf, pDHext1, pear(0:MAXNE)
    DOUBLE PRECISION, save :: pkTs2, pkTe2, pgam2, pLsize2, peta2, pDHext2, pphi

    DOUBLE PRECISION Lsize1, kTe1, kTs1, gam1, eta1, DHext1
    DOUBLE PRECISION Lsize2, kTe2, kTs2, gam2, eta2, DHext2
    DOUBLE PRECISION :: tau, qpo_freq, phi, af, reflag, sss_norm
    DOUBLE PRECISION :: Tsss(mesh_size+4), Ssss(mesh_size+4), Treal(mesh_size+4), Sreal(mesh_size+4)
    DOUBLE PRECISION :: Timag(mesh_size+4), Simag(mesh_size+4), X(NX+1), XM(NX), DX(NX)
    DOUBLE PRECISION :: Xmin, Xmax

    DOUBLE PRECISION :: dTe1_mod, dTs1_mod, dTe1_arg, dTs1_arg, Hexo01_out, eta_int1, dflux1
    DOUBLE PRECISION :: dTe2_mod, dTs2_mod, dTe2_arg, dTs2_arg, Hexo02_out, eta_int2, dflux2
    character(128) dTe1mod,dTs1mod,etaint1,flux1
    character(128) dTe2mod,dTs2mod,etaint2,flux2

    INTEGER Nsss, Nreal, Nimag, dim_int
    DOUBLE PRECISION bwRef(ne+1, 2), fracrms(ne+1), plag_scaled(ne+1)
    DOUBLE PRECISION bwX(NX, 2), fracrmsT(NX), plag_scaledT(NX)
    DOUBLE PRECISION SSS_band(ne+1), Re_band(ne+1), Im_band(ne+1)
    DOUBLE PRECISION SSS_band1(NX), Re_band1(NX), Im_band1(NX)
    DOUBLE PRECISION SSS_band2(NX), Re_band2(NX), Im_band2(NX)
    DOUBLE PRECISION SSS_bandT(NX), Re_bandT(NX), Im_bandT(NX)
    DOUBLE COMPLEX new_complex1(NX), new_complex2(NX), cj

    LOGICAL samecall

    REAL DGFILT
    LOGICAL :: DGQFLT
!    INTEGER :: DGNFLT

    DATA firstcall/.true./

    cj = (0.0,1.0)
    bwX = 0

    !This model does not return model variances.
    photer = 0

    !Header
    if(firstcall)then
        write(*,*) '   ====================================================================='
        write(*,*) '    This is vKompth VERSION the time-dependent Comptonization model'
        write(*,*) '                from Bellavita, Garcia, Mendez & Karpouzas (2022).'
        write(*,*) '    Original models: Karpouzas+(2020) and Garcia+(2021).'
        write(*,*) '    Please cite these papers if you use this model in your publications.'
        write(*,*) '    Feel free to contact us through email or vKompth GitHub page.'
        write(*,*) '   ====================================================================='
        write(*,*) '       Units: rms in fractional units. lags in radians.'
        write(*,*) '       Important: rms and lags normalizations must be fixed to unity.'
        write(*,*) '   ====================================================================='
        firstcall=.false.

        IF (ne.gt.MAXNE) THEN
          write(*,*) '    ERROR: energies must be shorter than ', MAXNE
          write(*,*) '    please recompile the module with a larger MAXNE'
          photar = 0
          return
        END IF

        !mode parameter (1:rms; 2:plag; 3:sss; 4:real; 5:imag)
        IF(.not.DGQFLT(ifl, 'mode')) THEN
          write(*,*) '    WARNING: Spectrum ', IFL, 'lacks mode value in header.'
          write(*,*) '             Assuming time-averaged spectrum (mode=0).'
          mode = 0
        ELSE IF(.not.DGQFLT(ifl, 'QPO')) THEN
          write(*,*) '    WARNING: Spectrum ', IFL, 'lacks QPO value in header'
          write(*,*) '             Assuming time-averaged spectrum (mode=0).'
          mode = 0
        ENDIF

        write(*,*) '     QPO frequency = ', DGFILT(ifl, 'QPO'), ' Hz'
        write(*,*) '   ====================================================================='
    end if

    IF(.not.DGQFLT(ifl, 'mode')) THEN
      mode = 0
    ELSE
      mode = int(DGFILT(ifl, 'mode'))
    ENDIF
    !write(*,*) 'mode=', mode, '  ifl=',ifl, '  ne=',ne, ' kTs=',param(1)
    if ((mode.lt.0.99).OR.(mode.gt.6.01)) mode=0
    reflag = param(15)
    qpo_freq = DGFILT(ifl, 'QPO')

    ! If kTs<0 then, it is -kTs, and the model should return the sss
    if (param(1).lt.0) then
        param(1) = -param(1)
        mode = 0
    endif

    samecall = .FALSE.
    IF (pkTs1.eq.param(1).and.pkTs2.eq.param(2).and.pkTe1.eq.param(3).and.pkTe2.eq.param(4).and. &
        pgam1.eq.param(5).and.pgam2.eq.param(6).and.pLsize1.eq.param(7).and.pLsize2.eq.param(8).and. &
        peta1.eq.param(9).and.peta2.eq.param(10).and. &
        paf.eq.param(11).and.pDHext1.eq.param(12).and.pDHext2.eq.param(13).and.pphi.eq.param(14).and.&
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
            else if (mode.eq.3.or.mode.eq.0) then
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
            else if (mode.eq.6) then
                photar = photpow(1:ne)
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
    paf = param(11)
    pDHext1 = param(12)
    pDHext2 = param(13)
    pphi = param(14)
    pqpo_freq = qpo_freq

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
    af = param(11)
    DHext1 = param(12)
    DHext2 = param(13)
    phi = param(14)
    qpo_freq = pqpo_freq

    ! If gamma<0 then, it is -tau. Otherwise, we get tau from kTe and gamma.
    if (gam1.lt.0) then
        tau = -gam1
    else
        tau = sqrt(2.25+3./((kTe1 / 511.)*((gam1+0.5)**2-2.25))) - 1.5
    endif

    CALL sco_MODEL_LOGdskb(af, Lsize1, kTe1, kTs1, tau, qpo_freq, DHext1, eta1, Nsss, Ssss,Tsss, &
      Nreal, Sreal, Treal, Nimag, Simag, Timag, dTe1_mod, dTs1_mod, dTe1_arg, dTs1_arg, Hexo01_out, eta_int1)

    call PDBVAL('dTe1_mod', dTe1_mod)
    write(dTe1mod,*) dTe1_mod
    call fpmstr('dTe1_mod',dTe1mod)
    call PDBVAL('dTs1_mod', dTs1_mod)
    write(dTs1mod,*) dTs1_mod
    call fpmstr('dTs1_mod',dTs1mod)
    call PDBVAL('eta1_int', eta_int1)
    write(etaint1,*) eta_int1
    call fpmstr('eta1_int',etaint1)

    dim_int = mesh_size+4

    Xmin = log(1e-2)
    Xmax = log(1000.)
    CALL sco_linspace(Xmin, Xmax, NX+1, X)
    X = exp(X)

    IF (ne .LT. ENEMAX) THEN
        DO I= 1, NX
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
    af = param(11)
    DHext1 = param(12)
    DHext2 = param(13)
    phi = param(14)
    qpo_freq = pqpo_freq

    ! We obtain a2 (r_in,2) from a1, kTs1 and kTs2, consistently.
    af = af * (kTs1/kTs2)**(4./3.)

    ! If gamma<0 then, it is -tau. Otherwise, we get tau from kTe and gamma.
    if (gam2.lt.0) then
        tau = -gam2
    else
        tau = sqrt(2.25+3./((kTe2 / 511.)*((gam2+0.5)**2-2.25))) - 1.5
    endif

    CALL sco_MODEL_LOGdskb(af, Lsize2, kTe2, kTs2, tau, qpo_freq, DHext2, eta2, Nsss, Ssss, &
    Tsss, Nreal, Sreal, Treal, Nimag, Simag, Timag, dTe2_mod, dTs2_mod, dTe2_arg, dTs2_arg, Hexo02_out, eta_int2)

    call PDBVAL('dTe2_mod', dTe2_mod)
    write(dTe2mod,*) dTe2_mod
    call fpmstr('dTe2_mod',dTe2mod)
    call PDBVAL('dTs2_mod', dTs2_mod)
    write(dTs2mod,*) dTs2_mod
    call fpmstr('dTs2_mod',dTs2mod)
    call PDBVAL('eta2_int', eta_int2)
    write(etaint2,*) eta_int2
    call fpmstr('eta2_int',etaint2)

    IF (ne .LT. ENEMAX) THEN
        DO I= 1, NX
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
    Re_bandT = REALPART(new_complex1+new_complex2)/SSS_bandT
    Im_bandT = IMAGPART(new_complex1+new_complex2)/SSS_bandT

    nestsol = NX + 4
    ! We now fit Splines to the solution
    XM = 0.5*(X(1:NX)+X(2:NX+1))
    DX = X(2:NX+1)-X(1:NX)
    SSS_bandT = SSS_bandT/DX
    CALL sco_InterpolatedUnivariateSpline(NX,XM,SSS_bandT,nestsol,Nsss,Tsss,Ssss)
    CALL sco_InterpolatedUnivariateSpline(NX,XM,Re_bandT,nestsol,Nreal,Treal,Sreal)
    CALL sco_InterpolatedUnivariateSpline(NX,XM,Im_bandT,nestsol,Nimag,Timag,Simag)

    ! We obtain the total outgoing flux of each component
    CALL sco_SIMPSON(NX,SSS_band1/DX*XM*XM,log(XM),dflux1)
    CALL sco_SIMPSON(NX,SSS_band2/DX*XM*XM,log(XM),dflux2)

    call PDBVAL('flux1', dflux1)
    write(flux1,*) dflux1
    call fpmstr('flux1',flux1)
    call PDBVAL('flux2', dflux2)
    write(flux2,*) dflux2
    call fpmstr('flux2',flux2)

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


        IF (mode.eq.0) THEN
            CALL sco_splev(Tsss,Nsss,Ssss,3,1.D0,sss_norm,1,ier)
        ELSE
            sss_norm = 1.D0
        ENDIF

    ELSE
      write(*,*) 'energies is longer than ', ENEMAX
      write(*,*) 'please recompile using a larger ENEMAX'
    ENDIF

    photrms = fracrms(2:ne+1)*(ear(1:ne)-ear(0:ne-1))
    photpow = 0.5*fracrms(2:ne+1)**2*(ear(1:ne)-ear(0:ne-1))
    photlags = plag_scaled(2:ne+1)*(ear(1:ne)-ear(0:ne-1))
    photsss = SSS_band(2:ne+1)/sss_norm
    photreal = Re_band(2:ne+1)/sss_norm
    photimag = Im_band(2:ne+1)/sss_norm


    if (mode.eq.1) then
        photar = photrms(1:ne)
        return
    else if (mode.eq.2) then
        photar = photlags(1:ne)+reflag*(ear(1:ne)-ear(0:ne-1))
        return
    else if (mode.eq.3.or.mode.eq.0) then
        photar = photsss(1:ne)
        return
    else if (mode.eq.4) then
        photar = photreal(1:ne)
        return
    else if (mode.eq.5) then
        photar = photimag(1:ne)
        return
    else if (mode.eq.6) then
        photar = photpow(1:ne)
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
