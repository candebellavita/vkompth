SUBROUTINE vkompthbb(ear,ne,param,IFL,photar,photer)
    IMPLICIT NONE
    INTEGER ifl, ne, mesh_size, ENEMAX, MAXNE
    parameter(mesh_size=2999, ENEMAX=10000, MAXNE=10000)

    DOUBLE PRECISION :: ear(0:ne), param(8), photar(ne), photer(ne)
    LOGICAL, save :: firstcall
    INTEGER, save :: pne
    DOUBLE PRECISION, save :: photrms(MAXNE), photlags(MAXNE), photsss(MAXNE), photreal(MAXNE), photimag(MAXNE), photpow(MAXNE)
    DOUBLE PRECISION, save :: pkTs, pkTe, pgam, psize, peta_frac, pqpo_freq, paf, pDHext, pear(0:MAXNE)

    DOUBLE PRECISION :: af, Lsize, kTe, kTs, gam, reflag, sss_norm
    DOUBLE PRECISION :: disk_size, corona_size, Tcorona, Tdisk, tau, qpo_freq, DHext, eta_frac
    DOUBLE PRECISION :: Tsss(mesh_size+4), Ssss(mesh_size+4), Treal(mesh_size+4), Sreal(mesh_size+4)
    DOUBLE PRECISION :: Timag(mesh_size+4), Simag(mesh_size+4)
    DOUBLE PRECISION :: dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out, eta_int

    character(128) dTemod, dTsmod, etaint

    INTEGER Nsss, Nreal, Nimag, dim_int, i, j, mode, neRef, ier
    DOUBLE PRECISION :: bwRef(ne+1, 2), fracrms(ne+1), plag_scaled(ne+1), SSS_band(ne+1), Re_band(ne+1), Im_band(ne+1)
    LOGICAL samecall

    REAL DGFILT
    LOGICAL :: DGQFLT
!    INTEGER :: DGNFLT

    DATA firstcall/.true./
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
    reflag = param(8)
    qpo_freq = DGFILT(ifl, 'QPO')

    ! If kTs<0 then, it is -kTs, and the model should return the sss
    if (param(1).lt.0) then
        param(1) = -param(1)
        mode = 0
    endif

    samecall = .FALSE.
    IF (pkTs.eq.param(1).and.pkTe.eq.param(2).and.pgam.eq.param(3).and.psize.eq.param(4).and. &
        peta_frac.eq.param(5).and.pqpo_freq.eq.qpo_freq.and.paf.eq.param(6).and. &
        pDHext.eq.param(7).and.pne.eq.ne) then
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
    pkTs = param(1)
    pkTe = param(2)
    pgam = param(3)
    psize = param(4)
    peta_frac = param(5)
    pqpo_freq = qpo_freq
    paf = param(6)
    pDHext = param(7)

    kTs = param(1)
    kTe = param(2)
    gam = param(3)
    Lsize = param(4)
    eta_frac = param(5)
    af = param(6)
    DHext = param(7)

    disk_size = af
    corona_size = Lsize
    Tcorona = kTe
    Tdisk = kTs

    ! If gamma<0 then, it is -tau. Otherwise, we get tau from kTe and gamma.
    if (gam.lt.0) then
        tau = -gam
    else
        tau = sqrt(2.25+3./((kTe / 511.)*((gam+0.5)**2-2.25))) - 1.5
    endif

    CALL sco_MODEL_LOGbb(disk_size, corona_size, Tcorona, Tdisk, tau, qpo_freq, DHext, eta_frac, Nsss, Ssss,Tsss, &
      Nreal, Sreal, Treal, Nimag, Simag, Timag, dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out, eta_int)

    call PDBVAL('dTe_mod', dTe_mod)
    write(dTemod,*) dTe_mod
    call fpmstr('dTe_mod',dTemod)
    call PDBVAL('dTs_mod', dTs_mod)
    write(dTsmod,*) dTs_mod
    call fpmstr('dTs_mod',dTsmod)
    call PDBVAL('eta_int', eta_int)
    write(etaint,*) eta_int
    call fpmstr('eta_int',etaint)

    IF (ne .LT. ENEMAX) THEN
        neRef = ne+1
        bwRef(1,1) = 2.0
        bwRef(1,2) = 3.0
        DO I= 1, ne
          bwRef(i+1,1) = ear(i-1)
          bwRef(i+1,2) = ear(i)
        END DO

        dim_int = mesh_size+4
        CALL sco_band_integrated_amplitude(neRef,bwRef,dim_int, Nsss, Tsss, Ssss, Nreal, Treal, Sreal, Nimag, Timag, Simag, &
            fracrms, plag_scaled, SSS_band,Re_band, Im_band)

        IF (mode.eq.0) THEN
          CALL sco_splev(Tsss,Nsss,Ssss,3,1.D0,sss_norm,1,ier)
        ELSE
          sss_norm = 1.D0
        ENDIF

        photrms = fracrms(2:ne+1)*(ear(1:ne)-ear(0:ne-1))
        photpow = 0.5*fracrms(2:ne+1)**2*(ear(1:ne)-ear(0:ne-1))
        photlags = plag_scaled(2:ne+1)*(ear(1:ne)-ear(0:ne-1))
        photsss = SSS_band(2:ne+1)/sss_norm
        photreal = Re_band(2:ne+1)/sss_norm
        photimag = Im_band(2:ne+1)/sss_norm
    ELSE
        write(*,*) 'energies is longer than ', ENEMAX, ' so I will round it'
        neRef = ENEMAX
        bwRef(1,1) = 2.0
        bwRef(1,2) = 3.0
        DO I = 1, ENEMAX-1
          bwRef(i+1, 1) = float(i-1)/(ENEMAX-1)*(ear(ne) - ear(0)) + ear(0)
          bwRef(i+1, 2) = float(i)/(ENEMAX-1)*(ear(ne) - ear(0)) + ear(0)
        END DO

        dim_int = mesh_size+4
        CALL sco_band_integrated_amplitude(neRef,bwRef,dim_int, Nsss, Tsss, Ssss, Nreal, Treal, Sreal, Nimag, Timag, Simag, &
            fracrms, plag_scaled, SSS_band,Re_band, Im_band)

        IF (mode.eq.0) THEN
          CALL sco_splev(Tsss,Nsss,Ssss,3,1.D0,sss_norm,1,ier)
        ELSE
          sss_norm = 1.D0
        ENDIF

        DO I = 1, ne
            DO J = 2, ENEMAX-1
                IF (ear(I) .LE. bwRef(J, 2) .AND. ear(I) .GE. bwRef(J, 1)) THEN
                    photrms(i) = fracrms(J)*(ear(I)-ear(I-1))
                    photpow(i) = 0.5*fracrms(J)**2*(ear(I)-ear(I-1))
                    photlags(i) = plag_scaled(J)*(ear(I)-ear(I-1))
                    photsss(i) = SSS_band(J)/sss_norm
                    photreal(i) = Re_band(J)/sss_norm
                    photimag(i) = Im_band(J)/sss_norm
                END IF
            END DO
        END DO
    END IF

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
