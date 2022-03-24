SUBROUTINE vkompthdk(ear,ne,param,IFL,photar,photer)
    IMPLICIT NONE
    INTEGER ifl,ne,mesh_size,i,j, ENEMAX, MAXNE, mode, neRef
    parameter(mesh_size=2999, ENEMAX=10000, MAXNE=10000)

    REAL ear(0:ne), param(10), photar(ne), photer(ne), reflag
    INTEGER, save :: pne
    REAL, save :: photrms(MAXNE), photlags(MAXNE), photsss(MAXNE), photreal(MAXNE), photimag(MAXNE)
    REAL, save :: pkTs, pkTe, pgam, psize, peta_frac, pqpo_freq, paf, pDHext, pear(0:MAXNE)

    DOUBLE PRECISION af, Lsize, kTe, kTs, gam
    DOUBLE PRECISION :: disk_size, corona_size, Tcorona, Tdisk, tau, qpo_freq, DHext, eta_frac
    DOUBLE PRECISION :: Tsss(mesh_size+4), Ssss(mesh_size+4), Treal(mesh_size+4), Sreal(mesh_size+4)
    DOUBLE PRECISION :: Timag(mesh_size+4), Simag(mesh_size+4)
    DOUBLE PRECISION :: dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out

    character*(128) dTemod,dTsmod,pdTemod,pdTsmod

    INTEGER Nsss, Nreal, Nimag, dim_int
    DOUBLE PRECISION bwRef(ne+1, 2), fracrms(ne+1), plag_scaled(ne+1), SSS_band(ne+1), Re_band(ne+1), Im_band(ne+1)
    LOGICAL samecall

    DATA pdTemod,pdTsmod/'dTe_mod','dTs_mod'/

    !This model does not return model variances.
    photer = 0

    IF (ne.gt.MAXNE) THEN
        write(*,*) 'energies must be shorter than ', MAXNE
        write(*,*) 'please recompile the module with a larger MAXNE'
        photar = 0
        return
    END IF

    !mode parameter (1:rms; 2:plag; 3:sss; 4:real; 5:imag)
    mode = int(param(9))
    IF (mode.gt.5) THEN
        WRITE(*,*) 'Warning: parameter MODE should be between 1 and 5.'
        mode = 5
    ELSE IF (mode.lt.1) THEN
        WRITE(*,*) 'Warning: parameter MODE should be between 1 and 5.'
        mode = 1
    ENDIF
    reflag = param(10)
    samecall = .FALSE.

    IF (pkTs.eq.param(1).and.pkTe.eq.param(2).and.pgam.eq.param(3).and.psize.eq.param(4).and. &
        peta_frac.eq.param(5).and.pqpo_freq.eq.param(6).and.paf.eq.param(7).and. &
        pDHext.eq.param(8).and.pne.eq.ne) then
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
    pkTs = param(1)
    pkTe = param(2)
    pgam = param(3)
    psize = param(4)
    peta_frac = param(5)
    pqpo_freq = param(6)
    paf = param(7)
    pDHext = param(8)
    !write(*,*) pkTs, pkTe, ptau, psize, peta_frac, pqpo_freq, paf, pDHext

    kTs = param(1)
    kTe = param(2)
    gam = param(3)
    Lsize = param(4)
    eta_frac = param(5)
    qpo_freq = param(6)
    af = param(7)
    DHext = param(8)

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

    CALL sco_MODEL_LOGdskb(disk_size, corona_size, Tcorona, Tdisk, tau, qpo_freq, DHext, eta_frac, Nsss, Ssss,Tsss, &
      Nreal, Sreal, Treal, Nimag, Simag, Timag, dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out)

    write(dTemod,*) dTe_mod
    call fpmstr(pdTemod,dTemod)
    write(dTsmod,*) dTs_mod
    call fpmstr(pdTsmod,dTsmod)

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

        photrms = real(fracrms(2:ne+1)*(ear(1:ne)-ear(0:ne-1)))
        photlags = real(plag_scaled(2:ne+1)*(ear(1:ne)-ear(0:ne-1)))
        photsss = real(SSS_band(2:ne+1))
        photreal = real(Re_band(2:ne+1))
        photimag = real(Im_band(2:ne+1))
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

        DO I = 1, ne
            DO J = 2, ENEMAX-1
                IF (ear(I) .LE. bwRef(J, 2) .AND. ear(I) .GE. bwRef(J, 1)) THEN
                    photrms(i) = real(fracrms(J)*(ear(I)-ear(I-1)))
                    photlags(i) = real(plag_scaled(J)*(ear(I)-ear(I-1)))
                    photsss(i) = real(SSS_band(J))
                    photreal(i) = real(Re_band(J))
                    photimag(i) = real(Im_band(J))
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
