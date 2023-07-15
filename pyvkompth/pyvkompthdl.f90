SUBROUTINE pyvkompthdl(ear,ne,param,IFL,photar,photer, dTe_mod, dTs_mod, Hexo0_out, eta_int)
    IMPLICIT NONE
    INTEGER ifl,ne,mesh_size,i,j, ENEMAX, MAXNE
    parameter(mesh_size=2999, ENEMAX=10000, MAXNE=10000)

    DOUBLE PRECISION ear(0:ne), param(8), photar(ne), photer(ne)
    DOUBLE PRECISION, save :: photrms(MAXNE), photlags(MAXNE), photsss(MAXNE), photreal(MAXNE), photimag(MAXNE)
    DOUBLE PRECISION, save :: pkTs, pkTe, ptau, psize, peta_frac, pqpo_freq, paf, pdl

    DOUBLE PRECISION af, Lsize, kTe, kTs
    DOUBLE PRECISION :: disk_size, corona_size, Tcorona, Tdisk, tau, qpo_freq, dL, eta_frac
    DOUBLE PRECISION :: Tsss(mesh_size+4), Ssss(mesh_size+4), Treal(mesh_size+4), Sreal(mesh_size+4)
    DOUBLE PRECISION :: Timag(mesh_size+4), Simag(mesh_size+4)
    DOUBLE PRECISION :: dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out, eta_int

    INTEGER Nsss, Nreal, Nimag, io, rows, dim_int

    DOUBLE PRECISION bandwidth(ne, 2), fracrms(ne), plag_scaled(ne), SSS_band(ne), Re_band(ne), Im_band(ne)

    IF (ne.gt.MAXNE) THEN
        write(*,*) 'energies must be shorter than ', MAXNE
        write(*,*) 'please recompile the module with a larger MAXNE'
        photar = 0
        return
    END IF

    if (pkTs.eq.param(1).and.pkTe.eq.param(2).and.ptau.eq.param(3).and.psize.eq.param(4).and. &
        peta_frac.eq.param(5).and.pqpo_freq.eq.param(6).and.paf.eq.param(7).and.pdl.eq.param(8)) then
        if (ifl.eq.1) then
            !write(*,*) 'rms...'
            photar = photrms(1:ne)
            return
        else if (ifl.eq.2) then
            !write(*,*) 'Lags...'
            photar = photlags(1:ne)
            return
        else if (ifl.eq.3) then
            !write(*,*) 'SSS...'
            photar = photsss(1:ne)
            return
        else if (ifl.eq.4) then
            !write(*,*) 'Real...'
            photar = photreal(1:ne)
            return
        else if (ifl.eq.5) then
            !write(*,*) 'Imag...'
            photar = photimag(1:ne)
            return
        else if (ifl.eq.6) then
            photar = photrms(1:ne)**2/2
            return
        end if
    else
        !write(*,*) pkTs, pkTe, ptau, psize, peta_frac, pqpo_freq, paf
        pkTs = param(1)
        pkTe = param(2)
        ptau = param(3)
        psize = param(4)
        peta_frac = param(5)
        pqpo_freq = param(6)
        paf = param(7)
        pdl = param(8)
        !write(*,*) pkTs, pkTe, ptau, psize, peta_frac, pqpo_freq, paf
    end if

    !write(*,*) 'New call... calculating ', ifl

    kTs = param(1)
    kTe = param(2)
    tau = param(3)
    Lsize = param(4)
    eta_frac = param(5)
    qpo_freq = param(6)
    af = param(7)
    dL = param(8)

    disk_size = af
    corona_size = Lsize
    Tcorona = kTe
    Tdisk = kTs



    CALL sco_model_LOGdskb_dL(disk_size, corona_size, Tcorona, Tdisk, tau, qpo_freq, dL, eta_frac, Nsss, Ssss,Tsss, &
      Nreal, Sreal, Treal, Nimag, Simag, Timag, dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out, eta_int)

    IF (ne .LE. ENEMAX) THEN
        DO I= 1, ne
          bandwidth(i,1) = ear(i-1)
          bandwidth(i,2) = ear(i)
        END DO

        dim_int = mesh_size+4
        CALL sco_band_integrated_amplitude(ne,bandwidth,dim_int, Nsss, Tsss, Ssss, Nreal, Treal, Sreal, Nimag, Timag, Simag, &
            fracrms, plag_scaled, SSS_band,Re_band, Im_band)

        photrms = real(fracrms*(ear(1:ne)-ear(0:ne-1)))
        photlags = real(plag_scaled*(ear(1:ne)-ear(0:ne-1)))
        photsss = real(SSS_band)
        photreal = real(Re_band)
        photimag = real(Im_band)
    ELSE
        write(*,*) 'energies is longer than ', ENEMAX, ' so I will round it'
        DO I = 1, ENEMAX
          bandwidth(i, 1) = float(i-1)/ENEMAX*(ear(ne) - ear(0)) + ear(0)
          bandwidth(i, 2) = float(i)/ENEMAX*(ear(ne) - ear(0)) + ear(0)
        END DO

        dim_int = mesh_size+4
        CALL sco_band_integrated_amplitude(ne,bandwidth,dim_int, Nsss, Tsss, Ssss, Nreal, Treal, Sreal, Nimag, Timag, Simag, &
            fracrms, plag_scaled, SSS_band,Re_band, Im_band)

        DO I = 1, ne
            DO J = 1, ENEMAX
                IF (ear(I) .LE. bandwidth(J, 2) .AND. ear(I) .GE. bandwidth(J, 1)) THEN
                    photrms(i) = fracrms(J)*(ear(I)-ear(I-1))
                    photlags(i) = plag_scaled(J)*(ear(I)-ear(I-1))
                    photsss(i) = SSS_band(J)
                    photreal(i) = Re_band(J)
                    photimag(i) = Im_band(J)
                    EXIT
                END IF
            END DO
        END DO
    END IF


    if (ifl.eq.1) then
        photar = photrms(1:ne)
        return
    else if (ifl.eq.2) then
        photar = photlags(1:ne)
        return
    else if (ifl.eq.3) then
        photar = photsss(1:ne)
        return
    else if (ifl.eq.4) then
        photar = photreal(1:ne)
        return
    else if (ifl.eq.5) then
        photar = photimag(1:ne)
        return
    else if (ifl.eq.6) then
        photar = photrms(1:ne)**2/2
        return
    end if

END SUBROUTINE
