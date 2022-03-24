PROGRAM PRUEBA_DSKB
USE iso_fortran_env, ONLY : WP => REAL64
USE sco_global
USE sco_arrays
IMPLICIT NONE
    REAL(WP) :: disk_size, corona_size, Tcorona, Tdisk, tau, QPO_frequency, DHext, eta_frac
    REAL(WP) :: Tsss(meshlog+4), Ssss(meshlog+4), Treal(meshlog+4), Sreal(meshlog+4)
    REAL(WP) :: Timag(meshlog+4), Simag(meshlog+4)
    REAL(WP) :: dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out
    REAL(WP), ALLOCATABLE :: bandwidth(:,:), fracrms(:), plag_scaled(:), SSS_band(:),Re_band(:), Im_band(:)
    INTEGER :: Nsss, Nreal, Nimag, rows, dim_int

    disk_size     = 10.0 !km
    corona_size   = 7.00 !km
    Tcorona       = 6.00 !keV
    Tdisk         = 0.70 !keV
    tau           = 4.00 !
    QPO_frequency = 400. !Hz
    DHext         = 0.05 !
    eta_frac      = 0.40 !

    CALL sco_MODEL_LOGdskb(disk_size, corona_size, Tcorona, Tdisk, tau, QPO_frequency, DHext, eta_frac, Nsss, Ssss,Tsss, &
      Nreal, Sreal, Treal, Nimag,Simag, Timag, dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out)

    write(*,*) 'dTe_mod = ', dTe_mod
    write(*,*) 'dTs_mod = ', dTs_mod

    rows = 100
    ALLOCATE(bandwidth(rows,2),fracrms(rows), plag_scaled(rows), SSS_band(rows),Re_band(rows), Im_band(rows))
    DO I=1,rows
      !rows log steps from 0.1 to 100.0 keV
      bandwidth(i,1) = 10**( -1 + (i-1)/float(rows)*(3) )
      bandwidth(i,2) = 10**( -1 + i/float(rows)*(3) )
    ENDDO

    dim_int =mesh_size-4
    CALL sco_band_integrated_amplitude(rows,bandwidth,dim_int, Nsss, Tsss, Ssss, Nreal, Treal, Sreal, Nimag, Timag, Simag, &
        fracrms, plag_scaled, SSS_band,Re_band, Im_band)

    OPEN(UNIT=11, FILE= 'outputDSKB.dat')
    WRITE(11,*) '#  E_lo', '  E_hi', 'fractional rms ', ' phase lag ', ' SSS band ', ' Real band ', ' Imag band'
    DO I=1, ROWS
      WRITE(11,*) bandwidth(I,1),bandwidth(I,2),fracrms(I), plag_scaled(I), SSS_band(I), Re_band(I), Im_band(I)
    ENDDO
    CLOSE(11)
END PROGRAM
