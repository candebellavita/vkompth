SUBROUTINE sco_band_integrated_amplitude(rows,bandwidth,dim_int, Nsss, Tsss, Ssss, Nreal, Treal, Sreal, Nimag, Timag, Simag, &
   fracrms, plag_scaled, SSS_band,Re_band, Im_band)
USE iso_fortran_env, ONLY : WP => REAL64
USE sco_global
USE sco_arrays
IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nsss, Nreal, Nimag, dim_int, rows
    INTEGER, PARAMETER :: band_size = 250, k=3, lagrefband=0
    REAL(WP), PARAMETER :: NH = 1.0, NORM = 1.0, EXPTIME = 1.0
    REAL(WP), INTENT(IN) :: bandwidth(rows,2), Tsss(dim_int), Ssss(dim_int), Treal(dim_int), Sreal(dim_int)
    REAL(WP), INTENT(IN) :: Timag(dim_int), Simag(dim_int)
    REAL(WP), INTENT(OUT) :: Re_band(size(bandwidth,1)), Im_band(size(bandwidth,1)), SSS_band(size(bandwidth,1))
    REAL(WP) :: RMS_band(size(bandwidth,1))
    INTEGER :: band_number , J , ier       ! band_numner number of rows of bandwidth
    REAL(WP) :: sss_splev(band_size), real_splev(band_size), imag_splev(band_size), temp1, temp2, temp
    REAL(WP) :: band_mesh(band_size), Re_simps(band_size), Im_simps(band_size), RMS_simps(band_size)
    COMPLEX(WP) :: DNj_band(size(bandwidth,1))
    REAL(WP) :: fracrms(size(bandwidth,1)), angle(size(bandwidth,1)), plag(size(bandwidth,1)), plag_scaled(size(bandwidth,1))
    band_number = size(bandwidth,1)

    !if (.not. present(lagrefband)) lagrefband = 0


    ! Define the np.arrays for output with lengths NCH
    ! DNj_band = 0.0
    ! Re_band = 0.0
    ! Im_band = 0.0
    ! SSS_band = 0.0
    ! RMS_band = 0.0

    DO J=1, band_number

      temp1=bandwidth(j,1)
      temp2=bandwidth(j,2)

      CALL sco_linspace(temp1, temp2, band_size, band_mesh)


      CALL sco_splev(Tsss,Nsss,Ssss,k,band_mesh,sss_splev,band_size,ier)
      CALL sco_splev(Treal,Nreal,Sreal,k,band_mesh,real_splev,band_size,ier)
      CALL sco_splev(Timag,Nimag,Simag,k,band_mesh,imag_splev,band_size,ier)


      CALL sco_SIMPSON(band_size,SSS_splev, band_mesh, temp)
      SSS_band(J) = temp * NORM

      Re_simps = real_splev * sss_splev
      CALL sco_SIMPSON(band_size,Re_simps, band_mesh, Re_band(J))
      Re_band(J) = Re_band(J) * NORM

      Im_simps = imag_splev * sss_splev
      CALL sco_SIMPSON(band_size,Im_simps, band_mesh, Im_band(J))
      Im_band(J) = Im_band(J) * NORM

      DNj_band(J) = dcmplx (Re_band(J) / SSS_band(J) ,  Im_band(J) / SSS_band(J))

      RMS_simps = sss_splev * sqrt(real_splev**2 + imag_splev**2)
      CALL sco_SIMPSON(band_size,RMS_simps, band_mesh, RMS_band(J))
      RMS_band(J)= RMS_band(J)* NORM
    ENDDO

    angle = atan2(IMAGPART(DNj_band), REALPART(DNj_band))


    CALL sco_SMOOTHER(band_number,angle, plag)

    fracrms = abs(DNj_band)


    plag_scaled = plag - plag(lagrefband+1)

ENDSUBROUTINE



! Function smoother :
! input : phi is a numpy ndarray that contains the phase lags for the particular QPO at different energies
! (could be either the full grid or just the centers of the arbitrary energy bands)
SUBROUTINE sco_SMOOTHER(N,phi, nphi)
USE iso_fortran_env, ONLY : WP => REAL64
      REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)
      INTEGER, INTENT(IN) :: N
      REAL(WP), INTENT(IN) :: phi(N)
      REAL(WP), INTENT(OUT) :: nphi(N)
      REAL(WP) :: condition, corrlist(3), compcases(3), compref, patch, dist(3)
      INTEGER, ALLOCATABLE :: pivots(:)
      INTEGER :: ps, piv_dim, pivot, i, j, mindlist

      ps = size(phi)
      nphi = phi
      ! pivots dimention
      piv_dim = 0
      DO I = 1, ps-1
          condition = phi(i) * phi(i+1)
          IF (condition.LT.0) piv_dim = piv_dim +1
      ENDDO
      ! pivots
      IF (ALLOCATED(pivots)) DEALLOCATE(pivots)
      ALLOCATE(pivots(piv_dim))

      J=0
      DO I = 1, ps-1
          condition = phi(i) * phi(i+1)
          IF (condition.LT.0) then
            j=j+1
            pivots(j) = i
          ENDIF
      ENDDO

      corrlist(1) = 2. * pi
      corrlist(2) = -2 * pi
      corrlist(3) = 0.

      DO I=1, piv_dim
        pivot = pivots(i)
        patch = phi(pivot)
        compcases = patch + corrlist
        compref = phi(pivot+1)
        dist = abs(compcases - compref)
        mindlist = minloc(dist, 1)
        DO J = 1, pivot
            nphi(j) = nphi(j) + corrlist(mindlist)
        ENDDO
      ENDDO

ENDSUBROUTINE
