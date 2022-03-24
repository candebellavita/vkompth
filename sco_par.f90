SUBROUTINE sco_par(NS_radius, corona_size, Tcorona, Tdisc, tau, QPO_frequency,dim_vect,vector, c2, nc, c5, c6, c11, &
 Nesc, Vc, KN_corr_interpol)
USE iso_fortran_env, ONLY : WP => REAL64
USE sco_global
USE sco_arrays
IMPLICIT NONE
    REAL(WP), INTENT(IN) :: NS_radius, corona_size, Tcorona, Tdisc, tau, QPO_frequency, vector(dim_vect)! Units: km , km , keV , keV , - , Hz
    REAL(WP) Tcorona_keV, XX(dim_vect), tau_kn(dim_vect), X_interpol(dim_vect)
    REAL(WP)  phi(dim_vect), omega, tc
    REAL(WP), INTENT(OUT) :: Vc, nc, c5, c6, c11, Nesc(dim_vect), c2(dim_vect), KN_corr_interpol(dim_vect)
    INTEGER, INTENT(IN) :: dim_vect
    INTEGER, ALLOCATABLE :: X_intermed(:), X_high(:)
    INTEGER dim_intermed, dim_high, k, ier


    call sco_constants(dist, mass, time, energy_norm,  eV2J, keV2J, MeV2J, J2keV, Etrans, kbol, hplanck, c, cc2, me, sigma, stau)


    CALL sco_arrdef(N, T, S)
    omega = 2. * PI * QPO_frequency

    X_interpol = vector * Tcorona / (me * CC2)
    k=3     !degree of spline

    CALL sco_splev(t,n,s,k,X_interpol,KN_corr_interpol,dim_vect,ier)
!	write(*,*) 'KNCORR INTERPOL' , KN_corr_interpol
    tau_kn = (3. / 4.) * tau * KN_corr_interpol
!	write(*,*) 'TAU_KN', tau_kn
    Tcorona_keV = Tcorona / Etrans
    XX = vector * Tcorona_keV / 511.  ! (me * c ^ 2)

    ! In the next DO block, we define the phi array with 1 in all its components and the dimension of the X_intermed array
    ! X_intermed is an array whose components are the components index of the vect array for which 0.1 < XX < 1
    ! X_high is an array whose components are the components index of the vect array for which XX ≥ 1

    dim_intermed = 0
    dim_high = 0
    DO I = 1, dim_vect
        phi(I) = 1.
        IF ((XX(I) .GT. 0.1) .AND. (XX(I) .LT. 1.)) THEN
        dim_intermed = dim_intermed + 1
        ELSEIF (XX(I) .GE. 1.) THEN
        dim_high = dim_high + 1
        ENDIF
    ENDDO
!    write(*,*) dim_intermed, DIM_HIGH
    IF (ALLOCATED(X_intermed)) DEALLOCATE(X_intermed)
    ALLOCATE(X_intermed(dim_intermed))
    IF (ALLOCATED(X_high)) DEALLOCATE(X_high)
    ALLOCATE(X_high(dim_high))

    dim_intermed = 0
    dim_high = 0
    DO I = 1, dim_vect
        IF ((XX(I) .GT. 0.1) .AND. (XX(I) .LT. 1.)) THEN
        dim_intermed = dim_intermed + 1
        X_intermed(dim_intermed) = I
        ELSEIF (XX(I) .GE. 1.) THEN
        dim_high = dim_high + 1
        X_high(dim_high) = I
        ENDIF
    ENDDO

    ! we redefine some components of the phi array
    phi(X_intermed) = (1. - XX(X_intermed))/0.9
    phi(X_high) = 0.0

    Nesc = tau + (1. / 3.) * tau * tau_kn * phi
    Vc = (4. * PI / 3.) * ((NS_radius + corona_size) ** 3. - NS_radius ** 3.)  ! volume of the spherical cell in which the Comptonization occurs
    tc = corona_size / (c * tau)
    c2 = (me * cc2) / (Tcorona * Nesc)
    nc = (me * tc * 2. * PI * Tcorona * 3. * NS_radius**2) / (hplanck * hplanck * hplanck * ((NS_radius + corona_size) ** 3 - &
NS_radius ** 3))
    c5 = (me * cc2 * tc * omega) / Tcorona
    c6 = (Vc * (kbol ** 4) * nc * (Tcorona ** 2)) / (16. * PI * NS_radius**2 * sigma * (Tdisc ** 4.) * tc)
    c11 = (Tdisc ** 4) / ((Vc * (nc * (Tcorona ** 2)) * (kbol ** 4)) / (4. * PI * NS_radius**2 * sigma * tc))


ENDSUBROUTINE sco_par
