SUBROUTINE sco_MODEL(disk_size, corona_size, Tcorona, Tdisk, tau, QPO_frequency, DHext, eta_frac, Nsss, Ssss,Tsss, &
  Nreal, Sreal, Treal, Nimag,Simag, Timag, dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out, eta_int)
USE iso_fortran_env, ONLY : WP => REAL64
USE sco_global
USE sco_arrays
IMPLICIT NONE
    ! scalar arguments
    REAL(WP), INTENT(INOUT) :: disk_size, corona_size, Tcorona, Tdisk, tau, QPO_frequency, DHext, eta_frac
    REAL(WP) Emin_adim, Emax_adim, Emin, Emax
!    CHARACTER(4) method
    ! array arguments
    REAL(WP) :: x2(mesh_size) , x_use(mesh_size-2)
!    REAL :: param(5), photar(mesh_size-1), EAR(0:mesh_size-1), photer(mesh_size-1)
!    integer :: near, ifl
    ! parameters steady state solution
    REAL(WP) :: c2(mesh_size-2), nc, c5, c6, c11, Nesc(mesh_size-2), Vc, omega, dx, xtot_low, xtot_up
    REAL(WP) :: KN_corr_interpol(mesh_size-2)
    INTEGER  Ntri, columns_CC, INFO, nestsol, Nsss, Nreal, Nimag
    REAL(WP) :: L(mesh_size-3), U(mesh_size-3), D(mesh_size-2), CC(mesh_size-2), n0(mesh_size)
    !parameters perturbative solution
    REAL(WP) :: Nescp(mesh_size), c2p(mesh_size), dn0(mesh_size), dn02(mesh_size)
    REAL(WP) :: KNp_int(mesh_size), Hexo0(mesh_size-2), eta_max, eta, transf
    REAL(WP) :: x_use_square(mesh_size-2), powerfact(mesh_size-2), x2square(mesh_size), exp1(mesh_size-2), exp2(mesh_size-2)
    REAL(WP) :: L_subsol(mesh_size-3), U_subsol(mesh_size-3), L_dn0(mesh_size-3), U_dn0(mesh_size-3), Nescp_use(mesh_size-2)
    REAL(WP) :: Q1(mesh_size-2), Q2(mesh_size-2), Q3(mesh_size-2), stau_kn(mesh_size-2), factor1(mesh_size-2)
    REAL(WP) :: A1(mesh_size-2), A2(mesh_size-2), p1, x2_withunit(mesh_size), rad_sphere, surf, corona_simps
    REAL(WP) :: corona_Lum, corona_Lum_out, vect4(mesh_size), Iex01, Iex02, Iex03
    REAL(WP) :: area, tc, to_phys, vect1(mesh_size-2), vect2(mesh_size-2), vect3(mesh_size-2), x2_trans(mesh_size)
    REAL(WP) :: SOLsss(mesh_size), SOLreal(mesh_size), SOLimag(mesh_size), Tsss(mesh_size+4), Ssss(mesh_size+4)
    REAL(WP) :: Treal(mesh_size+4), Sreal(mesh_size+4), Timag(mesh_size+4), Simag(mesh_size+4)
    COMPLEX(WP), DIMENSION(mesh_size-3) :: Lp, Up
    COMPLEX(WP), DIMENSION(mesh_size-2) :: Dp, denom, k0, k1, k2, sol_ongrid
    COMPLEX(WP), DIMENSION(mesh_size) :: solution, auxiliar, auxiliar2, auxiliar3
    COMPLEX(WP) :: dTe, dTs, dTesimps1, dTesimps2, dTssimps
    REAL(WP) :: dTesimps1_real, dTesimps2_real, dTesimps1_imag, dTesimps2_imag, dTssimps_real, dTssimps_imag
    REAL(WP) :: dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out, eta_int
!    REAL(WP), DIMENSION(:), ALLOCATABLE :: x2 , x_use , L, U, D, CC, n0

    call sco_constants(dist, mass, time, energy_norm,  eV2J, keV2J, MeV2J, J2keV, Etrans, kbol, hplanck, c, cc2, me, sigma, stau)
!ccccccccccc QUIZAS ESTO CONVENGA PONERLO ANTES DE ESTA SUBRUTINA, CUANDO LE PASO EL MESH_SIZE PARA NO TENER COMPONENTES EXTRA EN X Y X_USE ccccccccccccc

!    !The Simpson method requires an odd number of mesh points
!    if (mod(mesh_size,2).EQ.0) THEN
!    mesh_size = mesh_size - 1
!    end if
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    ! We define the energy regime for the BVP solution
    Emin_adim = 1.e-3
    Emax_adim = 40.
    Emin = Emin_adim * Tcorona
    Emax = Emax_adim * Tcorona

    ! the output is the X array, whose components are evenly spaced numbers between Emin/Tcorona and Emax/Tcorona
    CALL sco_linspace(Emin_adim, Emax_adim, mesh_size, X2)

    ! param(1) = 2.872684
    ! param(2) = 6.
    ! param(3) = 0.7
    ! param(4) = 0
    ! param(5) = 0
    ! near = mesh_size -1
    ! ear = real(x2)*6.
    ! ifl = 0
    ! CALL donthcomp(ear,near,Param,Ifl,Photar,Photer)
    ! open(unit=32, file='photar.dat')
    ! do i=0, mesh_size-2
    !   write(32,*) 0.5*(ear(i)+ear(i+1)), photar(i+1), 0.5*(-ear(i)+ear(i+1))
    ! enddo
    ! close(32)

    ! We transform input parameters to the internal units
    Tcorona = Tcorona * Etrans
    Tdisk = Tdisk * Etrans
    disk_size = disk_size * 1000. * (1. / dist)
    corona_size = corona_size * 1000. * (1. / dist)
    QPO_frequency = QPO_frequency * time

    ! We define the energy step size for the numerical integration
    dx = x2(4) - x2(3)

    ! We define the integration limits for the energy averaged rms
    xtot_low = 2. * Etrans / Tcorona
    xtot_up = 60. * Etrans / Tcorona

    DO I = 1, mesh_size-2
        x_use(i) = x2(i+1)
    ENDDO

!--------BEGIN OF: construction of the steady state solution----------

    Ntri=mesh_size-2 !dimension of the x_use and then it will be the dimension of the tridiagonal matrix &
                     !the constant vector of the system that we need to solve

    CALL sco_par(disk_size, corona_size, Tcorona, Tdisk, tau, QPO_frequency, Ntri, x_use, c2, nc, c5, c6, c11, Nesc, &
     Vc, KN_corr_interpol)
    omega = 2.0 * PI * QPO_frequency

    ! preparing to solve the steady state Kompaneets equation (SS) after discretization
    L = 1. / (dx * dx) - 1. / (2.0 * dx)  ! sub-diagonal elements: constant
    D = -2. / (x_use * x_use) + 2. / x_use - c2 / (x_use * x_use) - 2. / (dx * dx)  ! diagonal elements: x-dependent
    U = 1. / (dx * dx) + 1. / (2. * dx)  ! super-diagonal elements: constant
    CC = -1. / (exp(x_use * (Tcorona / Tdisk)) - 1.)  ! vector of constants

    columns_CC=1
    CALL dgtsv(Ntri, columns_CC, L, D, U, CC, Ntri, INFO)     ! on exit, CC has the solution of (LDU)*X=CC
    n0(1) = 0.0
    n0(mesh_size) = 0.0

    DO I =1, Ntri
        n0(I+1) = CC(I)
    ENDDO
    transf = Etrans/Tcorona ! parameter used to transform the input grid from keV to internal units

!   Solution of the linearized equation

    CALL sco_par(disk_size, corona_size, Tcorona, Tdisk, tau, QPO_frequency, mesh_size, x2, c2p, &
    nc, c5, c6, c11, Nescp, Vc, KNp_int)

    ! We define the first and second order derivative of n0

    dn0(1) = (-n0(3) + 4. * n0(2) - 3. * n0(1)) / (2. * dx)
    dn0(mesh_size) = (3. * n0(mesh_size) - 4. * n0(mesh_size-1) + n0(mesh_size-2)) / (2. * dx)
    dn02(1) =  (2. * n0(1) - 5. * n0(2) + 4. * n0(3) - n0(4))
    dn02(mesh_size) = (2. * n0(mesh_size) + 5. * n0(mesh_size-1) + 4. * n0(mesh_size-2) + n0(mesh_size-3))

    DO I = 1, Ntri
        dn0(i+1) = (n0(i+2) - n0(i)) / (2.0 * dx)
        dn02(i+1) = (n0(i+2) - 2.0 * n0(i+1) + n0(i)) / (dx * dx)
        x_use_square(i) = (x_use(i))**2
        Nescp_use(i) = Nescp(i+1)
    ENDDO

    DO I = 1, mesh_size
        x2square(i) = (x2(i))**2
    ENDDO

    powerfact = Tcorona * x_use / Tdisk


    exp1 = 1. / (exp(powerfact) - 1)
    exp2 = 1. / (exp(powerfact) - 2. + exp(-powerfact))



    DO I = 1, Ntri-1
        L_subsol(I) = CC(I+1)
        U_subsol(I) = CC(I)
        L_dn0(I) = dn0(i+2)
        U_dn0(I) = dn0(i+1)
    ENDDO

    Dp = dcmplx(2. + (dx**2) * exp1 / CC , -c5 * dx**2 / x_use_square)
    Lp = dcmplx(-1. + dx/2. + (L_dn0 * dx)/L_subsol , 0)
    Up = dcmplx(-1. - dx/2. - (U_dn0 * dx)/U_subsol , 0)

    Q2 = x_use_square * CC
    Q1 = x_use * CC
    Q3 = Q1 / Nescp_use

    vect1 = CC * x_use
    CALL sco_SIMPSON(Ntri,vect1,x_use,Iex01)
    vect2 = CC * x_use_square
    CALL sco_SIMPSON(Ntri,vect2,x_use,Iex02)
    vect3 = CC * x_use/ Nescp_use
    CALL sco_SIMPSON(Ntri,vect3,x_use,Iex03)
!    write(*,*) Iex01, Iex02, Iex03


!    CALL arrdef(M, T, S)
!    X_interpol = x_use * Tcorona / (me * cc2)
!    CALL splev(t,M,s,k,X_interpol,KN_corr_interpol,Ntri,ier)
    stau_kn = (3. / 4.) * stau * KN_corr_interpol        !KN_corr(x_use, Tcorona)  # klein Nishina correction
    factor1 = ((Tcorona ** 3) * stau_kn * nc) / (me * c)
    Hexo0 = factor1 * (4. * Iex01 - Iex02)
    eta_max = c11 / Iex03
    eta = eta_frac * eta_max
    denom = dcmplx(4. * factor1 * Iex01 , - (3. / 2.) * omega * Tcorona)     ! denominator in eq (A6) of Karpouzas et al 2019

    eta_int = eta ! we record eta_int (\tilde\eta)

    k0 = DHext * Hexo0 / denom
    k1 = -4. * factor1 / denom
    k2 = factor1 / denom

    DO I = 1, Ntri
    A1(I) = (dx**2) * (-2. / x_use_square(I) + dn02(I+1) / CC(I))
    ENDDO

    A2 = (dx**2) * (powerfact * exp2 / CC)

    ! Calculation of the solution
    p1 = c6 * eta


    CALL sco_MPPINV(Lp, Dp, Up, p1, A1, A2, k0, k1, k2, Q1, Q2, Q3, dx, Ntri, sol_ongrid)

    solution(1) = (0,0)
    solution(mesh_size) = (0,0)
    DO I = 1, Ntri
         solution(i+1) = sol_ongrid(i)
    ENDDO
    SOLreal = REALPART(solution)
    SOLimag = IMAGPART(solution)

    x2_withunit = x2 * Tcorona     ! energy grid in internal units
    rad_sphere = disk_size + corona_size
    surf = 4. * pi * rad_sphere ** 2
    vect4 = x2 * n0/ Nescp
    CALL sco_SIMPSON(mesh_size,vect4,x2,corona_simps)
    corona_Lum = surf * (1. - eta) * c * nc * Tcorona ** 2 * corona_simps   ! this is probably wrong
    corona_Lum_out = corona_Lum * keV2J / (time * Etrans)

    ! Vector to go from grid units to physical units in ph cm^-2 s^-1 keV^-1 @ 1kpc
    area = 4. * pi * (3e19/dist) **2
    tc = disk_size / (c * tau)
    to_phys = (1.-eta)*Vc*nc/area/(tau+tau**2/3)/tc
    to_phys = to_phys* 1e-4 / dist**2 /time * Etrans

    transf = Etrans/Tcorona ! parameter used to transform the input grid from keV to internal units
    x2_trans = x2/transf
    nestsol= mesh_size + 4
    SOLsss = n0*transf*to_phys

    ! open(unit=88, file='LINLIN.dat')
    ! do i = 1, mesh_size
    !     write(88,*) x2(i)*6, SOLreal(i), solimag(i), SOLSSS(I)
    ! enddo
    ! close(88)

    !fractional amplitude of the corona temperature oscillation
    auxiliar = x2 * n0 * solution
    CALL sco_SIMPSON(size(x2),realpart(auxiliar),x2,dTesimps1_real)
    CALL sco_SIMPSON(size(x2),imagpart(auxiliar),x2,dTesimps1_imag)
    dTesimps1 = complex(dTesimps1_real,dTesimps1_imag)
    auxiliar2 = x2 * auxiliar
    CALL sco_SIMPSON(size(x2),realpart(auxiliar2),x2,dTesimps2_real)
    CALL sco_SIMPSON(size(x2),imagpart(auxiliar2),x2,dTesimps2_imag)
    dTesimps2 = complex(dTesimps2_real,dTesimps2_imag)

    dTe = 0.
    do i=1, size(k0)
        dTe = dTe + (k0(i) + k1(i) * dTesimps1 + k2(i) * dTesimps2) / size(k0)
    enddo
    dTe_mod = cdabs(dTe)
    dTe_arg = atan2(imagpart(dTe),realpart(dTe))

    ! fractional amplitude of the seed source temperature oscillation
    auxiliar3 = auxiliar / Nescp
    CALL sco_SIMPSON(size(x2),realpart(auxiliar3),x2,dTssimps_real)
    CALL sco_SIMPSON(size(x2),imagpart(auxiliar3),x2,dTssimps_imag)
    dTssimps=complex(dTssimps_real,dTssimps_imag)
    dTs = p1 * dTssimps
    dTs_mod = cdabs(dTs)
    dTs_arg = atan2(imagpart(dTs),realpart(dTs))

    Hexo0_out = 0.
    do i=1, size(Hexo0)
        Hexo0_out = Hexo0_out + Hexo0(i) * keV2J / (time * Etrans * size(Hexo0))   ! units : Joules per second
    enddo

    CALL sco_InterpolatedUnivariateSpline(mesh_size,x2_trans,SOLsss,nestsol,Nsss,Tsss,Ssss)
    CALL sco_InterpolatedUnivariateSpline(mesh_size,x2_trans,SOLreal,nestsol,Nreal,Treal,Sreal)
    CALL sco_InterpolatedUnivariateSpline(mesh_size,x2_trans,SOLimag,nestsol,Nimag,Timag,Simag)

END SUBROUTINE sco_MODEL
