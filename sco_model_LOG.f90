SUBROUTINE sco_MODEL_LOG(disk_size, corona_size, Tcorona, Tdisk, tau, QPO_frequency, DHext, eta_frac, Nsss, Ssss,Tsss, &
  Nreal, Sreal, Treal, Nimag,Simag, Timag, dTe_mod,dTs_mod, dTe_arg, dTs_arg, Hexo0_out, eta_int)
USE iso_fortran_env, ONLY : WP => REAL64
USE sco_global
USE sco_arrays
IMPLICIT NONE
    ! scalar arguments
    REAL(WP), INTENT(INOUT) :: disk_size, corona_size, Tcorona, Tdisk, tau, QPO_frequency, DHext, eta_frac
    REAL(WP) Emin_adim, Emax_adim, Emin, Emax
!    REAL :: param(5), photar(meshlog-1), EAR(0:meshlog-1), photer(meshlog-1)
!    integer :: near, ifl
!    CHARACTER(4) method
    ! array arguments
    REAL(WP) :: x2(meshlog) , x_use(meshlog-2), xlog(meshlog), xlog_use(meshlog-2)
    ! parameters steady state solution
    REAL(WP) :: c2(meshlog-2), nc, c5, c6, c11, Nesc(meshlog-2), Vc, omega, dxlog, xtot_low, xtot_up, dx
    REAL(WP) :: KN_corr_interpol(meshlog-2)
    INTEGER  Ntri, columns_CC, INFO, nestsol, Nsss, Nreal, Nimag
    REAL(WP) :: L(meshlog-3), U(meshlog-3), D(meshlog-2), CC(meshlog-2), n0(meshlog)
    ! !parameters perturbative solution
    REAL(WP) :: Nescp(meshlog), c2p(meshlog), dn0log(meshlog), dn02log(meshlog)
    REAL(WP) :: KNp_int(meshlog), eta_max, eta, transf, xlog_square(meshlog), xlog_use_square(meshlog-2)
    REAL(WP) :: x_use_square(meshlog-2), powerfact(meshlog-2), x2square(meshlog), exp1(meshlog-2), exp2(meshlog-2)
!     REAL(WP) :: L_subsol(mesh_size-3), U_subsol(mesh_size-3), L_dn0(mesh_size-3), U_dn0(mesh_size-3), Nescp_use(mesh_size-2)
    REAL(WP) :: Q1(meshlog-2), Q2(meshlog-2), Q3(meshlog-2), stau_kn(meshlog-2), factor1, Hexo0
    REAL(WP) :: A1(meshlog-2), A2(meshlog-2), p1, x_withunit(meshlog), rad_sphere, surf, corona_simps
    REAL(WP) :: corona_Lum, corona_Lum_out, vect4(meshlog), Iex01, Iex02, Iex03, Nescp_use(meshlog-2)
    REAL(WP) :: area, tc, to_phys, vect1(meshlog-2), vect2(meshlog-2), vect3(meshlog-2), xlog_trans(meshlog)
    REAL(WP) :: SOLsss(meshlog), SOLreal(meshlog), SOLimag(meshlog), Tsss(meshlog+4), Ssss(meshlog+4)
    REAL(WP) :: Treal(meshlog+4), Sreal(meshlog+4), Timag(meshlog+4), Simag(meshlog+4)
    COMPLEX(WP), DIMENSION(meshlog-3) :: Lp, Up
    COMPLEX(WP), DIMENSION(meshlog-2) :: Dp, sol_ongrid, k0n
    COMPLEX(WP), DIMENSION(meshlog) :: solution, auxiliar, auxiliar2, auxiliar3
    COMPLEX(WP) :: dTe, dTs, dTesimps1, dTesimps2, dTssimps, denom, k0, k1, k2
    REAL(WP) :: dTesimps1_real, dTesimps2_real, dTesimps1_imag, dTesimps2_imag, dTssimps_real, dTssimps_imag
    REAL(WP) :: dTe_mod, dTs_mod, dTe_arg, dTs_arg, Hexo0_out, eta_int
!    REAL(WP), DIMENSION(:), ALLOCATABLE :: x2 , x_use , L, U, D, CC, n0

    call sco_constants(dist, mass, time, energy_norm,  eV2J, keV2J, MeV2J, J2keV, Etrans, kbol, hplanck, c, cc2, me, sigma, stau)


!---We define the energy regime for the BVP solution
    Emin_adim = 1.e-3
    Emax_adim = 40.
    Emin = Emin_adim * Tcorona
    Emax = Emax_adim * Tcorona

!---The output is an array, whose components are evenly spaced numbers between log(Emin/Tcorona) and log(Emax/Tcorona)
    CALL sco_linspace(log(Emin_adim), log(Emax_adim), meshlog, xlog)
    
!---We define the grid of energy spaced evenly on a log scale
    x2 = exp(xlog)

!---We transform input parameters to the internal units
    Tcorona = Tcorona * Etrans
    Tdisk = Tdisk * Etrans
    disk_size = disk_size * 1000. * (1. / dist)
    corona_size = corona_size * 1000. * (1. / dist)
    QPO_frequency = QPO_frequency * time

!---We define the energy step size for the numerical integration
    dxlog = xlog(4) - xlog(3)
    dx = x2(4) - x2(3)
!---We define the integration limits for the energy averaged rms
    xtot_low = 2. * Etrans / Tcorona
    xtot_up = 60. * Etrans / Tcorona

    DO I = 1, meshlog-2
        xlog_use(i) = xlog(i+1)
        x_use(i) = x2(i+1)
    ENDDO

!--------BEGIN OF: construction of the steady state solution----------

    Ntri=meshlog-2 		!dimension of the x_use and then it will be the dimension of the tridiagonal matrix &
                     		!the constant vector of the system that we need to solve

    CALL sco_par(disk_size, corona_size, Tcorona, Tdisk, tau, QPO_frequency, Ntri, x_use, c2, nc, c5, c6, c11, Nesc, &
     Vc, KN_corr_interpol)
    omega = 2.0 * PI * QPO_frequency

!---Preparing to solve the steady state Kompaneets equation (SS) after discretization
    DO I=2, Ntri
      L(I-1) =  1. / (dxlog **2) - (x_use(i)-1.) / (2. * dxlog)  		! sub-diagonal elements
      U(I-1) = 1. / (dxlog **2) + (x_use(i-1)-1.) / (2. * dxlog)  		! super-diagonal elements
    ENDDO
    D = -2. + 2. * x_use - c2  - 2. / (dxlog **2)  				! diagonal elements: x-dependent
    CC = -(x_use)**2 / (exp(x_use * (Tcorona / Tdisk)) - 1.)  			! vector of constants


    columns_CC=1
    CALL dgtsv(Ntri, columns_CC, L, D, U, CC, Ntri, INFO)    			! on exit, CC has the solution of (LDU)*X=CC
    n0(1) = 0.0								! boundary conditions
    n0(meshlog) = 0.0

    DO I =1, Ntri
        n0(I+1) = CC(I)
    ENDDO

!--- n0 is the steady state solution

!---BEGIN OF: construction of the linearized equation

    CALL sco_par(disk_size, corona_size, Tcorona, Tdisk, tau, QPO_frequency, meshlog, x2, c2p, &
    nc, c5, c6, c11, Nescp, Vc, KNp_int)

!---We define the first and second order derivative of n0

    dn0log(1) = (-n0(3) + 4. * n0(2) - 3. * n0(1)) / (2. * dxlog)
    dn0log(meshlog) = (3. * n0(meshlog) - 4. * n0(meshlog-1) + n0(meshlog-2)) / (2. * dxlog)
    dn02log(1) =  (2. * n0(1) - 5. * n0(2) + 4. * n0(3) - n0(4))
    dn02log(meshlog) = (2. * n0(meshlog) + 5. * n0(meshlog-1) + 4. * n0(meshlog-2) + n0(meshlog-3))

    DO I = 1, Ntri
        dn0log(i+1) = (n0(i+2) - n0(i)) / (2.0 * dxlog)
        dn02log(i+1) = (n0(i+2) - 2.0 * n0(i+1) + n0(i)) / (dxlog * dxlog)
        x_use_square(i) = (x_use(i))**2
        xlog_use_square(i) = (xlog_use(i))**2
        Nescp_use(i) = Nescp(i+1)
    ENDDO

    DO I = 1, meshlog
        x2square(i) = (x2(i))**2
        xlog_square(i) = (xlog(i))**2
    ENDDO

!---Argument of the exponentials in eq (A9) of Karpouzas et al 2019
    powerfact = Tcorona * x_use / Tdisk

    exp1 = 1. / (exp(powerfact) - 1)
    exp2 = 1. / (exp(powerfact) - 2. + exp(-powerfact))


    Dp = dcmplx(2. + (dxlog**2) * exp1 * x_use_square/ CC , -c5 * dxlog**2 )

    DO I=1, Ntri -1
      Lp (i) = dcmplx(-1. + (x_use(i+1)-1) * dxlog/2. + (dn0log(i+2) * dxlog)/ CC(i+1) , 0)
      Up (i)= dcmplx(-1. - (x_use(i)-1) * dxlog/2. - (dn0log(i+1) * dxlog)/ CC(I) , 0)
    ENDDO

    stau_kn = (3. / 4.) * stau * KN_corr_interpol      			! klein Nishina correction
    
    Q2 = x_use_square * x_use * CC * stau_kn
    Q1 = CC * x_use_square * stau_kn
    Q3 = CC * x_use_square / Nescp_use

    vect1 = CC * x_use_square * stau_kn
    CALL sco_SIMPSON(Ntri,vect1,xlog_use,Iex01)
    vect2 = CC * x_use * x_use_square * stau_kn
    CALL sco_SIMPSON(Ntri,vect2,xlog_use,Iex02)
    vect3 = CC * x_use_square/ Nescp_use
    CALL sco_SIMPSON(Ntri,vect3,xlog_use,Iex03)

    factor1 = ((Tcorona ** 3) * nc) / (me * c)
    Hexo0 = factor1 * (4. * Iex01 - Iex02)
    eta_max = c11 / Iex03
    eta = eta_frac * eta_max
    denom = dcmplx(4. * factor1 * Iex01 , - (3. / 2.) * omega * Tcorona)   	! denominator in eq (A6) of Karpouzas et al 2019

    eta_int = eta 		! we record eta_int (\tilde\eta)

    k0 = DHext * Hexo0 / denom			! first term in eq (A6) of Karpouzas et al 2019
    k1 = -4. * factor1 / denom			! multiplicative factor of the 3rd term in eq (A6) of Karpouzas et al 2019
    k2 = factor1 / denom			! multiplicative factor of the 2nd term in eq (A6) of Karpouzas et al 2019

    DO I = 1, Ntri
    A1(I) = (dxlog**2) * (-2. -dn0log(i+1)/cc(i) + dn02log(I+1) / CC(I))
    ENDDO

    A2 = (dxlog**2) * (powerfact * x_use_square * exp2 / CC)


    p1 = c6 * eta		! multiplicative factor in eq (A7) of Karpouzas et al 2019 
    k0n = k0 			! the mppinv needs a vector
!---Calculation of the solution

    CALL sco_MPPINV(Lp, Dp, Up, p1, A1, A2, k0n, k1, k2, Q1, Q2, Q3, dxlog, Ntri, sol_ongrid)

    solution(1) = (0,0)			! boundary conditions
    solution(meshlog) = (0,0)
    DO I = 1, Ntri
         solution(i+1) = sol_ongrid(i)
    ENDDO

    SOLreal = REALPART(solution)
    SOLimag = IMAGPART(solution)



    x_withunit = x2 * Tcorona     ! energy grid in internal units
    rad_sphere = disk_size + corona_size
    surf = 4. * pi * rad_sphere ** 2
    vect4 = x2 * n0/ Nescp
    CALL sco_SIMPSON(meshlog,vect4,xlog,corona_simps)
    corona_Lum = surf * (1. - eta) * c * nc * Tcorona ** 2 * corona_simps   ! this is probably wrong
    corona_Lum_out = corona_Lum * keV2J / (time * Etrans)

!---Vector to go from grid units to physical units in ph cm^-2 s^-1 keV^-1 @ 1kpc
    area = 4. * pi * (3e19/dist) **2
    tc = disk_size / (c * tau)
    to_phys = (1.-eta)*Vc*nc/area/(tau+tau**2/3)/tc
    to_phys = to_phys* 1e-4 / dist**2 /time * Etrans

    transf = Etrans/Tcorona 		! parameter used to transform the input grid from keV to internal units
    xlog_trans = x2/transf
    nestsol= meshlog + 4
    SOLsss = n0*transf*to_phys

!---Calculation of the fractional amplitude of the corona temperature oscillation
    auxiliar = x2square * n0 * solution * (3. / 4.) * stau * KNp_int
    CALL sco_SIMPSON(size(xlog),realpart(auxiliar),xlog,dTesimps1_real)
    CALL sco_SIMPSON(size(xlog),imagpart(auxiliar),xlog,dTesimps1_imag)
    dTesimps1 = complex(dTesimps1_real,dTesimps1_imag)

    auxiliar2 = x2 * auxiliar
    CALL sco_SIMPSON(size(xlog),realpart(auxiliar2),xlog,dTesimps2_real)
    CALL sco_SIMPSON(size(xlog),imagpart(auxiliar2),xlog,dTesimps2_imag)
    dTesimps2 = complex(dTesimps2_real,dTesimps2_imag)

    dTe = k0 + k1 * dTesimps1 + k2 * dTesimps2

    dTe_mod = cdabs(dTe)
    dTe_arg = atan2(imagpart(dTe),realpart(dTe))
    

!---Calculation of the fractional amplitude of the seed source temperature oscillation
    auxiliar3 = x2square * n0 * solution / Nescp
    CALL sco_SIMPSON(size(xlog),realpart(auxiliar3),xlog,dTssimps_real)
    CALL sco_SIMPSON(size(xlog),imagpart(auxiliar3),xlog,dTssimps_imag)
    dTssimps=complex(dTssimps_real,dTssimps_imag)

    dTs = p1 * dTssimps
    dTs_mod = cdabs(dTs)
    dTs_arg = atan2(imagpart(dTs),realpart(dTs))

!---Calculation of the external heating rate at thermal equilibrium
    Hexo0_out = Hexo0 * keV2J / (time * Etrans)   ! units : Joules per second


!---Interpolating splines
    CALL sco_InterpolatedUnivariateSpline(meshlog,xlog_trans,SOLsss,nestsol,Nsss,Tsss,Ssss)
    CALL sco_InterpolatedUnivariateSpline(meshlog,xlog_trans,SOLreal,nestsol,Nreal,Treal,Sreal)
    CALL sco_InterpolatedUnivariateSpline(meshlog,xlog_trans,SOLimag,nestsol,Nimag,Timag,Simag)

ENDSUBROUTINE
