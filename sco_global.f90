MODULE sco_global
USE iso_fortran_env, ONLY : WP => REAL64
    REAL(WP) dist, mass, time, energy_norm, kbol_SI, hplanck_SI, c_SI, me_SI, sigma_SI, stau_SI, eV2J, keV2J, MeV2J, J2keV, Etrans
    REAL(WP) kbol, hplanck, c, cc2, me, sigma, stau
    REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)
    INTEGER, PARAMETER :: mesh_size = 2999 , meshlog = 494

CONTAINS
SUBROUTINE sco_constants(dist, mass, time, energy_norm,  eV2J, keV2J, MeV2J, J2keV, Etrans, kbol, hplanck, c, cc2, me, sigma, stau)
USE iso_fortran_env, ONLY : WP => REAL64
    REAL(WP) dist, mass, time, energy_norm, kbol_SI, hplanck_SI, c_SI, me_SI, sigma_SI, stau_SI, eV2J, keV2J, MeV2J, J2keV, Etrans
    REAL(WP) kbol, hplanck, c, cc2, me, sigma, stau
    ! Here the physical units are arbitrarily defined:
    ! We do that so we can avoid very small or large numbers in the coefficients of the equation
    dist = 1.e-6  ! distance unit
    mass = 1.6021766208e-1  ! mass unit
    time = 1e3  ! time unit
    energy_norm = (1. / mass) * ((1. / dist) ** 2) / ((1. / time) ** 2)  ! energy unit transformation: J -> defined system

    ! Various physical constants in SI

    kbol_SI = 1.3806488e-23  ! Boltzman constant [SI]
    hplanck_SI = 6.62607004e-34  ! Planck constant [SI]
    c_SI = 299792458.0  ! speed of light [m/s]
    me_SI = 9.10938356e-31  ! electron mass [kg]
    sigma_SI = 5.67037321e-8  ! Stefan Boltzman constant [SI]
    stau_SI = 6.6524e-29  ! Thomson cross section [SI]

    ! some usefull transformation tools
    eV2J = 1.6021766208e-19  ! eV to J
    keV2J = 1.6021766208e-16  ! keV to J
    MeV2J = 1.6021766208e-13  ! MeV to J
    J2keV = 6.2415096471204e15
    Etrans = keV2J * energy_norm  ! transformation parameter for keV energy input

    ! Physical constants in the arbitrarily defined unit system
    kbol = kbol_SI * ((1. / mass) * ((1. / dist) ** 2) / ((1. / time) ** 2))
    hplanck = hplanck_SI * ((1. / mass) * ((1. / dist) ** 2) / (1. / time))
    c = c_SI * ((1. / dist) / (1. / time))
    cc2 = c * c  ! speed of light squared
    me = me_SI * (1. / mass)
    sigma = sigma_SI * (1. / mass) * ((1. / time) ** (-3))
    stau = stau_SI * ((1. / dist) ** 2.)

END SUBROUTINE sco_constants

ENDMODULE sco_global
