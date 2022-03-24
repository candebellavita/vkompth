
      SUBROUTINE XSDSKB(EAR,NE,PARAM,IDT,PHOTAR,PHOTER)

      INTEGER NE, IDT
      REAL    EAR(0:NE), PARAM(1), PHOTAR(NE), PHOTER(NE)

Cf2py intent(in) NE, IDT
Cf2py intent(in) EAR, PARAM
Cf2py intent(out) PHOTAR, PHOTER

c     Multicolour disk blackbody model used in ISAS, Japan.
c     See Mitsuda et al. PASJ, 36, 741 (1984)
c     & Makishima et al. ApJ 308, 635 1986)
c     Ken Ebisawa 1992/12/22

C     Modified to use double precision and to make numerical
C     integration faster.
C     Ken Ebisawa 1993/07/29

C     Modified the algorithm. by Kazuhisa MITSUDA May 31, 1994
C     A numerical calculation was first done within an accruacy of 0.01 %.
C     Then an interpolation formula was found.
C     The interpolation is precise within 1e-5 level.

      INTEGER I, J
      REAL XN, XH

      double precision TIN, E, photon

C     These coefficients are taken from the spectral fitting program SPFD
C     in ISAS.
C     These are used for Gaussian Integral in the given energy band

      DOUBLE PRECISION  GAUSS(5,2)

      SAVE GAUSS

      DATA  GAUSS /
     $     0.236926885,  0.478628670, 0.568888888, 0.478628670,
     $     0.236926885, -0.906179846, -0.538469310, 0.0,  0.538469310,
     $     0.906179846 /

c suppress a warning message from the compiler
      i = idt

c this model has no errors

      DO I = 1, NE
         PHOTER(I) = 0.0
      ENDDO


      TIN = DBLE(PARAM(1))
      DO 100 I = 1, NE
         XN = (EAR(I)-EAR(I-1))/2.0
         PHOTAR(I) = 0.0
         XH = XN + EAR(I-1)
         DO 200 J = 1, 5
            E = DBLE(XN) * GAUSS(J,2) + DBLE(XH)
            CALL mcdspc(E, TIN, 1.0d0, PHOTON)
            PHOTAR(I) = PHOTAR(I) + REAL(GAUSS(J,1) * PHOTON)
 200     CONTINUE
         PHOTAR(I) = PHOTAR(I) * XN
 100  CONTINUE
      END

C     Multi-Color Disk SPECTRUM
C     SEE MITSUDA ET AL. 1984 PASJ 36, 741.
C     & MAKISHIMA ET AL. APJ 308, 635 1986)
C     NORMALIZATION={RIN(KM)/(D/10KPC)}^2*COS(THETA)
C     TIN=     INNER TEMPERATURE OF THE DISK
C     E  =     ENERGY
C     Rin2 = Normalization factor in the unit of above normalization.
C     PHOTON = PHOTONS/S/KEV/CM2

      subroutine mcdspc(E, Tin, Rin2, Flux)
      DOUBLE PRECISION E, Tin, Rin2, Flux

C  E = X-ray energy (keV)
C  Tin = inner edge color-temperature (keV)
C  Rin2 = inner edge radius(km) ^2 * cos(inclination)/ [D (10kpc)]^2
C  Flux = photon flux, photons/sec/cm^2/keV

      DOUBLE PRECISION normfact
      PARAMETER (normfact=361.0d0)

      DOUBLE PRECISION et
      DOUBLE PRECISION value

      if ( Tin .EQ. 0.0d0 ) then
         Flux = 0.0d0
         return
      endif

      et = E/Tin

      call mcdint(et, value)

      Flux = value *Tin*Tin * Rin2/normfact

      end

      subroutine mcdint(et,value)

      IMPLICIT NONE

      INTEGER NRES
      PARAMETER(NRES=98)
      DOUBLE PRECISION VALUE0, ET0, STEP, A, B, BEKI
      PARAMETER(VALUE0=0.19321556D+03, ET0=0.10000000D-02)
      PARAMETER(STEP=0.60000000D-01, A=0.52876731D+00)
      PARAMETER(B=0.16637530D+01, BEKI=-2.0d0/3.0d0)

      DOUBLE PRECISION et, value
      DOUBLE PRECISION gc(3),gw(3),gn(3),res(NRES)
      DOUBLE PRECISION loget, z, pos
      DOUBLE PRECISION resfact, gaufact
      integer j

      save gc, gw, gn, res

      data gc/  0.78196667D-01, -0.10662020D+01,  0.11924180D+01/
      data gw/  0.52078740D+00,  0.51345700D+00,  0.40779830D+00/
      data gn/  0.37286910D+00,  0.39775528D-01,  0.37766505D-01/
      data (res(j),j=1,40)/
     1   0.96198382D-03, 0.10901181D-02, 0.12310012D-02, 0.13841352D-02,
     1   0.15481583D-02, 0.17210036D-02, 0.18988943D-02, 0.20769390D-02,
     1   0.22484281D-02, 0.24049483D-02, 0.25366202D-02, 0.26316255D-02,
     1   0.26774985D-02, 0.26613059D-02, 0.25708784D-02, 0.23962965D-02,
     1   0.21306550D-02, 0.17725174D-02, 0.13268656D-02, 0.80657672D-03,
     1   0.23337584D-03,-0.36291778D-03,-0.94443569D-03,-0.14678875D-02,
     1  -0.18873741D-02,-0.21588493D-02,-0.22448371D-02,-0.21198179D-02,
     1  -0.17754602D-02,-0.12246034D-02,-0.50414167D-03, 0.32507078D-03,
     1   0.11811065D-02, 0.19673402D-02, 0.25827094D-02, 0.29342526D-02,
     1   0.29517083D-02, 0.26012166D-02, 0.18959062D-02, 0.90128649D-03/
      data (res(j),j=41,80)/
     1  -0.26757144D-03,-0.14567885D-02,-0.24928550D-02,-0.32079776D-02,
     1  -0.34678637D-02,-0.31988217D-02,-0.24080969D-02,-0.11936240D-02,
     1   0.26134145D-03, 0.17117758D-02, 0.28906898D-02, 0.35614435D-02,
     1   0.35711778D-02, 0.28921374D-02, 0.16385898D-02, 0.49857464D-04,
     1  -0.15572671D-02,-0.28578151D-02,-0.35924212D-02,-0.36253044D-02,
     1  -0.29750860D-02,-0.18044436D-02,-0.37796664D-03, 0.10076215D-02,
     1   0.20937327D-02, 0.27090854D-02, 0.28031667D-02, 0.24276576D-02,
     1   0.17175597D-02, 0.81030795D-03,-0.12592304D-03,-0.94888491D-03,
     1  -0.15544816D-02,-0.18831972D-02,-0.19203142D-02,-0.16905849D-02,
     1  -0.12487737D-02,-0.66789911D-03,-0.27079461D-04, 0.59931935D-03/
      data (res(j),j=81,98)/
     1   0.11499748D-02, 0.15816521D-02, 0.18709224D-02, 0.20129966D-02,
     1   0.20184702D-02, 0.19089181D-02, 0.17122289D-02, 0.14583770D-02,
     1   0.11760717D-02, 0.89046768D-03, 0.62190822D-03, 0.38553762D-03,
     1   0.19155022D-03, 0.45837109D-04,-0.49177834D-04,-0.93670762D-04,
     1  -0.89622968D-04,-0.401538532D-04/

      loget = log10(et)
      pos=(loget - log10(ET0))/STEP+1
      j=INT(pos)
      if(j .lt. 1) then
        resfact=res(1)
      elseif(j .ge. NRES) then
        resfact=res(NRES)
      else
        pos=pos-j
        resfact = res(j)*(1.0d0 - pos) + res(j+1)*pos
      endif
      gaufact=1.0d0
      do j=1, 3
        z= (loget - gc(j))/gw(j)
        gaufact = gaufact + gn(j)*exp(-z*z/2.0d0)
      enddo
      value = VALUE0 * (et/ET0)**BEKI
     1  * (1.0d0+a*et**b)*exp(-et)
     1  * gaufact * (1.0d0 + resfact)
      end
