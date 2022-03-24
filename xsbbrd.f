
      SUBROUTINE xsbbrd(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(1), photar(ne), photer(ne)

C---
C XSPEC model subroutine
C black body, normalization proportional to surface area,
C uses simple 2 point approximation to the integral
C---
C see ADDMOD for parameter descriptions
C number of model parameters: 1
C      1      kT      (keV), must be > 0.
C intrinsic energy range:
C      Emin=epsilon(>0),Emax=infinity
C algorithm:
C      n(E) = K 1.0344E-3  E**2 dE / (exp(E/kt)-1)
C normalization:
C      The above normalization is such that K = (RKM)**2 /(D10)**2
C      that is, the norm is the source radius squared (assumeing a spherical
C      surface in units of km)
C      divided by the distance to the source (in units of 10 kpc) squared.
C---
C rashafer 18 Feb 1984
C---
      REAL t, anorm, elow, x, tinv, anormh, alow, ehi, ahi
      INTEGER ie, je

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO


C---
      t = param(1)
      tinv = 1./t
      anorm = 1.0344E-3
      anormh = 0.5*anorm
      elow = ear(0)
      x = elow*tinv
      IF (x.LE.1.e-4) THEN
         alow = elow*t
      ELSEIF (x.GT.60.) THEN
         DO ie = 1, ne
            photar(ie) = 0.0
         ENDDO
         RETURN
      ELSE
         alow = elow*elow/sngl(dble(exp(x))-1.D0)
      ENDIF
      DO ie = 1, ne
         ehi = ear(ie)
         x = ehi*tinv
         IF (x.LE.1.e-4) THEN
            ahi = ehi*t
         ELSEIF (x.GT.60.) THEN
            DO je = ie, ne
               photar(je) = 0.0
            ENDDO
            RETURN
         ELSE
            ahi = ehi*ehi/sngl(dble(exp(x))-1.D0)
         ENDIF
         photar(ie) = anormh*(alow+ahi)*(ehi-elow)
         elow = ehi
         alow = ahi
      ENDDO
      RETURN
      END
