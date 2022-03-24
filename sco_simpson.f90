SUBROUTINE sco_SIMPSON(N,Y,X,value)
USE iso_fortran_env, ONLY : WP => REAL64
IMPLICIT NONE
      INTEGER, INTENT(IN) :: N        ! Dimension of the array and number of subintervals (IT HAS TO BE EVEN!!!!!)
      REAL(WP), INTENT(IN) :: X(N), Y(N)   ! Y array to be integrated, X points at which y is sampled
      REAL(wp), INTENT(OUT) :: value
      REAL(wp) :: sum_odd, sum_even, step
      INTEGER :: I, M, Nsimp


      step= x(2)-x(1)

      IF (mod(N,2).EQ.0) THEN
          Nsimp = N-1       ! Simpson needs an odd number of dots so we discard the last one
      ELSE
          Nsimp = N
      ENDIF
      M = (Nsimp-1)/2
      sum_even = 0.
      DO I = 1, M
          sum_even = sum_even + Y(2*I)
      ENDDO


      sum_odd = 0.
      DO I = 2, M
          sum_odd = sum_odd + Y(2*I-1)
      ENDDO

      value = (step/3.) * (Y(1) + Y(Nsimp) + 4 * sum_even + 2 * sum_odd)

ENDSUBROUTINE
