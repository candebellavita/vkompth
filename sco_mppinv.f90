SUBROUTINE sco_MPPINV(L, D, U, p1, C1, C2, k0,k1, k2, Q1, Q2, Q3, dx, N, sol_ongrid)
USE iso_fortran_env, ONLY : WP => REAL64
IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    COMPLEX(WP), DIMENSION(N-1) :: L, U
    COMPLEX(WP), DIMENSION(N) :: D, sol_ongrid, k0
    COMPLEX(WP), DIMENSION(N,N) :: A
    COMPLEX(WP) :: k1, k2
    REAL(WP), INTENT(INOUT) :: C1(N), C2(N), p1, Q1(N), Q2(N), Q3(N), dx
    INTEGER I, J, IPIV(N), INFO, NRHS, M
    INTEGER, DIMENSION(:), ALLOCATABLE :: JJ1, JJ2
    CHARACTER TRANS

    A = (0.,0.)
    do i = 1, N
        A(i,i) = D(i)
    enddo
    do i = 1, N-1
        A(i,i+1) = U(I)
    enddo
    do i = 1, N-1
        A(i+1,i) = L(I)
    enddo


    DO i = 1, N
        A(i,1) = A(i,1) - (dx / 3.) * (c1(i) * k2 * Q2(1) + c2(i) * p1 * Q3(1) + c1(I) * k1 * Q1(1))
    ENDDO

    DO I = 1,N
        A(i,N) = A(i,N) -(dx / 3.) * (c1(i) * k2 * Q2(N) + c2(i) * p1 * Q3(N) + c1(i) * k1 * Q1(N))
    ENDDO

!        DO J = 2, N-2, 2
!             A(i,j) = A(i,j) - (4. * dx / 3.) * (c1(i) * k2(i) * Q2(J) + c2(I) * p1 * Q3(J) + c1(I) * k1(I) * Q1(J))
!             IF (j.lt.(n-2)) then
!                 A(i,j+1) = A(i,j+1) - (2. * dx / 3.) * (c1(I) * k2(I) * Q2(J+1) + c2(I) * p1 * Q3(J+1) + c1(I) * k1(I) * Q1(J+1))
!             ENDIF
!        ENDDO
    IF (mod(N,2).EQ.2) THEN
      M = (N-2)/2
    ELSE
      M =(N-1)/2
    ENDIF
    IF (ALLOCATED(JJ1)) DEALLOCATE(JJ1)
    IF (ALLOCATED(JJ2)) DEALLOCATE(JJ2)
    ALLOCATE (JJ1(M), JJ2(M-1))
    j=0
    DO I = 1, M-1
      J = 2*I
      JJ1(i) = J
      JJ2(I) = J+1
    ENDDO
      JJ1(M) = 2*M
    ! j=0
    ! DO I = 1, M-1
    !   J = 2*I+1
    !   JJ2(i) = J
    ! ENDDO


    DO I = 1, N
      A(I, JJ1) = A(i,JJ1) - (4. * dx / 3.) * (c1(i) * k2 * Q2(Jj1) + c2(I) * p1 * Q3(JJ1) + c1(I) * k1 * Q1(JJ1))
      A(i, JJ2) = A(i,JJ2) - (2. * dx / 3.) * (c1(I) * k2 * Q2(JJ2) + c2(I) * p1 * Q3(JJ2) + c1(I) * k1 * Q1(JJ2))
    ENDDO

    CALL ZGETRF( N, N, A, N, IPIV, INFO )


    sol_ongrid = C1 * k0
    TRANS = 'N'      ! Specifies the form of the system of equations: in this case No Transpose
    NRHS = 1

    CALL ZGETRS(TRANS, N, NRHS, A, N, IPIV, sol_ongrid, N, INFO)

    !write(*,*) 'INFO', INFO




ENDSUBROUTINE
