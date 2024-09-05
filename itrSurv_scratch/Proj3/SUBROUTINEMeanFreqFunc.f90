
SUBROUTINE MeanFreqFunc(nt_endpoint, nt_survival, tp_survival, tp_endpoint, Nj_endpoint, Oj_endpoint, Nj_survival, Oj_survival, survRE, dRhat, mu, dmu)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nt_endpoint
  INTEGER, INTENT(IN) :: nt_survival
  REAL(dp), DIMENSION(1:nt_survival), INTENT(IN) :: tp_survival
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(OUT) :: tp_endpoint
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(IN) :: Nj_endpoint
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(IN) :: Oj_endpoint
  REAL(dp), DIMENSION(1:nt_survival), INTENT(IN) :: Nj_survival
  REAL(dp), DIMENSION(1:nt_survival), INTENT(IN) :: Oj_survival
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(OUT) :: survRE
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(OUT) :: dRhat
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(OUT) :: mu
  REAL(dp), DIMENSION(1:nt_endpoint), INTENT(OUT) :: dmu

  INTEGER :: i, j

  REAL(dp), DIMENSION(1:nt_endpoint) :: muTmp
  REAL(dp), DIMENSION(1:nt_survival) :: survD

  PRINT *, "******************** mean frequency function ********************"
  mu = 0.d0
  dmu = 0.d0
  
  PRINT *, "******************** obtain survival KM (deaths) ********************"
  call kaplan(nt_survival, Nj_survival, Oj_survival, survD)
  
  DO i = 1, nt_endpoint ! Loop over all elements in tp_endpoint (or another array of length nt_endpoint)

    ! to get the time points that have individuals at risk
    ! if number at risk is equal to 0, then skip to next time point index (i)
    IF (Nj_endpoint(i) .LT. 1d-8) CYCLE
    PRINT *, "******************** dRhat ********************"
    dRhat(i) = Oj_endpoint(i) / Nj_endpoint(i) ! dRhat at each timepoint
    
    PRINT *, "******************** Recurrent KM (using RE times) ********************"
    DO j = 1, tp_survival-1  ! Loop over elements in tp_survival, up to tp_survival-1
      IF (tp_survival(j) <= tp_endpoint(i) .AND. tp_endpoint(i) < tp_survival(j+1)) THEN
        survRE(i) = survD(j) ! Assign survD(j) to survRE(i) if condition is met
        EXIT                 ! Exit the inner loop if condition is met
      END IF
    END DO
  END DO
  
  mu(1) = 0
  if (nt_endpoint .LT. 2) RETURN
  DO i = 2, nt_endpoint
    mu(i) = mu(i-1) + survRE(i) * dRhat(i)
  END DO

END SUBROUTINE MeanFreqFunc




