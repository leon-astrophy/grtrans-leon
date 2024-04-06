!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Leon's module file for inverse bulk compton scattering 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE COMPTON
IMPLICIT NONE

! Integer !
INTEGER :: nr

! Real !
REAL*8 :: rbh_leon, mbh_leon, aspin_leon

! grmhd variables !
REAL*8, ALLOCATABLE, DIMENSION(:) :: tau_leon
REAL*8, ALLOCATABLE, DIMENSION(:) :: rad_leon
REAL*8, ALLOCATABLE, DIMENSION(:) :: uph_leon
REAL*8, ALLOCATABLE, DIMENSION(:) :: ut_leon
REAL*8, ALLOCATABLE, DIMENSION(:) :: ur_leon
REAL*8, ALLOCATABLE, DIMENSION(:) :: um_leon
REAL*8, ALLOCATABLE, DIMENSION(:) :: up_leon

! Solid angles !
REAL*8, ALLOCATABLE, DIMENSION(:) :: rhop_dir
REAL*8, ALLOCATABLE, DIMENSION(:) :: alphap_dir
REAL*8, ALLOCATABLE, DIMENSION(:) :: cosalpha_p
REAL*8, ALLOCATABLE, DIMENSION(:) :: sinalpha_p
REAL*8, ALLOCATABLE, DIMENSION(:) :: cosrho_p
REAL*8, ALLOCATABLE, DIMENSION(:) :: sinrho_p

! Optical Depth !
REAL*8, PARAMETER :: tau_ics = 1.0d0

! Pi !
REAL*8, PARAMETER :: pi_ics = 4.D0*DATAN(1.D0)

! Use grmhd snapshot ? !
LOGICAL, PARAMETER :: grmhd = .false.

! Check transport ? !
LOGICAL, PARAMETER :: check_transport = .false.

! Use which radiation, 1 = blackbody, 2 = power law !
INTEGER, PARAMETER :: ics_radiation = 2

! Double integral resolution !
INTEGER, PARAMETER ::  n_integral = 32

! Limits !
REAL*8, PARAMETER :: x_low = 0 
REAL*8, PARAMETER :: x_high = pi_ics
REAL*8, PARAMETER :: y_low = 0 
REAL*8, PARAMETER :: y_high = 2*pi_ics

! Step size !
REAL*8, PARAMETER :: drhop = (x_high - x_low)/DBLE(n_integral)
REAL*8, PARAMETER :: dalphap = (y_high - y_low)/DBLE(n_integral)

! Power law index !
REAL*8, PARAMETER :: a_pow = 2.0d0
REAL*8, PARAMETER :: c_pow = 1.0d0
REAL*8, PARAMETER :: s_pow = 1.0d0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Build solid angle arrays 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BUILD_ARRAYS
IMPLICIT NONE

! Integer !
INTEGER :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Allocate !
ALLOCATE(rhop_dir(0:n_integral))
ALLOCATE(alphap_dir(0:n_integral))
ALLOCATE(cosalpha_p(0:n_integral))
ALLOCATE(sinalpha_p(0:n_integral))
ALLOCATE(cosrho_p(0:n_integral))
ALLOCATE(sinrho_p(0:n_integral))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get x-coordinate !
DO i = 0, n_integral
  rhop_dir(i) = x_low + i*drhop
  cosrho_p(i) = DCOS(rhop_dir(i))
  sinrho_p(i) = DSIN(rhop_dir(i))
END DO

! Get y-coordinate !
DO i = 0, n_integral
  alphap_dir(i) = y_low + i*dalphap
  cosalpha_p(i) = DCOS(alphap_dir(i))
  sinalpha_p(i) = DSIN(alphap_dir(i))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Load hdf5 variables into the compton module of grtrans
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HDF5_LOAD
IMPLICIT NONE

! Integer !
INTEGER :: i, nlines
INTEGER :: j

! Real !
REAL*8 :: dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read the number of lines in the file !
nlines = 0 
OPEN (999, file = 'grrmhd') 
DO 
  READ (999,*, END=10) 
  nlines = nlines + 1 
END DO 
10 CLOSE (999)

! Assign nr !
nr = nlines - 1

! Allocate !
ALLOCATE(tau_leon(1:nr))
ALLOCATE(rad_leon(1:nr))
ALLOCATE(uph_leon(1:nr))
ALLOCATE(ut_leon(1:nr))
ALLOCATE(ur_leon(1:nr))
ALLOCATE(um_leon(1:nr))
ALLOCATE(up_leon(1:nr))

! Now read over the file !
OPEN (999, file = 'grrmhd') 
DO i = 0, nlines - 1
  IF(i == 0) THEN
    READ(999,*) nr, dummy, dummy, dummy, dummy
  ELSE
    READ(999,*) rad_leon(i), ut_leon(i), ur_leon(i), um_leon(i), up_leon(i), tau_leon(i), uph_leon(i)
  ENDIF
END DO
CLOSE (999)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given a dimensionless radius r, interpolate Uph and tau
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTERPOLATE_GRMHD (r_in, c_out, tau_out)
IMPLICIT NONE 

! Input !
REAL*8, INTENT(IN) :: r_in
REAL*8, INTENT(OUT) :: c_out, tau_out

! Integer !
INTEGER :: i

! Real !
REAL*8 :: xm, xp
REAL*8 :: ym, yp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Zero everything beyond within the horizon !
IF(r_in <= rbh_leon) THEN 

  tau_out = 0.0d0
  c_out = 0.0d0

ELSE

  ! Loop across the array !
  DO i = 1, nr
    xm = rad_leon(i)
    xp = rad_leon(i+1)
    IF(xm < r_in .AND. xp >= r_in) THEN
      ym = tau_leon(i)
      yp = tau_leon(i+1)
      tau_out = linear_intp (xm, xp, ym, yp, r_in)
      ym = uph_leon(i)
      yp = uph_leon(i+1)
      c_out = linear_intp (xm, xp, ym, yp, r_in)
      EXIT
    END IF
  END DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

REAL*8 function linear_intp (xm, xp, ym, yp, x_in)
implicit none
REAL*8 :: xm, xp, ym, yp, x_in
linear_intp = ym + (x_in - xm)*(yp - ym)/(xp - xm)
end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given a dimensionless radius r, interpolate v1, v2, v3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTERPOLATE_VEL (r_in, th_in, a_in, ut_out, ur_out, um_out, up_out)
IMPLICIT NONE 

! Input !
REAL*8, INTENT(IN) :: r_in, th_in, a_in
REAL*8, INTENT(OUT) :: ut_out, ur_out, um_out, up_out

! Integer !
INTEGER :: i

! Real !
REAL*8 :: norm
REAL*8 :: sigma
REAL*8 :: delta
REAL*8 :: costh
REAL*8 :: sinth
REAL*8 :: xm, xp
REAL*8 :: ym, yp

! Real !
REAL*8, DIMENSION(0:3) :: tmp4
REAL*8, DIMENSION(0:3) :: ucon_ks
REAL*8, DIMENSION(0:3) :: ucon_bl
REAL*8, DIMENSION(0:3,0:3) :: gcov_ks
REAL*8, DIMENSION(0:3,0:3) :: gcov_bl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Constants !
costh = DCOS(th_in)
sinth = DSIN(th_in)
sigma = r_in*r_in + a_in*a_in*costh*costh
delta = r_in*r_in + a_in*a_in - 2*r_in

! Assign metric !
gcov_ks(:,:) = 0.0d0
gcov_ks(0,0) = -1.0d0 + 2.0d0*r_in/sigma
gcov_ks(1,1) = 1.0d0 + 2.0d0*r_in/sigma
gcov_ks(2,2) = sigma
gcov_ks(3,3) = (r_in*r_in + a_in*a_in + 2*a_in*a_in*r_in*sinth*sinth/sigma)*sinth*sinth
gcov_ks(0,1) = 2*r_in/sigma
gcov_ks(0,3) = -2*a_in*r_in*sinth*sinth/sigma
gcov_ks(1,3) = -(1.0d0 + 2.0d0*r_in/sigma)*a_in*sinth*sinth
gcov_ks(1,0) = gcov_ks(0,1)
gcov_ks(3,0) = gcov_ks(0,3)
gcov_ks(3,1) = gcov_ks(1,3)

! Assign metric !
gcov_bl(:,:) = 0.0d0
gcov_bl(0,0) = -1.0d0 + 2.0d0*r_in/sigma
gcov_bl(1,1) = sigma/delta
gcov_bl(2,2) = sigma
gcov_bl(3,3) = (r_in*r_in + a_in*a_in + 2*a_in*a_in*r_in*sinth*sinth/sigma)*sinth*sinth
gcov_bl(0,3) = -2*a_in*r_in*sinth*sinth/sigma
gcov_bl(3,0) = gcov_bl(0,3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Zero everything beyond within the horizon !
IF(r_in <= rbh_leon) THEN

  ! Zero insinde the horizon !
  ucon_ks(0) = 0.0d0
  ucon_ks(1) = 0.0d0
  ucon_ks(2) = 0.0d0
  ucon_ks(3) = 0.0d0

ELSE

  ! Loop across the array !
  DO i = 1, nr
    xm = rad_leon(i)
    xp = rad_leon(i+1)
    IF(xm < r_in .AND. xp >= r_in) THEN
      ym = ut_leon(i)
      yp = ut_leon(i+1)
      ucon_ks(0) = linear_intp (xm, xp, ym, yp, r_in)
      ym = ur_leon(i)
      yp = ur_leon(i+1)
      ucon_ks(1) = linear_intp (xm, xp, ym, yp, r_in)
      ym = um_leon(i)
      yp = um_leon(i+1)
      ucon_ks(2) = linear_intp (xm, xp, ym, yp, r_in)
      ym = up_leon(i)
      yp = up_leon(i+1)
      ucon_ks(3) = linear_intp (xm, xp, ym, yp, r_in)
      EXIT
    END IF
  END DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Normalize in case the norm is not -1 !

tmp4(0) = SUM(gcov_ks(:,0)*ucon_ks(:))
tmp4(1) = SUM(gcov_ks(:,1)*ucon_ks(:))
tmp4(2) = SUM(gcov_ks(:,2)*ucon_ks(:))
tmp4(3) = SUM(gcov_ks(:,3)*ucon_ks(:))
norm = SUM(ucon_ks*tmp4)
ucon_ks = -ucon_ks/SQRT(ABS(norm))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert from ks to bl !

ucon_bl(1) = ucon_ks(1)
ucon_bl(2) = ucon_ks(2)
ucon_bl(0) = ucon_ks(0)-2*r_in/delta*ucon_ks(1)
ucon_bl(3) = ucon_ks(3)-a_in/delta*ucon_ks(1)
norm = SUM(ucon_bl*tmp4)
ucon_bl = -ucon_bl/SQRT(ABS(norm))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ut_out = ucon_bl(0)
ur_out = ucon_bl(1)
um_out = ucon_bl(2)
up_out = ucon_bl(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

REAL*8 function linear_intp (xm, xp, ym, yp, x_in)
implicit none
REAL*8 :: xm, xp, ym, yp, x_in
linear_intp = ym + (x_in - xm)*(yp - ym)/(xp - xm)
end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate ICS with a given nu0z, cosrho0, and gamma_vel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SCATTERING (nu0z, temp, r_in, cosrho0, gamma_vel, stokesi, stokesq, stokesu, i0)
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: nu0z, temp, gamma_vel, r_in
REAL*8, INTENT(OUT) :: stokesi, stokesq, stokesu, i0

! Real variables !
REAL*8 :: sinrho0, D2_lorentz, prefactor, beta_vel
REAL*8 :: cosrho0, cosrho0p, sinrho0p
REAL*8 :: i_local, q_local

! For GRMHD model !
REAL*8 :: c_local, tau_local, rlim

! Dummies !
REAL*8 :: num, den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate beta_vel !
beta_vel = SQRT(1.0d0 - 1.0d0/gamma_vel**2)

! Caluclate angle between electrons and photons 
cosrho0 = MIN(cosrho0, 1.0d0)
cosrho0 = MAX(cosrho0, -1.0d0)
sinrho0 = DSQRT(1.0d0 - cosrho0**2)

! Calculate lorentz factor !
D2_lorentz = 1.0d0/(gamma_vel*(1.0d0 - beta_vel*cosrho0))

! abberration !
num = cosrho0 - beta_vel
den = 1.0d0 - beta_vel*cosrho0
IF(den == 0) THEN
  cosrho0p = 0.0d0
ELSE
  cosrho0p = num/den
END IF        
cosrho0p = MIN(cosrho0p, 1.0d0)
cosrho0p = MAX(cosrho0p, -1.0d0)
sinrho0p = DSQRT(1.0d0 - cosrho0p**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(grmhd) THEN

  ! Interpolate !
  CALL INTERPOLATE_GRMHD(r_in, c_local, tau_local)

  ! Initial field !
  CALL RADIATION_FIELD(nu0z, temp, c_local, i0)

  ! Prefactor !
  prefactor = (3.0d0/16.0d0/pi_ics)*(1.0d0 - beta_vel*cosrho0)*tau_local*D2_lorentz**3

  ! Scattering !
  CALL DBINTEGRAL(cosrho0p, sinrho0p, nu0z, temp, beta_vel, gamma_vel, D2_lorentz, c_local, i_local, q_local)

  ! bad interpolation data !
  IF(isnan(gamma_vel)) THEN
    i_local = 0.0d0
    q_local = 0.0D0
    prefactor = 0.0d0
  END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE

  ! Set prefactor !
  c_local = c_pow*(r_in/rbh_leon)**(-a_pow)

  ! Initial field !
  CALL RADIATION_FIELD(nu0z, temp, c_local, i0)

  ! Prefactor !
  prefactor = (3.0d0/16.0d0/pi_ics)*(1.0d0 - beta_vel*cosrho0)*tau_ics*D2_lorentz**3

  ! Scattering !
  CALL DBINTEGRAL(cosrho0p, sinrho0p, nu0z, temp, beta_vel, gamma_vel, D2_lorentz, c_local, i_local, q_local)

END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! output !
stokesi = prefactor*i_local
stokesq = prefactor*q_local
stokesu = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fortran subroutine for performing double integral
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DBINTEGRAL(cosrho0p_in, sinrho0p_in, nu0z_in, temp_in, betavel_in, gamma_in, D2_in, c_fac, i_out, q_out)
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: cosrho0p_in, sinrho0p_in, nu0z_in, temp_in, betavel_in, gamma_in, D2_in, c_fac
REAL*8, INTENT(OUT) :: i_out, q_out

! INTEGER ! 
INTEGER :: i, j, k

! REAL !
REAL*8 :: fin_m1, fin_c, fin_p1
REAL*8 :: fout_m1, fout_c, fout_p1
REAL*8 :: gin_m1, gin_c, gin_p1
REAL*8 :: gout_m1, gout_c, gout_p1

! REAL ! 
REAL*8 :: d1_m1, d1_c, d1_p1
REAL*8 :: nu_m1, nu_c, nu_p1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! evaluate the integral !

! Initialize !
i_out = 0
q_out = 0

! Loop over the domain !
DO i = 1, n_integral, 2

  ! Initialize !
  fout_m1 = 0
  fout_c = 0
  fout_p1 = 0
  gout_m1 = 0
  gout_c = 0
  gout_p1 = 0

  ! Compute doppler shifted frequency and the d1 factor !
  CALL doppler_d1(gamma_in, betavel_in, cosrho_p(i-1), d1_m1)
  CALL doppler_d1(gamma_in, betavel_in, cosrho_p(i), d1_c)
  CALL doppler_d1(gamma_in, betavel_in, cosrho_p(i+1), d1_p1)
  CALL frequency_shift(nu0z_in, d1_m1, D2_in, nu_m1)
  CALL frequency_shift(nu0z_in, d1_c, D2_in, nu_c)
  CALL frequency_shift(nu0z_in, d1_p1, D2_in, nu_p1)

  ! Loop over the domain !
  DO j = 1, n_integral, 2

    ! For the scattering intensity !
    CALL INTEGRAND_IP(cosalpha_p(j-1), cosrho_p(i-1), sinrho_p(i-1), & 
                      cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, temp_in, c_fac, fin_m1)
    CALL INTEGRAND_IP(cosalpha_p(j), cosrho_p(i-1), sinrho_p(i-1), & 
                      cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, temp_in, c_fac, fin_c)
    CALL INTEGRAND_IP(cosalpha_p(j+1), cosrho_p(i-1), sinrho_p(i-1), & 
                      cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, temp_in, c_fac, fin_p1)
    fout_m1 = fout_m1 + dalphap*(fin_m1 + 4.0d0*fin_c + fin_p1)/3.0d0

    CALL INTEGRAND_IP(cosalpha_p(j-1), cosrho_p(i), sinrho_p(i), & 
                      cosrho0p_in, sinrho0p_in, d1_c, nu_c, temp_in, c_fac, fin_m1)
    CALL INTEGRAND_IP(cosalpha_p(j), cosrho_p(i), sinrho_p(i), & 
                      cosrho0p_in, sinrho0p_in, d1_c, nu_c, temp_in, c_fac, fin_c)
    CALL INTEGRAND_IP(cosalpha_p(j+1), cosrho_p(i), sinrho_p(i), & 
                      cosrho0p_in, sinrho0p_in, d1_c, nu_c, temp_in, c_fac, fin_p1)
    fout_c = fout_c + dalphap*(fin_m1 + 4.0d0*fin_c + fin_p1)/3.0d0

    CALL INTEGRAND_IP(cosalpha_p(j-1), cosrho_p(i+1), sinrho_p(i+1), & 
                      cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, temp_in, c_fac, fin_m1)
    CALL INTEGRAND_IP(cosalpha_p(j), cosrho_p(i+1), sinrho_p(i+1), & 
                      cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, temp_in, c_fac, fin_c)
    CALL INTEGRAND_IP(cosalpha_p(j+1), cosrho_p(i+1), sinrho_p(i+1), & 
                      cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, temp_in, c_fac, fin_p1)
    fout_p1 = fout_p1 + dalphap*(fin_m1 + 4.0d0*fin_c + fin_p1)/3.0d0

    ! For the stokes U and Q !
    CALL INTEGRAND_QP(alphap_dir(j-1), cosalpha_p(j-1), cosrho_p(i-1), sinrho_p(i-1), &
                      cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, temp_in, c_fac, gin_m1)
    CALL INTEGRAND_QP(alphap_dir(j), cosalpha_p(j), cosrho_p(i-1), sinrho_p(i-1), &
                      cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, temp_in, c_fac, gin_c)
    CALL INTEGRAND_QP(alphap_dir(j+1), cosalpha_p(j+1), cosrho_p(i-1), sinrho_p(i-1), &
                      cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, temp_in, c_fac, gin_p1)
    gout_m1 = gout_m1 + dalphap*(gin_m1 + 4.0d0*gin_c + gin_p1)/3.0d0

    CALL INTEGRAND_QP(alphap_dir(j-1), cosalpha_p(j-1), cosrho_p(i), sinrho_p(i), &
                      cosrho0p_in, sinrho0p_in, d1_c, nu_c, temp_in, c_fac, gin_m1)
    CALL INTEGRAND_QP(alphap_dir(j), cosalpha_p(j), cosrho_p(i), sinrho_p(i), &
                      cosrho0p_in, sinrho0p_in, d1_c, nu_c, temp_in, c_fac, gin_c)
    CALL INTEGRAND_QP(alphap_dir(j+1), cosalpha_p(j+1), cosrho_p(i), sinrho_p(i), &
                      cosrho0p_in, sinrho0p_in, d1_c, nu_c, temp_in, c_fac, gin_p1)
    gout_c = gout_c + dalphap*(gin_m1 + 4.0d0*gin_c + gin_p1)/3.0d0    

    CALL INTEGRAND_QP(alphap_dir(j-1), cosalpha_p(j-1), cosrho_p(i+1), sinrho_p(i+1), & 
                      cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, temp_in, c_fac, gin_m1)
    CALL INTEGRAND_QP(alphap_dir(j), cosalpha_p(j), cosrho_p(i+1), sinrho_p(i+1), & 
                      cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, temp_in, c_fac, gin_c)
    CALL INTEGRAND_QP(alphap_dir(j+1), cosalpha_p(j+1), cosrho_p(i+1), sinrho_p(i+1), & 
                      cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, temp_in, c_fac, gin_p1)
    gout_p1 = gout_p1 + dalphap*(gin_m1 + 4.0d0*gin_c + gin_p1)/3.0d0  

  END DO
  
  ! Do the outer integration !
  i_out = i_out + drhop*(fout_m1 + 4.0d0*fout_c + fout_p1)/3.0d0
  q_out = q_out + drhop*(gout_m1 + 4.0d0*gout_c + gout_p1)/3.0d0

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing the inner integrand of the scattering intensity function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTEGRAND_IP(cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, & 
                        d1_in, nu_in, temp_in, c_fac, int_out)
IMPLICIT NONE

! REAL !
REAL*8, INTENT(IN) :: cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, d1_in, nu_in, temp_in, c_fac
REAL*8, INTENT(OUT) :: int_out

! REAL !
REAL*8 :: cosw_p, sinw_p, inu_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
CALL W_P (cosalphap_in, cosrho0p_in, sinrho0p_in, cosrhop_in, sinrhop_in, cosw_p, sinw_p)

! Get intensity !
CALL RADIATION_FIELD(nu_in, temp_in, c_fac, inu_shift)

! Assign !
int_out = sinrhop_in*(1.0d0+cosw_p**2)*inu_shift*d1_in**(3.0d0) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing the inner integrand of the Stokes Q and U parameter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTEGRAND_QP(alphap_in, cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, & 
                        d1_in, nu_in, temp_in, c_fac, int_out)
IMPLICIT NONE

! REAL !
REAL*8, INTENT(IN) :: alphap_in, cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, d1_in, nu_in, temp_in, c_fac
REAL*8, INTENT(OUT) :: int_out

! REAL !
REAL*8 :: inu_shift
REAL*8 :: cosw_p, sinw_p
REAL*8 :: cos2eta, sin2eta
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
CALL W_P (cosalphap_in, cosrho0p_in, sinrho0p_in, cosrhop_in, sinrhop_in, cosw_p, sinw_p)
CALL TWOETA(alphap_in, cosrhop_in, cosrho0p_in, sinrho0p_in, cosw_p, sinw_p, cos2eta, sin2eta)

! Get intensity !
CALL RADIATION_FIELD(nu_in, temp_in, c_fac, inu_shift)

! Assign !
int_out = sinrhop_in*(1.0d0-cosw_p**2)*cos2eta*inu_shift*d1_in**(3.0d0)  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Choose radiation law according to user input 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RADIATION_FIELD(nu_in, temp_in, c_fac, inu_out)
use phys_constants, only: h,c2,k
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: nu_in, temp_in, c_fac

! Output !
REAL*8, INTENT(OUT) :: inu_out

! Choose !
IF(ics_radiation == 1) THEN
  IF(h*nu_in/k/temp_in.lt.1d-6) THEN
    inu_out = 2d0*nu_in*nu_in*k*temp_in/c2
  ELSE
    inu_out = 2d0*h*nu_in/c2*nu_in*nu_in/(exp(h*nu_in/k/temp_in)-1d0)
  END IF
ELSE
  ! The proportionality constant C is 1 !
  ! So all stokes parameters are in units of C !
  inu_out = c_fac*nu_in**(-s_pow)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing the doppler factor d1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE doppler_d1(gamma_in, betavel_in, cosrhop_in, d1_out)
IMPLICIT NONE

! REAL !
REAL*8, INTENT(IN) :: gamma_in, betavel_in, cosrhop_in
REAL*8, INTENT(OUT) :: d1_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
d1_out = 1.0d0/(gamma_in*(1.0d0 + betavel_in*cosrhop_in))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing doppler shifted photon frequency 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE frequency_shift(nu_in, d1_in, d2_in, nu_out)
IMPLICIT NONE

! REAL !
REAL*8, INTENT(IN) :: nu_in, d1_in, d2_in
REAL*8, INTENT(OUT) :: nu_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
nu_out = nu_in/d1_in/d2_in

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate cow_p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE W_P(cosalphap_in, cosrho0p_in, sinrho0p_in, cosrhop_in, sinrhop_in, coswp_out, sinwp_out)
IMPLICIT NONE
  
! Input !
REAL*8, INTENT(IN) :: cosalphap_in, cosrho0p_in, sinrho0p_in, cosrhop_in, sinrhop_in

! Output !
REAL*8, INTENT(OUT) :: coswp_out, sinwp_out
  
! Calculate 
coswp_out = cosalphap_in*sinrhop_in*sinrho0p_in + cosrhop_in*cosrho0p_in
coswp_out = MIN(coswp_out, 1.0d0)
coswp_out = MAX(coswp_out, -1.0d0)
sinwp_out = DSQRT(1.0d0 - coswp_out**2)
  
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate cos2eta and sin2eta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TWOETA(alphap_in, cosrhop_in, cosrho0p_in, sinrho0p_in, coswp_in, sinwp_in, cos2eta_out, sin2eta_out)
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: alphap_in, cosrhop_in, cosrho0p_in, sinrho0p_in, coswp_in, sinwp_in

! Output !
REAL*8, INTENT(OUT) :: cos2eta_out, sin2eta_out

! Local !
REAL*8 :: coseta, sineta, eta, num, den

! Calculate 
num = (cosrhop_in - cosrho0p_in*coswp_in)
den = sinrho0p_in*sinwp_in
IF(den == 0) THEN
  coseta = 0.0d0
ELSE
  coseta = num/den
END IF
coseta = MIN(coseta, 1.0d0)
coseta = MAX(coseta, -1.0d0)
eta = DACOS(coseta)

! Angle !
IF(alphap_in >= pi_ics) THEN
  eta = 2.0D0*pi_ics - eta
END IF
sineta = DSIN(eta)

! Get angles !
cos2eta_out = 2.0d0*coseta**2 - 1.0d0
sin2eta_out = 2.0d0*sineta*coseta

END SUBROUTINE

END MODULE
