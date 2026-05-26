module imf_commons
  use amr_parameters, only: dp
  implicit none

  ! Lookup table: eta_sn(M_cl) and M_SNII(M_cl) via Chabrier + Weidner-Kroupa
  integer, parameter :: nimf = 200
  real(dp) :: imf_mcl  (nimf)   ! cluster mass grid [Msun], log-spaced
  real(dp) :: imf_eta  (nimf)   ! N_SN per Msun of stars formed
  real(dp) :: imf_msnii(nimf)   ! average SN progenitor mass [Msun]

  ! Chabrier (2003) IMF constants
  real(dp), parameter :: mstar_low  = 0.08d0   ! IMF lower stellar mass limit [Msun]
  real(dp), parameter :: mstar_high = 150.0d0  ! IMF upper stellar mass limit [Msun]
  real(dp), parameter :: msn_thresh = 8.0d0    ! minimum SN progenitor mass [Msun]
  real(dp), parameter :: chab_mc    = 0.22d0   ! lognormal characteristic mass [Msun]
  real(dp), parameter :: chab_sig   = 0.57d0   ! lognormal sigma (in log10)
  real(dp), parameter :: chab_alpha = 2.35d0   ! power-law slope for m > 1 Msun

contains

  !------------------------------------------------------------------
  ! Chabrier (2003) IMF: dN/dm (unnormalized)
  ! Lognormal for m <= 1 Msun, power-law for m > 1 Msun,
  ! continuous at m = 1 Msun.
  !------------------------------------------------------------------
  pure function xi_chabrier(m) result(xi)
    real(dp), intent(in) :: m
    real(dp) :: xi
    real(dp), parameter :: cont = &
         exp(-(log10(0.22d0))**2 / (2.0d0*0.57d0**2))   ! continuity at m=1
    if(m <= 1.0d0) then
      xi = (1.0d0/m) * exp(-(log10(m) - log10(chab_mc))**2 &
           & / (2.0d0*chab_sig**2))
    else
      xi = cont * m**(-chab_alpha)
    endif
  end function xi_chabrier

  !------------------------------------------------------------------
  ! Number integral: int_{m1}^{m2} xi(m) dm  (midpoint, log-spaced)
  !------------------------------------------------------------------
  function integ_n(m1, m2) result(res)
    real(dp), intent(in) :: m1, m2
    real(dp) :: res
    integer, parameter :: np = 4000
    integer :: k
    real(dp) :: lm1, lm2, dlm, mk
    res = 0.0d0
    if(m2 <= m1) return
    lm1 = log(m1);  lm2 = log(m2)
    dlm = (lm2 - lm1) / np
    do k = 1, np
      mk  = exp(lm1 + (k - 0.5d0)*dlm)
      res = res + xi_chabrier(mk) * mk * dlm   ! dm = m d(ln m)
    enddo
  end function integ_n

  !------------------------------------------------------------------
  ! Mass integral: int_{m1}^{m2} m * xi(m) dm  (midpoint, log-spaced)
  !------------------------------------------------------------------
  function integ_m(m1, m2) result(res)
    real(dp), intent(in) :: m1, m2
    real(dp) :: res
    integer, parameter :: np = 4000
    integer :: k
    real(dp) :: lm1, lm2, dlm, mk
    res = 0.0d0
    if(m2 <= m1) return
    lm1 = log(m1);  lm2 = log(m2)
    dlm = (lm2 - lm1) / np
    do k = 1, np
      mk  = exp(lm1 + (k - 0.5d0)*dlm)
      res = res + mk * xi_chabrier(mk) * mk * dlm
    enddo
  end function integ_m

  !------------------------------------------------------------------
  ! Weidner-Kroupa (2004): m_max(M_cl) via bisection in log-space.
  ! Solves: integ_n(mmax, mstar_high) = integ_m(mstar_low, mstar_high) / mcl
  ! i.e., one expected star more massive than mmax in a cluster of mass mcl.
  !------------------------------------------------------------------
  function wk_mmax(mcl) result(mmax)
    real(dp), intent(in) :: mcl
    real(dp) :: mmax
    real(dp) :: itot, target_n, nmid, mlo, mhi, mmid
    integer :: iter

    itot     = integ_m(mstar_low, mstar_high)
    target_n = itot / mcl

    ! Cluster large enough: mmax saturates at stellar cap
    if(target_n <= 0.0d0) then
      mmax = mstar_high; return
    endif
    ! Cluster too small: no star can be expected above mstar_low
    if(target_n >= integ_n(mstar_low, mstar_high)) then
      mmax = mstar_low; return
    endif

    ! Bisect in log-space
    mlo = mstar_low
    mhi = mstar_high
    do iter = 1, 80
      mmid = sqrt(mlo * mhi)
      nmid = integ_n(mmid, mstar_high)
      if(nmid > target_n) then
        mlo = mmid
      else
        mhi = mmid
      endif
      if((mhi - mlo) / mhi < 1.0d-7) exit
    enddo
    mmax = min(sqrt(mlo * mhi), mstar_high)
  end function wk_mmax

  !------------------------------------------------------------------
  ! Build the IMF lookup table.
  ! Called once at startup when sf_cluster_sampling=.true.
  !------------------------------------------------------------------
  subroutine init_imf_table(mcl_min, mcl_max)
    use amr_commons, only: myid
    real(dp), intent(in) :: mcl_min, mcl_max   ! cluster mass range [Msun]
    integer :: i
    real(dp) :: mcl, mmax, nsn, mtot, msnii_tot, dlm, eta_lo, eta_hi

    dlm = log(mcl_max / mcl_min) / dble(nimf - 1)

    do i = 1, nimf
      mcl  = mcl_min * exp(dble(i-1) * dlm)
      mmax = wk_mmax(mcl)

      imf_mcl(i) = mcl

      if(mmax > msn_thresh) then
        nsn       = integ_n(msn_thresh, mmax)
        mtot      = integ_m(mstar_low,  mmax)
        msnii_tot = integ_m(msn_thresh, mmax)
        imf_eta  (i) = nsn       / max(mtot,  1.0d-30)
        imf_msnii(i) = msnii_tot / max(nsn,   1.0d-30)
      else
        imf_eta  (i) = 0.0d0
        imf_msnii(i) = 0.0d0
      endif
    enddo

    if(myid == 1) then
      call eta_sn_cluster(1.0d4, eta_lo, mmax)
      call eta_sn_cluster(1.0d6, eta_hi, mmax)
      write(*,*) '>>> CbC: IMF table initialized (Chabrier + Weidner-Kroupa)'
      write(*,'(A,2ES12.4)') '    M_cl range [Msun]        : ', mcl_min, mcl_max
      write(*,'(A,ES12.4)')  '    eta_sn at M_cl=1e4 Msun  : ', eta_lo
      write(*,'(A,ES12.4)')  '    eta_sn at M_cl=1e6 Msun  : ', eta_hi
    endif
  end subroutine init_imf_table

  !------------------------------------------------------------------
  ! Interpolate eta_sn and M_SNII for a given cluster mass [Msun].
  !------------------------------------------------------------------
  subroutine eta_sn_cluster(mcl, eta, msnii)
    real(dp), intent(in)  :: mcl
    real(dp), intent(out) :: eta, msnii
    integer :: i
    real(dp) :: t

    if(mcl <= imf_mcl(1)) then
      eta = imf_eta(1);  msnii = imf_msnii(1);  return
    endif
    if(mcl >= imf_mcl(nimf)) then
      eta = imf_eta(nimf);  msnii = imf_msnii(nimf);  return
    endif
    i = 1
    do while(imf_mcl(i+1) < mcl .and. i < nimf-1)
      i = i + 1
    enddo
    t     = log(mcl / imf_mcl(i)) / log(imf_mcl(i+1) / imf_mcl(i))
    eta   = imf_eta  (i) + t*(imf_eta  (i+1) - imf_eta  (i))
    msnii = imf_msnii(i) + t*(imf_msnii(i+1) - imf_msnii(i))
  end subroutine eta_sn_cluster

end module imf_commons
