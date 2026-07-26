module imf_commons
  use amr_parameters, only: dp
  implicit none

  ! Lookup table: eta_sn(M_cl) and M_SNII(M_cl) via Chabrier + Weidner-Kroupa
  integer, parameter :: nimf = 200
  real(dp) :: imf_mcl  (nimf)   ! cluster mass grid [Msun], log-spaced
  real(dp) :: imf_eta  (nimf)   ! N_SN per Msun of stars formed
  real(dp) :: imf_msnii(nimf)   ! average SN progenitor mass [Msun]

  ! Maximum clusters drawn from CMF per SF event (safety ceiling)
  integer, parameter :: max_clusters_per_event = 10000

  ! Per-cell characteristic cluster mass [Msun] and lognormal width, for
  ! sf_cluster_kernel='lognormal'. Filled in starform2 (sf_model=1 branch),
  ! consumed by the CbC pre-pass in star_formation. Indexed like flag2.
  real(dp), allocatable :: cbc_mclchar_buf(:), cbc_sigmalnm_buf(:)

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

  !------------------------------------------------------------------
  ! Decompose a SF event into clusters from CMF (dN/dM ~ M^{-2}).
  ! Environment-dependent caps privilege massive clusters at high Mj:
  !   mmin_eff = max(mmin0, mmin0*(Mj/Mj_ref)^delta)
  !   mmax_eff = min(mmax,  f_cap*Mj)
  ! Mass is exactly conserved: remainder added to last cluster.
  ! If M_sf < mmin_eff, returns one cluster with the full mass.
  ! All masses in Msun.
  !------------------------------------------------------------------
  subroutine sample_cmf_clusters(M_sf, Mj, mmin0, mmax, f_cap, Mj_ref, delta, &
                                  seed, cluster_masses, n_cl)
    use random
    real(dp), intent(in)    :: M_sf, Mj, mmin0, mmax, f_cap, Mj_ref, delta
    integer,  intent(inout) :: seed(IRandNumSize)
    real(dp), intent(out)   :: cluster_masses(max_clusters_per_event)
    integer,  intent(out)   :: n_cl

    real(dp)     :: M_rem, mmin_eff, mmax_eff, mmax_draw, inv_mmin, inv_mmax, mcl
    real(kind=8) :: RandNum

    ! Lower cap: raises floor at high Mj, privileges massive clusters
    mmin_eff = mmin0
    if(Mj > 0.0d0 .and. Mj_ref > 0.0d0) &
      mmin_eff = mmin0 * (Mj / Mj_ref)**delta
    mmin_eff = max(mmin_eff, mmin0)

    ! Upper cap
    mmax_eff = mmax
    if(f_cap * Mj < mmax_eff) mmax_eff = f_cap * Mj

    n_cl     = 0
    M_rem    = M_sf
    inv_mmin = 1.0d0 / mmin_eff

    do while(M_rem >= mmin_eff .and. n_cl < max_clusters_per_event)
      mmax_draw = min(mmax_eff, M_rem)
      if(mmax_draw <= mmin_eff) exit

      ! Inverse CDF for dN/dM ~ M^{-2}:
      ! m = 1 / (1/mmin - u*(1/mmin - 1/mmax))
      inv_mmax = 1.0d0 / mmax_draw
      call ranf(seed, RandNum)
      mcl = 1.0d0 / (inv_mmin - dble(RandNum) * (inv_mmin - inv_mmax))

      n_cl = n_cl + 1
      cluster_masses(n_cl) = mcl
      M_rem = M_rem - mcl
    end do

    ! Exact mass conservation: remainder into last cluster
    if(n_cl > 0) then
      cluster_masses(n_cl) = cluster_masses(n_cl) + M_rem
    else
      ! M_sf below mmin_eff: single cluster with full mass
      n_cl = 1
      cluster_masses(1) = M_sf
    end if

  end subroutine sample_cmf_clusters

  !------------------------------------------------------------------
  ! Decompose a SF event into clusters drawn from a lognormal kernel
  ! tied to the unresolved turbulent density PDF (subgrid fragmentation
  ! model): ln M_cl ~ N(ln Mcl_char, sigma_lnM^2).
  ! Mcl_char and sigma_lnM are precomputed per SF event (see starform2,
  ! sf_model=1 branch) from the same multi-ff sigs/scrit already used
  ! for the star formation rate.
  ! Mass is exactly conserved: remainder added to last cluster.
  ! mmin/mmax bound both the sampling center Mcl_char and every draw,
  ! matching the caps already enforced by sample_cmf_clusters (and the
  ! range the eta_sn_cluster IMF table is built over).
  ! If M_event < mmin, returns one cluster with the full mass.
  ! All masses in Msun.
  !------------------------------------------------------------------
  subroutine sample_lognormal_clusters(M_event, Mcl_char, sigma_lnM, mmin, mmax, &
                                        seed, cluster_masses, n_cl)
    use random
    real(dp), intent(in)    :: M_event, Mcl_char, sigma_lnM, mmin, mmax
    integer,  intent(inout) :: seed(IRandNumSize)
    real(dp), intent(out)   :: cluster_masses(max_clusters_per_event)
    integer,  intent(out)   :: n_cl

    real(dp)     :: M_rem, mcl, lnMcl_char
    real(kind=8) :: GaussNum

    lnMcl_char = log(min(max(Mcl_char, mmin), mmax))

    n_cl  = 0
    M_rem = M_event

    do while(M_rem >= mmin .and. n_cl < max_clusters_per_event)
      call gaussdev(seed, GaussNum)
      mcl = exp(lnMcl_char + sigma_lnM*dble(GaussNum))
      mcl = min(max(mcl, mmin), mmax)
      mcl = min(mcl, M_rem)

      n_cl = n_cl + 1
      cluster_masses(n_cl) = mcl
      M_rem = M_rem - mcl
    end do

    ! Exact mass conservation: remainder into last cluster
    if(n_cl > 0) then
      cluster_masses(n_cl) = cluster_masses(n_cl) + M_rem
    else
      ! M_event below mmin: single cluster with full mass
      n_cl = 1
      cluster_masses(1) = M_event
    end if

  end subroutine sample_lognormal_clusters

end module imf_commons
