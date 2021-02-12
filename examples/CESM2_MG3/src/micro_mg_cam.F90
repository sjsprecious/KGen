module micro_mg_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for MG microphysics
!
!---------------------------------------------------------------------------------
!
! How to add new packed MG inputs to micro_mg_cam_tend:
!
! If you have an input with first dimension [psetcols, pver], the procedure
! for adding inputs is as follows:
!
! 1) In addition to any variables you need to declare for the "unpacked"
!    (CAM format) version, you must declare an array for the "packed" 
!    (MG format) version.
!
! 2) Add a call similar to the following line (look before the
!    micro_mg_tend calls to see similar lines):
!
!      packed_array = packer%pack(original_array)
!
!    The packed array can then be passed into any of the MG schemes.
!
! This same procedure will also work for 1D arrays of size psetcols, 3-D
! arrays with psetcols and pver as the first dimensions, and for arrays of
! dimension [psetcols, pverp]. You only have to modify the allocation of
! the packed array before the "pack" call.
!
!---------------------------------------------------------------------------------
!
! How to add new packed MG outputs to micro_mg_cam_tend:
!
! 1) As with inputs, in addition to the unpacked outputs you must declare
!    an array for packed data. The unpacked and packed arrays must *also* 
!    be targets or pointers (but cannot be both).
!
! 2) Add the field to post-processing as in the following line (again,
!    there are many examples before the micro_mg_tend calls):
!
!      call post_proc%add_field(p(final_array),p(packed_array))
!  
!    *** IMPORTANT ** If the fields are only being passed to a certain version of
!    MG, you must only add them if that version is being called (see
!    the "if (micro_mg_version >1)" sections below
!
!    This registers the field for post-MG averaging, and to scatter to the
!    final, unpacked version of the array.
!
!    By default, any columns/levels that are not operated on by MG will be
!    set to 0 on output; this value can be adjusted using the "fillvalue"
!    optional argument to post_proc%add_field.
!
!    Also by default, outputs from multiple substeps will be averaged after
!    MG's substepping is complete. Passing the optional argument
!    "accum_method=accum_null" will change this behavior so that the last
!    substep is always output.
!
! This procedure works on 1-D and 2-D outputs. Note that the final,
! unpacked arrays are not set until the call to
! "post_proc%process_and_unpack", which sets every single field that was
! added with post_proc%add_field.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8

use netcdf

implicit none
private
save

public  :: micro_mg_cam_tend_pack

integer :: num_steps = 3 ! Number of MG substeps

!===============================================================================
contains
!===============================================================================

subroutine micro_mg_cam_tend_pack(dtime, mgncol, nlev, ncid)

   use micro_mg_utils, only: size_dist_param_basic_vect, size_dist_param_liq_vect, &
        mg_liq_props, mg_ice_props, avg_diameter, rhoi, rhosn, rhow, rhows, &
        mg_graupel_props, rhog, &
        qsmall, mincld

   use micro_mg3_0, only: micro_mg_tend3_0 => micro_mg_tend

   real(r8), intent(in) :: dtime
   integer,  intent(in) :: nlev
   integer,  intent(in) :: mgncol
   integer,  intent(in) :: ncid

   ! Local variables
   integer :: lchnk, ncol, psetcols, ngrdcol

   integer :: i, k, itim_old, it

   ! Packed versions of inputs.
   real(r8) :: packed_t(mgncol,nlev)
   real(r8) :: packed_q(mgncol,nlev)
   real(r8) :: packed_qc(mgncol,nlev)
   real(r8) :: packed_nc(mgncol,nlev)
   real(r8) :: packed_qi(mgncol,nlev)
   real(r8) :: packed_ni(mgncol,nlev)
   real(r8) :: packed_qr(mgncol,nlev)
   real(r8) :: packed_nr(mgncol,nlev)
   real(r8) :: packed_qs(mgncol,nlev)
   real(r8) :: packed_ns(mgncol,nlev)
   real(r8) :: packed_qg(mgncol,nlev)
   real(r8) :: packed_ng(mgncol,nlev)   

   real(r8) :: packed_relvar(mgncol,nlev)
   real(r8) :: packed_accre_enhan(mgncol,nlev)

   real(r8) :: packed_p(mgncol,nlev)
   real(r8) :: packed_pdel(mgncol,nlev)

   real(r8) :: packed_cldn(mgncol,nlev)
   real(r8) :: packed_liqcldf(mgncol,nlev)
   real(r8) :: packed_icecldf(mgncol,nlev)
   real(r8) :: packed_qsatfac(mgncol,nlev)

   real(r8) :: packed_naai(mgncol,nlev)
   real(r8) :: packed_npccn(mgncol,nlev)

   real(r8) :: packed_rndst(mgncol,nlev,4)
   real(r8) :: packed_nacon(mgncol,nlev,4)

   ! Optional outputs.
   real(r8) :: packed_tnd_qsnow(mgncol,nlev)
   real(r8) :: packed_tnd_nsnow(mgncol,nlev)
   real(r8) :: packed_re_ice(mgncol,nlev)

   real(r8) :: packed_frzimm(mgncol,nlev)
   real(r8) :: packed_frzcnt(mgncol,nlev)
   real(r8) :: packed_frzdep(mgncol,nlev)

   ! Packed versions of outputs.
   real(r8), target :: packed_rate1ord_cw2pr_st(mgncol,nlev)
   real(r8), target :: packed_tlat(mgncol,nlev)
   real(r8), target :: packed_qvlat(mgncol,nlev)
   real(r8), target :: packed_qctend(mgncol,nlev)
   real(r8), target :: packed_qitend(mgncol,nlev)
   real(r8), target :: packed_nctend(mgncol,nlev)
   real(r8), target :: packed_nitend(mgncol,nlev)

   real(r8), target :: packed_qrtend(mgncol,nlev)
   real(r8), target :: packed_qstend(mgncol,nlev)
   real(r8), target :: packed_nrtend(mgncol,nlev)
   real(r8), target :: packed_nstend(mgncol,nlev)
   real(r8), target :: packed_qgtend(mgncol,nlev)
   real(r8), target :: packed_ngtend(mgncol,nlev)

   real(r8), target :: packed_prect(mgncol)
   real(r8), target :: packed_preci(mgncol)
   real(r8), target :: packed_nevapr(mgncol,nlev)
   real(r8), target :: packed_am_evp_st(mgncol,nlev)
   real(r8), target :: packed_evapsnow(mgncol,nlev)
   real(r8), target :: packed_prain(mgncol,nlev)
   real(r8), target :: packed_prodsnow(mgncol,nlev)
   real(r8), target :: packed_cmeout(mgncol,nlev)
   real(r8), target :: packed_qsout(mgncol,nlev)
   real(r8), target :: packed_cflx(mgncol,nlev+1)
   real(r8), target :: packed_iflx(mgncol,nlev+1)
   real(r8), target :: packed_rflx(mgncol,nlev+1)
   real(r8), target :: packed_sflx(mgncol,nlev+1)
   real(r8), target :: packed_gflx(mgncol,nlev+1)
   real(r8), target :: packed_qrout(mgncol,nlev)
   real(r8), target :: packed_qcsevap(mgncol,nlev)
   real(r8), target :: packed_qisevap(mgncol,nlev)
   real(r8), target :: packed_qvres(mgncol,nlev)
   real(r8), target :: packed_cmei(mgncol,nlev)
   real(r8), target :: packed_vtrmc(mgncol,nlev)
   real(r8), target :: packed_vtrmi(mgncol,nlev)
   real(r8), target :: packed_qcsedten(mgncol,nlev)
   real(r8), target :: packed_qisedten(mgncol,nlev)
   real(r8), target :: packed_qrsedten(mgncol,nlev)
   real(r8), target :: packed_qssedten(mgncol,nlev)
   real(r8), target :: packed_qgsedten(mgncol,nlev)
   real(r8), target :: packed_umg(mgncol,nlev)
   real(r8), target :: packed_umr(mgncol,nlev)
   real(r8), target :: packed_ums(mgncol,nlev)
   real(r8), target :: packed_pra(mgncol,nlev)
   real(r8), target :: packed_prc(mgncol,nlev)
   real(r8), target :: packed_mnuccc(mgncol,nlev)
   real(r8), target :: packed_mnucct(mgncol,nlev)
   real(r8), target :: packed_msacwi(mgncol,nlev)
   real(r8), target :: packed_psacws(mgncol,nlev)
   real(r8), target :: packed_bergs(mgncol,nlev)
   real(r8), target :: packed_berg(mgncol,nlev)
   real(r8), target :: packed_melt(mgncol,nlev)
   real(r8), target :: packed_homo(mgncol,nlev)
   real(r8), target :: packed_qcres(mgncol,nlev)
   real(r8), target :: packed_prci(mgncol,nlev)
   real(r8), target :: packed_prai(mgncol,nlev)
   real(r8), target :: packed_qires(mgncol,nlev)
   real(r8), target :: packed_mnuccr(mgncol,nlev)
   real(r8), target :: packed_mnuccri(mgncol,nlev)
   real(r8), target :: packed_pracs(mgncol,nlev)
   real(r8), target :: packed_meltsdt(mgncol,nlev)
   real(r8), target :: packed_frzrdt(mgncol,nlev)
   real(r8), target :: packed_mnuccd(mgncol,nlev)
   real(r8), target :: packed_nrout(mgncol,nlev)
   real(r8), target :: packed_nsout(mgncol,nlev)
   real(r8), target :: packed_refl(mgncol,nlev)
   real(r8), target :: packed_arefl(mgncol,nlev)
   real(r8), target :: packed_areflz(mgncol,nlev)
   real(r8), target :: packed_frefl(mgncol,nlev)
   real(r8), target :: packed_csrfl(mgncol,nlev)
   real(r8), target :: packed_acsrfl(mgncol,nlev)
   real(r8), target :: packed_fcsrfl(mgncol,nlev)
   real(r8), target :: packed_rercld(mgncol,nlev)
   real(r8), target :: packed_ncai(mgncol,nlev)
   real(r8), target :: packed_ncal(mgncol,nlev)
   real(r8), target :: packed_qrout2(mgncol,nlev)
   real(r8), target :: packed_qsout2(mgncol,nlev)
   real(r8), target :: packed_nrout2(mgncol,nlev)
   real(r8), target :: packed_nsout2(mgncol,nlev)
   real(r8), target :: packed_freqs(mgncol,nlev)
   real(r8), target :: packed_freqr(mgncol,nlev)
   real(r8), target :: packed_freqg(mgncol,nlev)
   real(r8), target :: packed_nfice(mgncol,nlev)
   real(r8), target :: packed_prer_evap(mgncol,nlev)
   real(r8), target :: packed_qcrat(mgncol,nlev)

   real(r8), target :: packed_rel(mgncol,nlev)
   real(r8), target :: packed_rei(mgncol,nlev)
   real(r8), target :: packed_sadice(mgncol,nlev)
   real(r8), target :: packed_sadsnow(mgncol,nlev)
   real(r8), target :: packed_lambdac(mgncol,nlev)
   real(r8), target :: packed_mu(mgncol,nlev)
   real(r8), target :: packed_des(mgncol,nlev)
   real(r8), target :: packed_dei(mgncol,nlev)

!Hail/Graupel Output
   real(r8), target :: packed_qgout(mgncol,nlev)   
   real(r8), target :: packed_ngout(mgncol,nlev)   
   real(r8), target :: packed_dgout(mgncol,nlev)                  
   real(r8), target :: packed_qgout2(mgncol,nlev) 
   real(r8), target :: packed_ngout2(mgncol,nlev) 
   real(r8), target :: packed_dgout2(mgncol,nlev) 
!Hail/Graupel Process Rates                
   real(r8), target :: packed_psacr(mgncol,nlev)   
   real(r8), target :: packed_pracg(mgncol,nlev)   
   real(r8), target :: packed_psacwg(mgncol,nlev)  
   real(r8), target :: packed_pgsacw(mgncol,nlev)
   real(r8), target :: packed_pgracs(mgncol,nlev) 
   real(r8), target :: packed_prdg(mgncol,nlev)   
   real(r8), target :: packed_qmultg(mgncol,nlev)  
   real(r8), target :: packed_qmultrg(mgncol,nlev)   
   real(r8), target :: packed_npracg(mgncol,nlev)
   real(r8), target :: packed_nscng(mgncol,nlev)
   real(r8), target :: packed_ngracs(mgncol,nlev)
   real(r8), target :: packed_nmultg(mgncol,nlev)
   real(r8), target :: packed_nmultrg(mgncol,nlev)
   real(r8), target :: packed_npsacwg(mgncol,nlev)

   ! Dummy arrays for cases where we throw away the MG version and
   ! recalculate sizes on the CAM grid to avoid time/subcolumn averaging
   ! issues.
   real(r8) :: rel_fn_dum(mgncol,nlev)
   real(r8) :: dsout2_dum(mgncol,nlev)
   real(r8) :: drout_dum(mgncol,nlev)
   real(r8) :: reff_rain_dum(mgncol,nlev)
   real(r8) :: reff_snow_dum(mgncol,nlev)
   real(r8) :: reff_grau_dum(mgncol,nlev)   !not used for now or passed to COSP.

   character(128) :: errstring   ! return status (non-blank for error return)

   ! For rrtmg optics. specified distribution.
   real(r8), parameter :: dcon   = 25.e-6_r8         ! Convective size distribution effective radius (meters)
   real(r8), parameter :: mucon  = 5.3_r8            ! Convective size distribution shape parameter
   real(r8), parameter :: deicon = 50._r8            ! Convective ice effective diameter (meters)

   real(r8), pointer :: pckdptr(:,:)

   ! This is for netCDF input only
   integer               :: varid
!   integer               :: len
!   integer, dimension(2) :: dimids

   ! read in variables from the netcdf file
   call check( nf90_inq_varid(ncid, "packed_t", varid) )
!   call check( nf90_inquire_variable(ncid, varid, dimids = dimids) )
!   call check( nf90_inquire_dimension(ncid, dimids(2), len=len) )
   call check( nf90_get_var(ncid, varid, packed_t, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_q", varid) )
   call check( nf90_get_var(ncid, varid, packed_q, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_qc", varid) )
   call check( nf90_get_var(ncid, varid, packed_qc, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_qi", varid) )
   call check( nf90_get_var(ncid, varid, packed_qi, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_nc", varid) )
   call check( nf90_get_var(ncid, varid, packed_nc, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_ni", varid) )
   call check( nf90_get_var(ncid, varid, packed_ni, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_qr", varid) )
   call check( nf90_get_var(ncid, varid, packed_qr, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_qs", varid) )
   call check( nf90_get_var(ncid, varid, packed_qs, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_nr", varid) )
   call check( nf90_get_var(ncid, varid, packed_nr, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_ns", varid) )
   call check( nf90_get_var(ncid, varid, packed_ns, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_qg", varid) )
   call check( nf90_get_var(ncid, varid, packed_qg, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_ng", varid) )
   call check( nf90_get_var(ncid, varid, packed_ng, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_relvar", varid) )
   call check( nf90_get_var(ncid, varid, packed_relvar, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_accre_enhan", varid) )
   call check( nf90_get_var(ncid, varid, packed_accre_enhan, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_p", varid) )
   call check( nf90_get_var(ncid, varid, packed_p, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_pdel", varid) )
   call check( nf90_get_var(ncid, varid, packed_pdel, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_cldn", varid) )
   call check( nf90_get_var(ncid, varid, packed_cldn, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_liqcldf", varid) )
   call check( nf90_get_var(ncid, varid, packed_liqcldf, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_icecldf", varid) )
   call check( nf90_get_var(ncid, varid, packed_icecldf, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_qsatfac", varid) )
   call check( nf90_get_var(ncid, varid, packed_qsatfac, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_naai", varid) )
   call check( nf90_get_var(ncid, varid, packed_naai, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_npccn", varid) )
   call check( nf90_get_var(ncid, varid, packed_npccn, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_tnd_qsnow", varid) )
   call check( nf90_get_var(ncid, varid, packed_tnd_qsnow, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_tnd_nsnow", varid) )
   call check( nf90_get_var(ncid, varid, packed_tnd_nsnow, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_re_ice", varid) )
   call check( nf90_get_var(ncid, varid, packed_re_ice, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_frzimm", varid) )
   call check( nf90_get_var(ncid, varid, packed_frzimm, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_frzcnt", varid) )
   call check( nf90_get_var(ncid, varid, packed_frzcnt, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_frzdep", varid) )
   call check( nf90_get_var(ncid, varid, packed_frzdep, &
                     start=(/1, 1/), count=(/mgncol,nlev/)) )

   call check( nf90_inq_varid(ncid, "packed_rndst", varid) )
   call check( nf90_get_var(ncid, varid, packed_rndst, &
                     start=(/1, 1, 1/), count=(/mgncol,nlev,4/)) )

   call check( nf90_inq_varid(ncid, "packed_nacon", varid) )
   call check( nf90_get_var(ncid, varid, packed_nacon, &
                     start=(/1, 1, 1/), count=(/mgncol,nlev,4/)) )

   !-------------------------------------------------------------------------------
         !$kgen begin_callsite mg3_tend 
         call micro_mg_tend3_0( &
              mgncol,         nlev,           dtime/num_steps,&
              packed_t,               packed_q,               &
              packed_qc,              packed_qi,              &
              packed_nc,              packed_ni,              &
              packed_qr,              packed_qs,              &
              packed_nr,              packed_ns,              &
              packed_qg,              packed_ng,              &
              packed_relvar,          packed_accre_enhan,     &
              packed_p,               packed_pdel,            &
              packed_cldn, packed_liqcldf, packed_icecldf, packed_qsatfac, &
              packed_rate1ord_cw2pr_st,                       &
              packed_naai,            packed_npccn,           &
              packed_rndst,           packed_nacon,           &
              packed_tlat,            packed_qvlat,           &
              packed_qctend,          packed_qitend,          &
              packed_nctend,          packed_nitend,          &
              packed_qrtend,          packed_qstend,          &
              packed_nrtend,          packed_nstend,          &
              packed_qgtend,          packed_ngtend,          &
              packed_rel,     rel_fn_dum,     packed_rei,     &
              packed_sadice,          packed_sadsnow,         &
              packed_prect,           packed_preci,           &
              packed_nevapr,          packed_evapsnow,        &
              packed_am_evp_st,                               &
              packed_prain,           packed_prodsnow,        &
              packed_cmeout,          packed_dei,             &
              packed_mu,              packed_lambdac,         &
              packed_qsout,           packed_des,             &
              packed_qgout,   packed_ngout,   packed_dgout,   &
              packed_cflx,    packed_iflx,                    &
              packed_gflx,                                    &
              packed_rflx,    packed_sflx,    packed_qrout,   &
              reff_rain_dum,          reff_snow_dum,   reff_grau_dum,       &
              packed_qcsevap, packed_qisevap, packed_qvres,   &
              packed_cmei,    packed_vtrmc,   packed_vtrmi,   &
              packed_umr,             packed_ums,             &
              packed_umg,             packed_qgsedten,        &
              packed_qcsedten,        packed_qisedten,        &
              packed_qrsedten,        packed_qssedten,        &
              packed_pra,             packed_prc,             &
              packed_mnuccc,  packed_mnucct,  packed_msacwi,  &
              packed_psacws,  packed_bergs,   packed_berg,    &
              packed_melt,            packed_homo,            &
              packed_qcres,   packed_prci,    packed_prai,    &
              packed_qires,   packed_mnuccr,  packed_mnuccri, packed_pracs,   &
              packed_meltsdt, packed_frzrdt,  packed_mnuccd,  &
              packed_pracg,   packed_psacwg,  packed_pgsacw,  &
              packed_pgracs,  packed_prdg,   &
              packed_qmultg,  packed_qmultrg, packed_psacr,   &
              packed_npracg,  packed_nscng,   packed_ngracs,  &
              packed_nmultg,  packed_nmultrg, packed_npsacwg, & 
              packed_nrout,           packed_nsout,           &
              packed_refl,    packed_arefl,   packed_areflz,  &
              packed_frefl,   packed_csrfl,   packed_acsrfl,  &
              packed_fcsrfl,          packed_rercld,          &
              packed_ncai,            packed_ncal,            &
              packed_qrout2,          packed_qsout2,          &
              packed_nrout2,          packed_nsout2,          &
              drout_dum,              dsout2_dum,             &
              packed_qgout2, packed_ngout2, packed_dgout2, packed_freqg,   &
              packed_freqs,           packed_freqr,           &
              packed_nfice,           packed_qcrat,           &
              errstring, &
              packed_tnd_qsnow,packed_tnd_nsnow,packed_re_ice,&
              packed_prer_evap,                                     &
              packed_frzimm,  packed_frzcnt,  packed_frzdep   )
         !$kgen end_callsite mg3_tend 

end subroutine micro_mg_cam_tend_pack

subroutine check(status)
   
   integer, intent(in) :: status
   
   if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
   end if

end subroutine check

end module micro_mg_cam
