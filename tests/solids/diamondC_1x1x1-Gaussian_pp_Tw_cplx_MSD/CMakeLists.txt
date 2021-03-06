
IF(ENABLE_SOA)
 IF (NOT QMC_CUDA)
   IF (QMC_COMPLEX)

# VMC long run with no Jastrows 
# Compare directly to sCI Energy of: E=-11.65519310
#Data generated using utils/make_ref_data.sh  diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG qmc_ref/qmc_long.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat 
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "totenergy" "-11.65519310     0.006077601667")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "variance" "1.539759784521 0.027855789336")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "kinetic" "7.857471129807 0.019397520629")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "potential" "-19.513384734772 0.020252648227")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "eeenergy" "-2.835644692415 0.007764684716")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "localecp" "-5.360598387703 0.029116710254")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "nonlocalecp" "1.458525910786 0.010167449064")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "ionion" "-12.775667565000 0.100000000000")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD"
                    vmc_long  
                    vmc_long.in.xml
                    1 16
                    TRUE
                    0 diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS #VMC 
                    )

# VMC short run with no Jastrows 
# Compare directly to sCI Energy of: E=-11.65519310
#Data generated using utils/make_ref_data.sh  diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT qmc_ref/qmc_short.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat 
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "totenergy" "-11.65519310     0.011805237730")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "variance" "1.539759784521 0.211329760575")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "kinetic" "7.857471129807 0.053834125981")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "potential" "-19.513384734772 0.054910387369")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "eeenergy" "-2.835644692415 0.017280154413")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "localecp" "-5.360598387703 0.096994263787")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "nonlocalecp" "1.458525910786 0.029846773207")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "ionion" "-12.775667565000 0.100000000000")
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD"
                    vmc_short 
                    vmc_short.in.xml
                    1 16
                    TRUE
                    0 diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS #VMC 
                    )




    ELSE()
      MESSAGE_VERBOSE("Skipping Complex Periodic LCAO  if Complex code not build (QMC_COMPLEX=0)")
    ENDIF()
  ELSE()
    MESSAGE_VERBOSE("Skipping Periodic LCAO as not supported by CUDA build (QMC_CUDA=1)")
  ENDIF()
ELSE()
  MESSAGE_VERBOSE("Skipping Periodic LCAO as not supported by AoS build (ENABLE_SOA=0)")
ENDIF()
