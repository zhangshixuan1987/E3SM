#!/bin/bash -fe

# E3SM Water Cycle v2.1 run_e3sm script template.
#
# Bash coding style inspired by:
# http://kfirlavi.herokuapp.com/blog/2012/11/14/defensive-bash-programming

main() {

# For debugging, uncomment line below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=pm-cpu
readonly PROJECT="m4506" 
readonly compiler="gnu"

# Simulation
readonly COMPSET="20TRSOI_EAM%CMIP6_ELM%SPBC_MPASSI_MPASO_MOSART_SGLC_WW3"
readonly RESOLUTION="atlanticx4v1pg2_r05_EC30to60E2r2_wQU225EC30to60E2r2"
# Code and compilation
readonly CHECKOUT="20240130"
readonly BRANCH="8042a73973e73ecf48fae969a8ae3cf434b1b8c0"
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# BEFORE RUNNING : CHANGE the following CASE_NAME to desired value
readonly CASE_NAME="v2_1.WCYCL20TR-WW3.${RESOLUTION}"

# Run options
readonly MODEL_START_TYPE="initial"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="1978-01-01"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE
readonly RUN_REFCASE="v2_1.LR.historical_0301"
readonly RUN_REFDATE="1978-01-01"   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'
readonly RUN_REFTOD="00000"
readonly RUN_REFDIR="/pscratch/sd/z/zhan391/WW3/initial_condition/${RUN_REFDATE}-${RUN_REFTOD}"

#readonly atm_in="${RUN_REFDIR}/v2_1.LR.historical.nco.ne32x4_ATL.eam.i.1978-01-01-00000.nc"
readonly atm_in="${RUN_REFDIR}/v2_1.LR.historical_0301.ne32x4_ATL.eam.i.1978-01-01-00000.nc"
readonly lnd_in="${RUN_REFDIR}/v2.LR.BGC-LNDATM.CONTRL.elm.r.1980-01-01-00000.nc"
readonly lnd_surdat="/global/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr1850_c211019.nc"
readonly lnd_use="/global/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/landuse.timeseries_0.5x0.5_HIST_simyr1850-2015_c211019.nc"
readonly lnd_soidat="/global/cfs/cdirs/e3sm/inputdata/lnd/clm2/paramdata/CNP_parameters_c180529.nc"
readonly drydep_in="/pscratch/sd/z/zhan391/WW3/ne32x4_mesh3/meshes/atmsrf_ne32x4_ATLpg2_20240531.nc"
readonly topo_in="/pscratch/sd/z/zhan391/WW3/ne32x4_mesh3/meshes/USGS-gtopo30_ne32x4_ATLpg2_12xdel2.nc"

readonly OCN_INIT_FILE="${RUN_REFDIR}/v2_1.LR.historical_0301.mpaso.rst.${RUN_REFDATE}_${RUN_REFTOD}.nc"
readonly ICE_INIT_FILE="${RUN_REFDIR}/v2_1.LR.historical_0301.mpassi.rst.${RUN_REFDATE}_${RUN_REFTOD}.nc"
readonly RIV_INIT_FILE="${RUN_REFDIR}/v2.LR.BGC-LNDATM.CONTRL.mosart.r.1980-01-01-00000.nc"

echo $atm_in 
echo $lnd_in
echo $OCN_INIT_FILE
echo $ICE_INIT_FILE
echo $RIV_INIT_FILE

# Set paths
readonly CODE_ROOT="/pscratch/sd/z/zhan391/WW3/E3SM"
readonly CASE_ROOT="/pscratch/sd/z/zhan391/WW3/simulation/${CASE_NAME}"
# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'XS_2x5_ndays', 'XS_1x10_ndays', 'S_1x10_ndays',
#               'M_1x10_ndays', 'ML_1x10_ndays', 'L_1x10_ndays'
#  or 'production' for full simulation
readonly run='production'
#readonly run='custom_20x1_nmonths'
#readonly run='custom_20x1_production'

if [[ "${run}" != *"production"* ]] ;then
  # Short test simulations
  tmp=($(echo $run | tr "_" " "))
  layout=${tmp[0]}
  units=${tmp[2]}
  resubmit=$(( ${tmp[1]%%x*} -1 ))
  length=${tmp[1]##*x}
  readonly custom_nodes=20
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/tests/${run}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/tests/${run}/run
  readonly PELAYOUT=${layout}
  readonly WALLTIME="01:00:00"
  readonly STOP_OPTION=${units}
  readonly STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=false
else
  # Production simulation
  readonly custom_nodes=20 #64
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT="L"
  readonly WALLTIME="48:00:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="10"
  readonly REST_OPTION="nyears"
  readonly REST_N="1"
  readonly RESUBMIT="2"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history
readonly HIST_OPTION="nyears"
readonly HIST_N="1"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

# --- Toggle flags for what to do ----
do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true

# --- Now, do the work ---
if [ "${do_create_newcase,,}" == "true" ]; then
  if [ -d ${CASE_ROOT}/case_scripts ];then
    rm -rvf ${CASE_ROOT}/case_scripts
  fi
fi 

if [ "${do_case_build,,}" == "true" ]; then
  if [ -d ${CASE_BUILD_DIR} ];then
    rm -rvf ${CASE_BUILD_DIR}
  fi
fi

# Make directories created by this script world-readable
umask 022

# Fetch code from Github
fetch_code

# Create case
create_newcase

# Setup
case_setup

# Build
case_build

# Configure runtime options
runtime_options

# Copy script into case_script directory for provenance
copy_script

# Submit
case_submit

# All done
echo $'\n----- All done -----\n'

}

# =======================
# Custom user_nl settings
# =======================

user_nl() {

cat << EOF >> user_nl_eam
 dtime                 = 1800
 se_tstep              = 75  ! dt_dyn    = se_tstep = 75s  
 dt_remap_factor       = 2   ! dt_remap  = se_tstep*dt_remap_factor      = 150s
 dt_tracer_factor      = 6   ! dt_tracer = se_tstep*dt_tracer_factor     = 450s 
 hypervis_subcycle     = 2   ! dt_vis    = dt_dyn/hypervis_subcycle      = 75s
 hypervis_subcycle_q   = 6   ! dt_vis_q  = dt_tracer/hypervis_subcycle_q = 75s
 hypervis_subcycle_tom = 3

 nhtfrq =   0,-24,-6,-6,-3,-24,0
 mfilt  = 1,30,120,120,240,30,1
 ncdata = "${atm_in}"
 drydep_srf_file="${drydep_in}"
 bnd_topo= "${topo_in}"

 avgflag_pertape = 'A','A','I','A','A','A','I'
 fexcl1 = 'CFAD_SR532_CAL', 'LINOZ_DO3', 'LINOZ_DO3_PSC', 'LINOZ_O3CLIM', 'LINOZ_O3COL', 'LINOZ_SSO3', 'hstobie_linoz'
 fincl1 = 'extinct_sw_inp','extinct_lw_bnd7','extinct_lw_inp','CLD_CAL', 'TREFMNAV', 'TREFMXAV'
 fincl2 = 'FLUT','PRECT','U200','V200','U850','V850','Z500','OMEGA500','UBOT','VBOT','TREFHT','TREFHTMN:M','TREFHTMX:X','QREFHT','TS','PS','TMQ','TUQ','TVQ','TOZ', 'FLDS', 'FLNS', 'FSDS', 'FSNS', 'SHFLX', 'LHFLX', 'TGCLDCWP', 'TGCLDIWP', 'TGCLDLWP', 'CLDTOT', 'T250', 'T200', 'T150', 'T100', 'T050', 'T025', 'T010', 'T005', 'T002', 'T001', 'TTOP', 'U250', 'U150', 'U100', 'U050', 'U025', 'U010', 'U005', 'U002', 'U001', 'UTOP', 'FSNT', 'FLNT'
 fincl3 = 'PSL','T200','T500','U850','V850','UBOT','VBOT','TREFHT', 'Z700', 'TBOT:M'
 fincl4 = 'FLUT','U200','U850','PRECT','OMEGA500'
 fincl5 = 'PRECT','PRECC','TUQ','TVQ','QFLX','SHFLX','U90M','V90M'
 fincl6 = 'CLDTOT_ISCCP','MEANCLDALB_ISCCP','MEANTAU_ISCCP','MEANPTOP_ISCCP','MEANTB_ISCCP','CLDTOT_CAL','CLDTOT_CAL_LIQ','CLDTOT_CAL_ICE','CLDTOT_CAL_UN','CLDHGH_CAL','CLDHGH_CAL_LIQ','CLDHGH_CAL_ICE','CLDHGH_CAL_UN','CLDMED_CAL','CLDMED_CAL_LIQ','CLDMED_CAL_ICE','CLDMED_CAL_UN','CLDLOW_CAL','CLDLOW_CAL_LIQ','CLDLOW_CAL_ICE','CLDLOW_CAL_UN'
 fincl7 = 'O3', 'PS', 'TROP_P'
EOF

cat << EOF >> user_nl_elm
 hist_dov2xy = .true.,.true.
 hist_fincl2 = 'H2OSNO', 'FSNO', 'QRUNOFF', 'QSNOMELT', 'FSNO_EFF', 'SNORDSL', 'SNOW', 'FSDS', 'FSR', 'FLDS', 'FIRE', 'FIRA'
 hist_mfilt = 1,365
 hist_nhtfrq = 0,-24
 hist_avgflag_pertape = 'A','A'

 check_finidat_fsurdat_consistency = .false.
 check_finidat_year_consistency = .false.
 check_dynpft_consistency = .false.
 check_finidat_pct_consistency   = .false.

 do_budgets = .true.

 create_crop_landunit = .true.
 do_transient_crops = .false.
 do_transient_pfts = .true.

 finidat="${lnd_in}"
!Unless using the above finidat for 1850, also set the following, esp. the 2nd one
 fsurdat="${lnd_surdat}"
 fsoilordercon="${lnd_soidat}"
 flanduse_timeseries="${lnd_use}"
EOF

cat << EOF >> user_nl_mosart
 rtmhist_fincl2 = 'RIVER_DISCHARGE_OVER_LAND_LIQ'
 rtmhist_mfilt = 1,365
 rtmhist_ndens = 2
 rtmhist_nhtfrq = 0,-24
EOF

cat << EOF >> user_nl_mpassi
 config_initial_condition_type = 'restart'
 config_do_restart = .true.
 config_do_restart_bgc = .true.
 config_do_restart_hbrine = .true.
 config_do_restart_snow_density = .true.
 config_do_restart_snow_grain_radius = .true.
 config_do_restart_zsalinity = .true.
 config_restart_timestamp_name = 'rpointer.ice'
EOF

#cat << EOF >> user_nl_ww3
# fldcou = "USS USP HS FP DP"
# fldout = "WND HS FP DP USS"
# grd_out_freq = 0
# pnt_out_freq = 0
# stafile = '/lcrc/group/e3sm/data/inputdata/wav/ww3/stations.txt'
#EOF

}

patch_mpas_streams() {
echo
echo 'Modifying MPAS streams files'
pushd ../run
# change streams.ocean file
patch streams.ocean << EOF
--- streams.ocean
+++ streams.ocean
@@ -12,1 +12,1 @@
-                  filename_template="/global/cfs/cdirs/e3sm/inputdata/ocn/mpas-o/EC30to60E2r2/mpaso.EC30to60E2r2.rstFromG-anvil.201001.nc"
+                  filename_template="${OCN_INIT_FILE}"
EOF
# change streams.seaice file
patch streams.seaice << EOF
@@ -11,1 +11,1 @@
-                  filename_template="/global/cfs/cdirs/e3sm/inputdata/ice/mpas-cice/EC30to60E2r2/mpassi.EC30to60E2r2.rstFromG-anvil.201001.nc"
+                  filename_template="${ICE_INIT_FILE}"
@@ -38,1 +38,1 @@
-                  filename_template="/global/cfs/cdirs/e3sm/inputdata/ice/mpas-cice/EC30to60E2r2/mpassi.EC30to60E2r2.rstFromG-anvil.201001.nc"
+                  filename_template="${ICE_INIT_FILE}"
EOF

# copy to SourceMods
cp streams.ocean ../case_scripts/SourceMods/src.mpaso/
cp streams.seaice ../case_scripts/SourceMods/src.mpassi/

popd
}

######################################################
### Most users won't need to change anything below ###
######################################################

#-----------------------------------------------------
fetch_code() {

    if [ "${do_fetch_code,,}" != "true" ]; then
        echo $'\n----- Skipping fetch_code -----\n'
        return
    fi

    echo $'\n----- Starting fetch_code -----\n'
    local path=${CODE_ROOT}
    local repo=e3sm

    echo "Cloning $repo repository branch $BRANCH under $path"
    if [ -d "${path}" ]; then
        echo "ERROR: Directory already exists. Not overwriting"
        exit 20
    fi
    mkdir -p ${path}
    pushd ${path}

    # This will put repository, with all code
    git clone git@github.com:E3SM-Project/${repo}.git .

    # Setup git hooks
    rm -rf .git/hooks
    git clone git@github.com:E3SM-Project/E3SM-Hooks.git .git/hooks
    git config commit.template .git/hooks/commit.template

    # Check out desired branch
    git checkout ${BRANCH}

    # Custom addition
    if [ "${CHERRY}" != "" ]; then
        echo ----- WARNING: adding git cherry-pick -----
        for commit in "${CHERRY[@]}"
        do
            echo ${commit}
            git cherry-pick ${commit}
        done
        echo -------------------------------------------
    fi

    # Bring in all submodule components
    git submodule update --init --recursive

    popd
}

#-----------------------------------------------------
create_newcase() {

    if [ "${do_create_newcase,,}" != "true" ]; then
        echo $'\n----- Skipping create_newcase -----\n'
        return
    fi

    echo $'\n----- Starting create_newcase -----\n'

    # Base arguments
    args=" --case ${CASE_NAME} \
        --output-root ${CASE_ROOT} \
        --script-root ${CASE_SCRIPTS_DIR} \
        --handle-preexisting-dirs u \
        --compset ${COMPSET} \
        --res ${RESOLUTION} \
        --machine ${MACHINE} \
        --walltime ${WALLTIME} \
        --pecount ${PELAYOUT}"

    # Oprional arguments
    if [ ! -z "${PROJECT}" ]; then
      args="${args} --project ${PROJECT}"
    fi
    if [ ! -z "${CASE_GROUP}" ]; then
      args="${args} --case-group ${CASE_GROUP}"
    fi
    if [ ! -z "${QUEUE}" ]; then
      args="${args} --queue ${QUEUE}"
    fi

    if [ ! -z "${compiler}" ]; then
      args="${args} --compiler ${compiler}"
    fi

    ${CODE_ROOT}/cime/scripts/create_newcase ${args}

    if [ $? != 0 ]; then
      echo $'\nNote: if create_newcase failed because sub-directory already exists:'
      echo $'  * delete old case_script sub-directory'
      echo $'  * or set do_newcase=false\n'
      exit 35
    fi

  echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'

  # Number of cores per node (machine specific)
  if [ "${MACHINE}" == "chrysalis" ]; then
      ncore=64
  elif [ "${MACHINE}" == "compy" ]; then
      ncore=40
  elif [ "${MACHINE}" == "pm-cpu" ]; then
      ncore=64
  else
      echo 'ERROR: MACHINE = '${MACHINE}' is not supported for custom PE layout.'
      exit 400
  fi

  # Extract number of nodes
  nnodes=${custom_nodes}

  # Customize
  pushd ${CASE_SCRIPTS_DIR}
  ./xmlchange NTASKS=$(( $nnodes * $ncore ))
  ./xmlchange NTHRDS=1
  ./xmlchange ROOTPE=0
  ./xmlchange MAX_MPITASKS_PER_NODE=$ncore
  ./xmlchange MAX_TASKS_PER_NODE=$ncore
  echo $nnodes $ncore 
  echo $(( $nnodes * $ncore ))

  ### Added by Shixuan for tri-grid
  ./xmlchange CPL_NTASKS=$(( $nnodes * $ncore ))
  ./xmlchange ATM_NTASKS=$(( $nnodes * $ncore ))
  ./xmlchange OCN_NTASKS=$(( $nnodes * $ncore ))
  ./xmlchange LND_NTASKS=$(( $nnodes * $ncore ))
  ./xmlchange ROF_NTASKS=$(( $nnodes * $ncore ))
  ./xmlchange ICE_NTASKS=$(( $nnodes * $ncore ))

  popd

}

#-----------------------------------------------------
case_setup() {

    if [ "${do_case_setup,,}" != "true" ]; then
        echo $'\n----- Skipping case_setup -----\n'
        return
    fi

    echo $'\n----- Starting case_setup -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Setup some CIME directories
    ./xmlchange EXEROOT=${CASE_BUILD_DIR}
    ./xmlchange RUNDIR=${CASE_RUN_DIR}

    # Short term archiving
    ./xmlchange DOUT_S=${DO_SHORT_TERM_ARCHIVING^^}
    ./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}

    # Build with COSP, except for a data atmosphere (datm)
    if [ `./xmlquery --value COMP_ATM` == "datm"  ]; then
      echo $'\nThe specified configuration uses a data atmosphere, so cannot activate COSP simulator\n'
    else
      echo $'\nConfiguring E3SM to use the COSP simulator\n'
      ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
    fi

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    # Custom user_nl
    user_nl

    # Finally, run CIME case.setup
    ./case.setup --reset

    popd
}

#-----------------------------------------------------
case_build() {

    pushd ${CASE_SCRIPTS_DIR}

    # do_case_build = false
    if [ "${do_case_build,,}" != "true" ]; then

        echo $'\n----- case_build -----\n'

        if [ "${OLD_EXECUTABLE}" == "" ]; then
            # Ues previously built executable, make sure it exists
            if [ -x ${CASE_BUILD_DIR}/e3sm.exe ]; then
                echo 'Skipping build because $do_case_build = '${do_case_build}
            else
                echo 'ERROR: $do_case_build = '${do_case_build}' but no executable exists for this case.'
                exit 297
            fi
        else
            # If absolute pathname exists and is executable, reuse pre-exiting executable
            if [ -x ${OLD_EXECUTABLE} ]; then
                echo 'Using $OLD_EXECUTABLE = '${OLD_EXECUTABLE}
                cp -fp ${OLD_EXECUTABLE} ${CASE_BUILD_DIR}/
            else
                echo 'ERROR: $OLD_EXECUTABLE = '$OLD_EXECUTABLE' does not exist or is not an executable file.'
                exit 297
            fi
        fi
        echo 'WARNING: Setting BUILD_COMPLETE = TRUE.  This is a little risky, but trusting the user.'
        ./xmlchange BUILD_COMPLETE=TRUE

    # do_case_build = true
    else

        echo $'\n----- Starting case_build -----\n'

        # Turn on debug compilation option if requested
        if [ "${DEBUG_COMPILE^^}" == "TRUE" ]; then
            ./xmlchange DEBUG=${DEBUG_COMPILE^^}
        fi

        # Run CIME case.build
        ./case.build

        # Some user_nl settings won't be updated to *_in files under the run directory
        # Call preview_namelists to make sure *_in and user_nl files are consistent.
        ./preview_namelists

    fi

    popd
}

#-----------------------------------------------------
runtime_options() {

    echo $'\n----- Starting runtime_options -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Set simulation start date
    ./xmlchange RUN_STARTDATE=${START_DATE}

    # Segment length
    ./xmlchange STOP_OPTION=${STOP_OPTION,,},STOP_N=${STOP_N}

    # Restart frequency
    ./xmlchange REST_OPTION=${REST_OPTION,,},REST_N=${REST_N}

    # Coupler history
    ./xmlchange HIST_OPTION=${HIST_OPTION,,},HIST_N=${HIST_N}

    # Coupler budgets (always on)
    ./xmlchange BUDGETS=TRUE

    # Set resubmissions
    if (( RESUBMIT > 0 )); then
        ./xmlchange RESUBMIT=${RESUBMIT}
    fi

    # Run type
    # Start from default of user-specified initial conditions
    if [ "${MODEL_START_TYPE,,}" == "initial" ]; then
        ./xmlchange RUN_TYPE="startup"
        ./xmlchange CONTINUE_RUN="FALSE"

    # Continue existing run
    elif [ "${MODEL_START_TYPE,,}" == "continue" ]; then
        ./xmlchange CONTINUE_RUN="TRUE"

    elif [ "${MODEL_START_TYPE,,}" == "branch" ] || [ "${MODEL_START_TYPE,,}" == "hybrid" ]; then
        ./xmlchange RUN_TYPE=${MODEL_START_TYPE,,}
        ./xmlchange GET_REFCASE=${GET_REFCASE}
	./xmlchange RUN_REFDIR=${RUN_REFDIR}
        ./xmlchange RUN_REFCASE=${RUN_REFCASE}
        ./xmlchange RUN_REFDATE=${RUN_REFDATE}
        echo 'Warning: $MODEL_START_TYPE = '${MODEL_START_TYPE}
	echo '$RUN_REFDIR = '${RUN_REFDIR}
	echo '$RUN_REFCASE = '${RUN_REFCASE}
	echo '$RUN_REFDATE = '${START_DATE}

    else
        echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
        exit 380
    fi

    # Patch mpas streams files
    patch_mpas_streams

    popd
}

#-----------------------------------------------------
case_submit() {

    if [ "${do_case_submit,,}" != "true" ]; then
        echo $'\n----- Skipping case_submit -----\n'
        return
    fi

    echo $'\n----- Starting case_submit -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Run CIME case.submit
    ./case.submit

    popd
}

#-----------------------------------------------------
copy_script() {

    echo $'\n----- Saving run script for provenance -----\n'

    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    local this_script_name=`basename $0`
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    cp -vp ${this_script_name} ${script_provenance_dir}/${script_provenance_name}

}

#-----------------------------------------------------
# Silent versions of popd and pushd
pushd() {
    command pushd "$@" > /dev/null
}
popd() {
    command popd "$@" > /dev/null
}

# Now, actually run the script
#-----------------------------------------------------
main
