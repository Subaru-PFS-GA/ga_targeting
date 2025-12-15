SUBARU_PFS_ROOT="~/Subaru-PFS"
SUBARU_PFS_GA_ROOT="~/Subaru-PFS-GA"

export PFSSPEC_ROOT="$SUBARU_PFS_GA_ROOT/ga_pfsspec-all"
export PFSSPEC_DATA="~/data/pfsspec"

export CMDFIT_ROOT="$SUBARU_PFS_GA_ROOT/ga_cmdfit"
export CMDFIT_DATA="~/data/cmdfit"

export PFS_UTILS="$SUBARU_PFS_ROOT/pfs_utils/python"
export PFS_INSTDATA="$SUBARU_PFS_ROOT/pfs_instdata/python"
export PFS_ICS_COBRAOPS="$SUBARU_PFS_ROOT/ics_cobraOps/python"
export PFS_ICS_COBRA_CHARMER="$SUBARU_PFS_ROOT/ics_cobraCharmer/python"
export PFS_ICS_FPSACTOR="$SUBARU_PFS_ROOT/ics_fpsActor/python"
export PFS_ETS_FIBER_ASSIGNER="$SUBARU_PFS_ROOT/ets_fiberalloc"
export PFS_SPT_OPDB="$SUBARU_PFS_ROOT/spt_operational_database/python"
export PFS_ISOCHRONES="$SUBARU_PFS_GA_ROOT/ga_isochrones/python"
export PFS_GA_COMMON="$SUBARU_PFS_GA_ROOT/ga_common/python"

export PFS_TARGETING_DEBUGPORT=5678
export PFS_TARGETING_ROOT="$SUBARU_PFS_GA_ROOT/ga_targeting"
export PFS_TARGETING_DATA="$CMDFIT_DATA"
export PFS_TARGETING_TEMP="$SUBARU_PFS_GA_ROOT/ga_targeting/tmp"
export PFS_TARGETING_CONDAPATH="~/miniconda3"
export PFS_TARGETING_CONDAENV="ga-targeting"

export PFS_TARGETING_MODULES="pfs_utils:$PFS_UTILS:python
pfs_instdata:$PFS_INSTDATA:python
ics_cobraOps:$PFS_ICS_COBRAOPS:python
ics_cobraCharmer:$PFS_ICS_COBRA_CHARMER:python
ics_fpsActor:$PFS_ICS_FPSACTOR:python
ets_fiberalloc:$PFS_ETS_FIBER_ASSIGNER:.
spt_operational_database:$PFS_SPT_OPDB:python
datamodel:$PFS_DATAMODEL:python
ga_isochrones:$PFS_GA_ISOCHRONES:python
ga_common:$PFS_GA_COMMON:python
ga_targeting:$PFS_TARGETING_ROOT:python
ga_pfsspec:$PFSSPEC_ROOT:python"