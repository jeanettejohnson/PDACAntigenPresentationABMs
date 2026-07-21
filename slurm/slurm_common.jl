# Generic SLURM helpers shared by every script in this directory -- no
# PhysiCellModelManager dependency, so this is safe to `include` from
# imc_spatial.jl (which doesn't load PCMM) as well as from hpc_setup.jl (which
# does, for the PCMM-based scripts).

# Determine the SLURM account to bill jobs to: pull this user's associated
# accounts from sacctmgr and use the first one listed. Set
# ENV["PCMM_SLURM_ACCOUNT"] to override (submit_driver.sh sets this so the
# account is computed once and passed through, rather than re-queried here).
function slurmAccount()
    haskey(ENV, "PCMM_SLURM_ACCOUNT") && return ENV["PCMM_SLURM_ACCOUNT"]
    username = readchomp(`whoami`)
    accounts = split(readchomp(`sacctmgr -n -P show assoc user=$username format=account`), "\n"; keepempty=false)
    isempty(accounts) && error("No SLURM account found for user $username. Set ENV[\"PCMM_SLURM_ACCOUNT\"] explicitly.")
    return accounts[1]
end

# Shared SLURM resource request for per-simulation / per-ROI jobs (not the
# driver job itself -- see submit_driver.sh for that one's own resources).
const JOB_CPUS_PER_TASK = "1"
const JOB_MEM = "2.5G"
const JOB_TIME = "7-00:00:00"
