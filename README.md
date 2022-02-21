# ARGsims

Code and simulated data used for the paper "Evaluation of methods for estimating coalescence times using ancestral recombination graphs"
https://www.biorxiv.org/content/10.1101/2021.11.15.468686v2.abstract

## trees
tree sequences simulated with msprime (only the ones within github file size limit)

## vcf
vcf files simulated with msprime (only the ones within github file size limit)

## scripts

### Simulations

Scripts for msprime simulations:

- `1_msprime_sim.py`
Simulations with msprime

- `2_msprimefinitesites.py`
Add mutations under JC finite sites model to previously simulated trees

### Analyses

- `calc_kld_ranks.R` calculate KLD from ranks distributions
- `plot_ranks.R` plot distribution of ranks
- `organizeQuantileOutput` get counts of ranks

### scripts/relate

- `2_run_relate2plots.sh` run Relate and call all plotting scripts:
    
    - `outputCoalTime` output pairwise coalescence times from anc/mut files
    - `thin.R` thin MCMC samples
    - `plot_pointest.R` plot point estimates
    - `plotTc_relate.R` plot distribution of coalescence times

- recombination maps: `sim_std.map`, `sim_r2e7.map`, `sim_r2e9.map`

### scripts/argweaver

- `2_vcf2sites.py` convert vcf to ARGweaver sites input format
- `3_run_argweaver_smcprime_hudsim.sh` run ARGweaver
- `4_arg-summarize2plot_allchunks.sh` summarize args and call plotting scripts:
    
    - `coalTimeRank.py` to get ranks of coalescence times
    - `MSPcoalTime2bed` convert coalescence times from msprime to bed format
    - `plotTc_argw.R` plot distribution of coalescence times

- `5_pointestimates.sh` plot point estimates of coalescence times
    
### scripts/tsinferdate

- `2_tsinfer_date.py` run tsinfer+tsdate
- `3_organizeTc.sh` match simulated to estimated coal times
