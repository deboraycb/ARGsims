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

### Other general scripts

- `calc_kld_ranks.R` calculate KLD from ranks distributions
- `plot_ranks.R` plot distribution of ranks
- `organizeQuantileOutput` get counts of ranks

### scripts/argweaver

- `2_vcf2sites.py` convert vcf to ARGweaver sites input format
- `3_run_argweaver_smcprime_hudsim.sh` run ARGweaver
- `4_arg-summarize2plot_allchunks.sh` summarize args to extract pw coalescence times and call plotting scripts:
    
    - `coalTimeRank.py` get ranks of coalescence times
    - `getPostMean_bedmap.sh` get point estimates of coalescence times (mean of the posterior) and match to true (simulated) values
    - `MSPcoalTime2bed` convert coalescence times from msprime to bed format
    - `plot_pointestimates.R` plot heatmaps of point estimates (estimated vs. simulated)
    - `plotTc_argw.R` plot distribution of coalescence times

- `5_pointestimates.sh` plot point estimates of coalescence times
    
### scripts/relate

- `2_run_relate2plots.sh` run Relate, extract pw coalescence times and call all plotting scripts:
    
    - `outputCoalTime` output pairwise coalescence times from anc/mut files
    - `plot_pointest.R` plot point estimates
    - `plotTc_relate.R` plot distribution of coalescence times
    - `thin.R` thin MCMC samples

- recombination maps: `sim_std.map`, `sim_r2e7.map`, `sim_r2e9.map`

### scripts/tsinferdate

- `2_tsinfer_date.py` run tsinfer+tsdate
- `3_organizeTc.sh` match simulated to estimated coal times
- `4_plot_pointestimates_tsinfer.R` plot point estimates of coalescence times
- `5_plotTc_tsdate.R` plot distribution of all pairwise coalescence times
