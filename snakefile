# Snakemake
##
## @Kevin Quinteros
##

#check minimum snakemake version utility 
from snakemake.utils import min_version

#####set minimum snakemakeversion####
min_version("5.30.1")

##--- Importing Configuration Files ---##
configfile: '01_config/config.yaml' 
        
##--- include rules ---##
include: '03_rules/common.smk' #contains input/output and helper functions. additional output arrays and libraries defined
include: '03_rules/validate.smk' #contains the main rules for bitacora pipeline
 
##--- universal rule that checks the output of every rule ---##
rule all:
    input:
        run_filter_overlapping_output
