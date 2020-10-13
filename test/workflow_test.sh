svteaser -h

# SURVIVOR SIM TEST
rm -rf mito_test.svt/
svteaser surv_sim chrM.fa mito_test --num_sv_regions 2 --debug
svteaser sim_reads mito_test.svt

# KNOWN SV SIM TEST
rm -rf known_sv_sim_test.svt
svteaser known_sv_sim chrM.fa sample_known_sv_chrM.vcf known_sv_sim_test
svteaser sim_reads known_sv_sim_test.svt
