rm -rf mito_test.svt/

svteaser -h

svteaser surv_sim chrM.fa mito_test --num_sv_regions 2 --debug
svteaser sim_reads mito_test.svt
