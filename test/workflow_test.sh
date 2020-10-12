rm -rf mito_test.svt/

svteaser -h

svteaser surv_sim chrM.fa mito_test --debug
svteaser sim_reads mito_test.svt
