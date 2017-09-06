
data(Ec_core, package="sybil")
ecoli <- as_GRBmodel(Ec_core)
ecoli$show_output(FALSE)

fva <- flux_variability(ecoli, obj_frac=NA)

fc1 <- flux_coupling(ecoli, min_fva_cor = 0.0)
fc2 <- flux_coupling(ecoli, min_fva_cor = 0.9)
