import numpy as np

import pyvinecopulib as pv

col_types = ["c", "c", "d", "d"]

# only 10 continuous cols
unifs = np.random.uniform(size=(100, len(col_types)))  # U = F(X)
unifs_with_lag = np.hstack([
    unifs.copy(), np.zeros((100, len(col_types) - 2))
])  #copy of data with the U^- = F(X-1) discrete columns
# "discretise" discrete columns (all bernoulli with p=thres_i)
for i in range(2, len(col_types)):
  thres = np.random.uniform()
  ind = unifs[:, i] <= thres
  unifs[ind, i] = 1
  unifs[~ind, i] = 1 - thres
  unifs_with_lag[:, i] = unifs[:, i]
  unifs_with_lag[ind, len(col_types) + i - 2] = 1 - thres
  unifs_with_lag[~ind, len(col_types) + i - 2] = 0

controls = pv.FitControlsVinecop(family_set=[
    pv.BicopFamily.clayton, pv.BicopFamily.frank, pv.BicopFamily.gaussian,
    pv.BicopFamily.indep, pv.BicopFamily.student, pv.BicopFamily.gumbel,
    pv.BicopFamily.joe, pv.BicopFamily.bb1, pv.BicopFamily.bb6,
    pv.BicopFamily.bb7, pv.BicopFamily.bb8
],
                                 num_threads=8,
                                 parametric_method="mle",
                                 selection_criterion="aic",
                                 show_trace=True)
cop = pv.Vinecop(unifs_with_lag, controls=controls, var_types=col_types)
cop.to_json("./test_save.json")

loaded = pv.Vinecop(filename="./test_save.json")  # Error on this line
