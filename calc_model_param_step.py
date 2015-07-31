import project_tools
import project_tools.parameter_fitting.FRET.compute_Jacobian as compute
import model_builder as mdb
import numpy as np

model, fitopts = mdb.inputs.load_model('1PB7',dry_run=True)
def_FRET_pairs = [[219, 371]]
defspacing = 0.1

sim_feature, sim_feature_err, Jacobian, Jacobian_err = compute.calculate_average_Jacobian(model,fitopts, FRET_pairs=def_FRET_pairs, spacing=defspacing)
target_feature, target_feature_err = compute.get_target_feature(model, fitopts)

np.savetxt('target_feature.dat',target_feature)
np.savetxt('sim_feature.dat',sim_feature)
np.savetxt('Jacobian.dat', Jacobian)

import project_tools.parameter_fitting.newton_solver.Truncated_SVD as truncated

truncated.find_solutions(model,scaling=False)
