import numpy as np
import pandas as pd
import date_routines as dr
import smoothing_routines as sm


import matplotlib.pyplot as plt
import gpflow as gp

#TODO: WEG
from importlib import reload


reload(dr)


def smooth_phi_estimates(phiPerBin_table, smoothingBandwidth=7):

	date_zero = dr.get_date_zero(phiPerBin_table['t'].iloc[0], phiPerBin_table['date'].iloc[0])

	smoothedPhi_table = sm.smooth_values(x=phiPerBin_table.t, y=phiPerBin_table.phi_d, bandwidth=smoothingBandwidth)

	# TODO WEG
	# X_ = np.array(
	# 	[
	# 		[0.865], [0.666], [0.804], [0.771], [0.147], [0.866], [0.007], [0.026],
	# 		[0.171], [0.889], [0.243], [0.028],
	# 	]
	# )
	# Y_ = np.array(
	# 	[
	# 		[1.57], [3.48], [3.12], [3.91], [3.07], [1.35], [3.80], [3.82], [3.49],
	# 		[1.30], [4.00], [3.82],
	# 	]
	# )

	# X = np.transpose(np.array([phiPerBin_table['t']]))
	# Y = np.transpose(np.array([phiPerBin_table['phi_d']]))

	# k = gp.kernels.Matern32(variance=1, lengthscales=1.2)
	# model = gp.models.GPR(
	# #	(phiPerBin_table['t'], phiPerBin_table['phi_d']),
	# 	(X,Y),
	# 	kernel = k
	# )

	# opt = gp.optimizers.Scipy()
	# opt.minimize(model.training_loss, model.trainable_variables)

	# Xplot = np.arange(min(X), max(X)+1)[:, None]
	# Xplot_ = np.linspace(-0.1, 1.1, 100)[:, None]

	# f_mean, f_var = model.predict_f(Xplot, full_cov=False)
	# y_mean, y_var = model.predict_y(Xplot)

	# f_lower = f_mean - 1.96 * np.sqrt(f_var)
	# f_upper = f_mean + 1.96 * np.sqrt(f_var)
	# y_lower = y_mean - 1.96 * np.sqrt(y_var)
	# y_upper = y_mean + 1.96 * np.sqrt(y_var)

	# plt.clf()
	# plt.plot(X, Y, "kx", mew=2, label="input data")
	# plt.plot(Xplot, f_mean, "-", color="C0", label="mean")
	# plt.plot(Xplot, f_lower, "--", color="C0", label="f 95% confidence")
	# plt.plot(Xplot, f_upper, "--", color="C0")
	# plt.fill_between(
	# 	Xplot[:, 0], f_lower[:, 0], f_upper[:, 0], color="C0", alpha=0.1
	# )
	# plt.plot(Xplot, y_lower, ".", color="C0", label="Y 95% confidence")
	# plt.plot(Xplot, y_upper, ".", color="C0")
	# plt.fill_between(
	# 	Xplot[:, 0], y_lower[:, 0], y_upper[:, 0], color="C0", alpha=0.1
	# )
	# plt.legend()
	# plt.show()
	
	return pd.DataFrame({'t': smoothedPhi_table.x.astype(int),
		'date': dr.add_days_to_date(smoothedPhi_table.x.to_list(), date_zero),
		'phiSmoothed': smoothedPhi_table.f_mean,
		'phi_lower': smoothedPhi_table.f_lower,
		'phi_upper': smoothedPhi_table.f_upper}) 