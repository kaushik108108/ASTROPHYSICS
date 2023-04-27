import phoebe
from phoebe import u
import numpy as np
from scipy import interpolate as interp
import matplotlib.pyplot as plt

logger = phoebe.logger()

r_band_width = phoebe.c.c * (1/(550 * u.nm) - 1/(700 * u.nm))
period = 0.2207213 * u.day
T_1 = 4750 * u.K
T_2_by_T_1 = 0.965935

observations = {
  'magnitudes' : np.array([14.7114, 15.0301, 14.822 , 14.6218, 14.8402, 
                           14.9702, 14.6378, 14.5788, 15.1366, 14.8615, 
                           14.5711, 14.4997, 14.5798, 14.5018])
  ,
  'times' : np.array([2460018.40143111, 2460018.41957866, 2460018.44037724, 2460020.44282839, 
                      2460026.35522657, 2460026.37828184, 2460026.39928253, 2460029.3165532,
                      2460031.33677605, 2460039.29658565, 2460039.31466435, 2460039.33430556, 
                      2460039.3552662, 2460042.31506944]) * u.day
  ,
  'mag_errors' : np.array([8.669e-3, 9.0969e-3, 6.929e-3, 1.122e-2, 
                           5.905e-3, 7.901e-3, 4.879e-3, 9.091e-3, 
                           1.015e-2, 1.1703e-2, 8.672e-3, 8.921e-3, 
                           9.516e-3, 8.829e-3])
}

observations['wrapped_times'] = observations['times'] % period
observations['phases'] = observations['wrapped_times'] / period
observations['fluxes'] = (3631 * u.Jy * 10**( - observations['magnitudes'] / 2.5 ) * r_band_width).to('W/m^2')
observations['flux_errors'] = observations['fluxes'] * np.log(10) / 2.5 * observations['mag_errors']

## Plot the Cubic Interpolant
# interpolated_light_curve = interp.CubicSpline(np.sort(observations['phases']), 
#                                               [observations['fluxes'][i].value for i in np.argsort(observations['phases'])],
#                                               bc_type='periodic')
# fake_phase = np.arange(0,1,0.01)
# fake_flux = interpolated_light_curve(fake_phase) * u.Unit('W/m^2')
# plt.plot(x= fake_phase, y=fake_flux)
# plt.scatter(x=observations['phases'], y=observations['fluxes'])
# plt.show()

b = phoebe.default_contact_binary()

b.flip_constraint('pot', solve_for='requiv@primary')
## mlp solver appeasement
# b.flip_constraint('requivsumfrac', solve_for='sma')
b.flip_constraint('ecosw', solve_for='per0')
b.flip_constraint('esinw', solve_for='ecc')

# Gravity darkening coeffcients
b.set_value('gravb_bol@primary'  , value = 0.32 )
b.set_value('gravb_bol@secondary', value = 0.32 )
# Known Temperature and ratio
b.set_value('teff@primary', value = T_1)
b.flip_constraint('teffratio@constraint', solve_for='teff@secondary')
b.set_value('teffratio@component', value = T_2_by_T_1)
# Known period
b.set_value('period@binary', value = period)


## from tutorial
b.add_dataset('lc', times=observations['times'],
                    fluxes=observations['fluxes'], 
                    sigmas=observations['flux_errors'], 
                    compute_phases=phoebe.linspace(0,1,101))
                    # compute_phases=observations['phases'])
                    # compute_times=observations['times'] % period)
b.set_value('passband@lc01', value='SDSS:rprime')
b.set_value_all('ld_mode', 'manual')
b.set_value_all('ld_mode_bol', 'manual')
b.set_value_all('atm', 'blackbody')
b.set_value('pblum_mode', 'dataset-scaled')

b.add_solver('estimator.ebai', ebai_method='knn', solver='ebai_knn')
b.run_solver('ebai_knn', solution='ebai_knn_solution')
# Remove q from solving domain
b.set_value('adopt_parameters', value=b['adopt_parameters'].value[:-1])

## Vary mass ratios and find minima
# q_0 = 0.6
# q_n = 0.8
# dq = 0.01
# q = q_0
# data = []
# sorted_fluxes = observations['fluxes'     ][ np.argsort(observations['wrapped_times']) ]
# sorted_errors = observations['flux_errors'][ np.argsort(observations['wrapped_times']) ]
# while (q <= q_n) :
#   b.set_value('q@binary@component', value = q)
#   b.adopt_solution('ebai_knn_solution')
#   b.run_compute(model='ebai_knn_model', overwrite=True)
#   data.append( (q, np.sum( ((sorted_fluxes.value - b['fluxes@model'].value) / sorted_errors.value)**2 ) ) )
#   q += dq
# plt.xlabel('Mass ratio q')
# plt.ylabel('chi^2 error')
# plt.plot(np.array(data)[:,0], np.array(data)[:,1], 'o--')
# plt.show()

# Value of minima
b.set_value('q', value = 0.7)
print(b.adopt_solution('ebai_knn_solution'))

b.run_compute(model='ebai_knn_model')
_ = b.plot(x='phase', ls='-', legend=True, show=True)

b.add_dataset('mesh', compute_phases=[0.25], dataset='mesh01')
b.run_compute(model='mesh01')
afig, mplfig = b['mesh01@model'].plot(x='ws', show=True)

## Converge to paper with this one neat trick
b.set_value('q', value = 1.95)
b.set_value('sma@binary@component', value = 1.45 * u.solRad)
b.set_value('pot@contact_envelope@envelope@component', value = 5.1)

# b.run_compute(model='default')
b.run_compute(model='paper')
_ = b.plot(x='phase', ls='-', legend=True, show=True)
