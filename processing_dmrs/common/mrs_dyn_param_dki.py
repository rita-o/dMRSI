# DKI parametrization
Parameters = {
    'Phi_0': 'variable',
    'Phi_1': 'fixed',
    'conc': {'other': {'dynamic': 'model_dki', 'params': ['c_amp', 'c_adc', 'c_kurtosis', ]},
             'Mac': {'dynamic': 'model_lin', 'params': ['c_amp', 'c_slope']}},
    'eps': 'fixed',
    'gamma': 'fixed',
    'baseline': {'dynamic': 'model_exp_offset', 'params': ['b_amp', 'b_adc', 'b_off']}
}
# Optionally define bounds on the parameters
Bounds = {
    'c_amp': (0, None),
    'c_adc': (0.01, 3),
    'c_kurtosis': (-4, 4),
    'gamma': (0, None),
    'b_amp': (None, None),
    'b_adc': (1E-5, 3),
    'b_off': (None, None)
}
