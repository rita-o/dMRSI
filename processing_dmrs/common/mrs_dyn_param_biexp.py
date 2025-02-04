# Biexp parametrization
Parameters = {
    'Phi_0': 'variable',
    'Phi_1': 'fixed',
    'conc': {'other': {'dynamic': 'model_biexp', 'params': ['c_amp', 'c_adc_slow', 'c_adc_fast', 'c_frac_slow']},
             'Mac': {'dynamic': 'model_lin', 'params': ['c_amp', 'c_slope']}},
    'eps': 'fixed',
    'gamma': 'fixed',
    'baseline': {'dynamic': 'model_exp_offset', 'params': ['b_amp', 'b_adc', 'b_off']}
}
# Optionally define bounds on the parameters
Bounds = {
    'c_amp': (0, None),
    'c_adc_slow': (0, .1),
    'c_adc_fast': (.1, 4),
    'c_frac_slow': (0, 1),
    'gamma': (0, None),
    'b_amp': (None, None),
    'b_adc': (1E-5, 3),
    'b_off': (None, None)
}