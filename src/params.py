'''Parameter file for perfusion and gas exchange model.
'''

__author__ = 'pzuritas'
__email__ = 'pzurita@uc.cl'

from dolfin import Constant


avg_speed_size = 1000  # (Normal is 1000)
mu = 3.5*0.00750062E-3  # mmHg*s
kappa = mu*avg_speed_size/4*300  # um^2s
p_max = 12  # mmHg
p_min = 8  # mmHg

d_pla_O2 = 1.62E-5*1E8  # um^2/s
d_ba_O2 = 1E-5*1E8  # um^2/s
beta_O2 = 0.9*1.41E-6/(1E15) # mol/um^3/mmHg
k_O2 = 20  # 1/s
k_prime_O2 = 66E-6*1E15  # um^3/mol/s

d_pla_CO2 = 1E-5*1E8  # um^2/s
d_ba_CO2 = 0.914E-5*1E8  # um^2/s
beta_CO2 = 0.9*28.2E-6/(1E15)  # mol/um^3/mmHg
k_CO2 = 0.008  # 1/s
k_prime_CO2 = 6E-6*1E15  # um^3/mol/s

c_t = 2.5E-3/(1E15)  # mol/um^3
h_ba = 0.3  # um

p_air_O2 = 100  # mmHg
p_air_CO2 = 0  # mmHg
p_O2_in = 40  # mmHg
p_CO2_in = 45  # mmHg

params = {
    "kappa": kappa,
    "mu": mu,
    "p_max": p_max,
    "p_min": p_min,
    "d_pla_O2": d_pla_O2,
    "d_ba_O2": d_ba_O2,
    "d_pla_CO2": d_pla_CO2,
    "d_ba_CO2": d_ba_CO2,
    "beta_O2": beta_O2,
    "k_O2": k_O2,
    "k_prime_O2": k_prime_O2,
    "beta_CO2": beta_CO2,
    "k_CO2": k_CO2,
    "k_prime_CO2": k_prime_CO2,
    "h_ba": h_ba,
    "p_air_O2": p_air_O2,
    "p_air_CO2": p_air_CO2,
    "p_O2_in": p_O2_in,
    "p_CO2_in": p_CO2_in,
    "c_t": c_t
}