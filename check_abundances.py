""" 
Could Kronos' abundances be due to incorrect surface gravity?
"""

import numpy as np
from astropy.table import Table 

from solar_abundances import asplund_2009 as solar_abundances

abundances = Table.read("brewer_2016_table9.txt", format="cds")

kronos = np.where(abundances["Name"] == "HD 240430")[0][0]

# Get the abundances.
kronos_abundances = {}
for col in abundances.dtype.names[2:]:
    element = col[1:].split("/")[0]
    
    x_h = abundances[col][kronos]
    log_eps = x_h + solar_abundances[element]
    kronos_abundances[element] = log_eps

# Calculate the corresponding equivalent widths using the
# abundances and the line list in Brewer et al. 2016
line_list = Table.read("brewer_2016_table7.txt", format="cds")

use = np.zeros(len(line_list), dtype=bool)
expected_log_eps = np.nan * np.ones(len(line_list))
for element, log_eps in kronos_abundances.items():
    match = (line_list["Element"] == element) \
          * (line_list["Masked"] != "x")

    use[match] = True
    expected_log_eps[match] = log_eps

# Discard other lines.
line_list = line_list[use]
line_list["expected_log_eps"] = expected_log_eps[use]

# Generate model photosphere for the reported stellar parameters.
# Assume v_micro = 1.1 km/s
MARCS = smh.photospheres.MARCSInterpolator()
expected_photosphere = MARCS(5803, 4.33, 0.20, 1.1)

expected_ew = smh.rt.moog.ewfind(line_list, expected_photosphere)


# Using those equivalent widths, calculate what the abundances
# would be if the logg = 4.43 instead of 4.33
alternative_photosphere = MARCS(5803, 4.43, 0.20, 1.1)
alternative_abundances = smh.rt.moog.abundances(expected_ew, alternative_photosphere)

# Summarise:


