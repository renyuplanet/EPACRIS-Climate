## python script to plot and compare TP profiles from EPACRIS output

import os
import numpy as np
import matplotlib.pyplot as plt

# Get all the files in the Results directory
results_dir = 'Results/K2-18b/'
cases = ['cloud_test_off_rad_new_test_10', 'cloud_test_on_rad_new_test_10']
tp_file_extension = 'NewTemperature.dat'

# Plot the TP profiles
for case in cases:
    tp_file = os.path.join(results_dir, case, tp_file_extension)
    data = np.loadtxt(tp_file)
    # Convert pressure from log10(Pascal) to bar
    pressure = 10**(data[:, 1] - 5)
    plt.plot(data[:, 2], pressure, label=case,lw=2)
# Ticks inwards
plt.tick_params(axis='both', direction='in')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (bar)')
plt.title('Temperature-Pressure Profile diff')
plt.yscale('log')
plt.ylim(1e-8, 1e3)
plt.gca().invert_yaxis()
plt.grid(False)
plt.savefig('TP_profile_new5_test_diff.png')
#plt.show()