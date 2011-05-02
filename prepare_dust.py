import numpy as np
np.seterr(all='ignore')

from hyperion.dust import BHDust

# Read in BHMie result files
d = BHDust('dust/d03_5.5_3.0_A/d03_5.5_3.0_A')

# Extrapolate
d.optical_properties._extrapolate(9.e-4, 1.e7)

# Make a plot
d.plot('dust/d03_5.5_3.0_A.png')

# Use slow dust sublimation
d.set_sublimation_temperature('slow', temperature=1600)

# Write out dust file for this minimum temperature    
d.write('mysg/data/dust/d03_5.5_3.0_A.hdf5')
