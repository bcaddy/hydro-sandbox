import h5py
import numpy as np
from numpy.core.defchararray import array


filePath0 = '/ihome/eschneider/rvc9/Code/cholla/bin/tDISKtemp_test/1/1_particles.h5.0'
filePath1 = '/ihome/eschneider/rvc9/Code/cholla/bin/tDISKtemp_test/1/1_particles.h5.1'

file0 = h5py.File(filePath0, 'r')
file1 = h5py.File(filePath1, 'r')

array0 = np.array(file0['particle_IDs'])
array1 = np.array(file1['particle_IDs'])

duplicates = np.intersect1d(array0, array1)

print(duplicates.shape)
# for num in duplicates:
#     print('Duplicated Number: ', num)

trash, counts = np.unique(array0, return_counts=True)

