"""
Awais .. Phase Unwrapping
"""

import numpy as np
import matplotlib.pyplot as plt
#from skimage import data, img_as_float, color, exposure
from skimage.restoration import unwrap_phase
from scipy.io import loadmat
import pandas as pd
import skimage.measure as sk
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy.io import savemat

'''
# Simulation

# Load an image as a floating-point grayscale
image = color.rgb2gray(img_as_float(data.chelsea()))
# Scale the image to [0, 4*pi]
image = exposure.rescale_intensity(image, out_range=(0, 4 * np.pi))
# Create a phase-wrapped image in the interval [-pi, pi)
image_wrapped = np.angle(np.exp(1j * image))

'''

data = loadmat('E:\\AWS_Research\\Phase Unwrap\\Data\\BoundaryLeaksData_phase_test.mat')


image_wrapped = data['phase_test']  # variable in mat file 

image_unwrapped = unwrap_phase(image_wrapped)


'''
# Perform phase unwrapping



x = np.array([[3, 2, 1, 2,4,3]])

img = np.vstack([np.zeros_like(x), x, x, x, np.zeros_like(x)])

LP = sk.profile_line(img, (2, 1), (2, 4))

fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)
fig.colorbar(ax.imshow(img, cmap='gist_rainbow'), ax=ax)
plt.show()


plt.plot(LP)
plt.ylabel('some numbers')
plt.show()
'''
# Line Profile

#LP = sk.profile_line(image_unwrapped, (600, 1), (600, 1000))

# 2D Image

fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)
fig.colorbar(ax.imshow(image_wrapped , cmap='jet'), ax=ax)
plt.show()
#plt.savefig('C:\\Users\\Awais\\Downloads\\path.eps', format='eps')





matEXP = {"data": image_unwrapped, "label": "experiment"}
savemat("matlab_python2.mat", matEXP)


'''
plt.plot(LP)
plt.ylabel('some numbers')
plt.show()
'''


