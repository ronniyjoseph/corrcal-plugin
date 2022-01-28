import numpy as np
from scipy.constants import c

def beam_width(frequency =150e6, diameter=4, epsilon=0.42):
    sigma = epsilon * c / (frequency * diameter)
    width = np.sin(0.5 * sigma)
    return width