import numpy as np 
import matplotlib as mpl

def complex_array_to_rgb(X, theme='dark', rmax=None):
    '''Takes an array of complex number and converts it to an array of [r, g, b], 
    where phase gives hue and saturaton/value are given by the absolute value.
    Especially for use with imshow for complex plots.'''
    absmax = rmax or np.abs(X).max()
    Y = np.zeros(X.shape + (3,), dtype='float')
    Y[..., 0] = np.angle(X) / (2 * np.pi) % 1 
    if theme == 'light':
        Y[..., 1] = np.clip(np.abs(X) / absmax, 0, 1)
        Y[..., 2] = 1 
    elif theme == 'dark':
        Y[..., 1] = 1 
        Y[..., 2] = np.clip(np.abs(X) / absmax, 0, 1)
    elif theme == 'color':
        Y[..., 1] = 1
        Y[..., 2] = 1
    elif theme == 'weird':
        Y[..., 1] = np.clip(np.abs(X) / absmax, 0, 1)
        Y[..., 2] = np.clip(np.abs(X) / absmax, 0, 1)    
    elif theme == 'wiggly':
        r = np.log2(1. + np.clip(np.abs(X) / absmax, 0, 1))
        Y[..., 1] = (1. + np.abs(np.sin(2. * np.pi * r))) / 2.
        Y[..., 2] = (1. + np.abs(np.cos(2. * np.pi * r))) / 2.

    Y = mpl.colors.hsv_to_rgb(Y)
    return Y
