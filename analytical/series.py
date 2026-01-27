import numpy as np
from typing import Callable
from numpy.typing import NDArray


def partial_sum(z: NDArray, N: int, f: Callable, max_order: int, r: float, *args, **kwargs) -> NDArray:
    coeffs = calculate_coeffs(N, f, max_order+1, r, *args, **kwargs)
    z = np.asarray(z)
    # Broadcast over the last axis so z can be a scalar or array
    powers = z[..., None] ** np.arange(max_order+1)
    return np.sum(coeffs * powers, axis=-1)


def calculate_coeffs(N: int, f: Callable, max_order: int, r: float, *args, **kwargs) -> NDArray:
    """Calculate power series coefficients for a function f(z).
    N is the number of points on the unit circle to sample.
    f is a function that calculates f(z)
    max_order is the maximum order of the power series to calculate.
    """

    omegas = np.exp(2j*np.pi*np.arange(N)/N)
    fks = f(r*omegas, *args, **kwargs)
    coeffs = np.sum( fks*1/(omegas**np.arange(max_order)[:,None]), axis=1)/(r**np.arange(max_order)*N)

    return coeffs



def estimate_radius(N:int, f:Callable, max_order:int, r:float, *args, **kwargs):

    coeffs = calculate_coeffs(N, f, max_order, r, *args, **kwargs)
    data= []
    for n,an in enumerate(coeffs):
        
        for m,am in enumerate(coeffs):

            if m!=n:
                if np.abs(an) < 1e-10 or np.abs(am) < 1e-10:
                    continue
                
                data.append(np.abs((an/am)**(1/(m-n))))

    return np.median(data)