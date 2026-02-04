"""
Robust Pade library
"""

import numpy as np
from scipy.linalg import toeplitz


def pade_approx(f, pade_orders, tol=1e-14):
    """Pade approximation to a function or Taylor series.
    a, b, mu, nu = pade_approx(f, m, n, tol)
    constructs a Pade approximant to F using the robust algorithm from
    [1] based on the SVD.
    F must be a vector of coefficients f_0, ..., f_{m + n}, and if F
    is a function handle, the function must be analytic in a
    neighborhood of the unit disc, since the coefficients are computed
    via FFT.  M and N are the desired numerator and denominator
    degrees, respectively, and must be nonnegative. The optional TOL
    argument specifies the relative tolerance; if omitted, it defaults
    to 1e-14. Set TOL to 0 to turn off robustness. The output is a
    function handle R of for an exact type (MU, NU) Pade approximant
    to F with coefficient vectors A and B and, optionally, the POLES
    and RESIDUES.
    This code is included in the Chebfun distribution for the
    convenience of readers of _Approximation Theory and Approximation
    Practice_, but it is not actually a Chebfun code. A Chebfun
    analogue is CHEBPADE.
    References:
    [1] P. Gonnet, S. Guettel, and L. N. Trefethen, "ROBUST PADE APPROXIMATION 
        VIA SVD", SIAM Rev., 55:101-117, 2013.
    See also CHEBPADE.
    Copyright 2014 by The University of Oxford and The Chebfun
    Developers.  See http://www.chebfun.org/ for Chebfun information.
    Translated to python by Jim Fowler.
    Added functionality by Fedor Simkovic
    """
    eps = np.finfo(float).eps
    
    if type(pade_orders[0])!=int:
        pade_orders=pade_orders[0]
    m,n = pade_orders

    def discard_trailing_epsilon(v):
        if np.abs(v)[-1] <= tol:
            return v[ : np.where( np.abs(v) > tol )[0][-1] + 1 ]
        else:
            return v
    
    #  Make sure c is long enough but not longer than necessary.
    c = f[:m + n + 1]
    if len(c) < m + n + 1:
        c = np.pad(f, (0, m + n + 1 - len(c)), 'constant')
    
    #  Compute absolute tolerance.
    ts = tol*np.linalg.norm(c)
    
    ## Compute the Pade approximation.

    # check for the special case r = 0
    if ( np.linalg.norm(c[0:m+1], np.inf) <= tol*np.linalg.norm(c, np.inf) ):
        a = 0
        b = 1
        mu = -np.inf
        nu = 0
    else:
        ## the general case
    
        #  First row/column of Toeplitz matrix.
        row = np.zeros(n+1, complex)
        row[0] = c[0]
        col = c
    
        #  Do diagonal hopping across block.
        while True:
            #  Special case n == 0.
            if n == 0:
                a = c[0:m+1]
                b = 1
                break
    
            #  Form Toeplitz matrix.
            Z = toeplitz(col[0:m+n+1], row[0:n+1])
    
            #  Compute numerical rank.
            C = Z[m+1:m+n+1,:]
            rho = np.sum(np.linalg.svd(C, compute_uv=False ) > ts)
            if rho == n:
                break
    
            #  Decrease mn, n if rank-deficient.
            m = m - (n - rho)
            n = rho
    
        #  Hopping finished. Now compute b and a.
        if n > 0:
            U, S, V = np.linalg.svd(C)

            #  Null vector gives b.
            b = V[:,n]
    
            #  Do final computation via reweighted QR for better zero preservation.
            D = np.diag(np.abs(b) + np.sqrt(eps))
            Q, R = np.linalg.qr( np.transpose( np.matmul(C,D) ), mode='complete' )
    
            #  Compensate for reweighting.
            b = np.matmul(D,Q)[:,n]
            b = b/np.linalg.norm(b)
    
            #  Multiplying gives a.
            a = np.dot( Z[0:m+1,0:n+1], b )
    
            #  Count leading zeros of b.
            lam = np.argmax( np.abs(b) > tol )
    
            #  Discard leading zeros of b and a.
            b = b[lam:]
            a = a[lam:]
    
            #  Discard trailing zeros of b.
            b = discard_trailing_epsilon(b)
    
        #  Discard trailing zero coefficients in a.
        a = discard_trailing_epsilon(a)
    
        #  Normalize.
        a = a/b[0]
        b = b/b[0]
    
        #  Exact numerator, denominator degrees.
        mu = len(a)
        nu = len(b)
        
        def r(U):
            val_a = 0
            val_b = 0
            for a_ in range(len(a)):
                val_a += a[a_]*U**a_
            for b_ in range(len(b)):
                val_b += b[b_]*U**b_    
            return val_a/val_b    
        
        # Compute zeros if requested.
        zeros = np.roots(a[::-1]);        
        
        # Compute poles if requested.
        poles = np.roots(b[::-1])

        # Estimate residues.
        t = max(tol, 1e-7);  # Perturbation for residue estimate.
        residues = t*(r(poles + t) - r(poles - t))/2

    return r, a, b, mu, nu, zeros, poles, residues
