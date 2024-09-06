import cupy as cp
import numpy as np
import datetime

a = np.random.random(10000).astype(np.complex64)
cp.allclose(a, a)
