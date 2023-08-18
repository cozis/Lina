import os
os.environ['OMP_NUM_THREADS'] = '1'

import numpy as np
import time

N = 1024

if __name__ == "__main__":
    A = np.random.randn(N,N).astype(np.float64)
    B = np.random.randn(N,N).astype(np.float64)

    start = time.monotonic()
    C = A @ B
    stop = time.monotonic()
    
    s = stop-start
    
    ops = 2*N*N*N
    
    print(f"NUMPY: {ops/s * 1e-9} GFLOPS\n")
    