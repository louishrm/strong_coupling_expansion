import h5py
import glob
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('U', type=float, help='Interaction strength U')
parser.add_argument('beta', type=float, help='Inverse temperature beta')
parser.add_argument('alpha', type=float, help='Alpha parameter for mcmc')

args = parser.parse_args()
U = args.U
beta = args.beta
alpha = args.alpha


filename = f"./results/full_lattice_Omega_scan_mu_U={U}_beta={beta}_alpha={alpha}.h5"

# Creates the folder path if it's missing
os.makedirs(os.path.dirname(filename), exist_ok=True)

# 1. Find all the tiny files
files = glob.glob("results/full_lattice_data_mu_*.h5")
files.sort() # Ensure sorted order

# 2. Create the master file
with h5py.File(filename, "w") as master:
    
    # storage arrays
    mus = []
    means = []
    errors = []

    for fname in files:
        with h5py.File(fname, "r") as f:
            # Read single file data
            mu = f["mu"][()]
            mean = f["mean"][()]
            err = f["error"][()]
            
            # Store in lists
            mus.append(mu)
            means.append(mean)
            errors.append(err)
            
            # Optional: Copy group into master if you want detailed structure
            # f.copy("/", master, name=f"mu_{mu:.4f}")

    # 3. Save combined arrays for easy plotting
    master.create_dataset("mu_list", data=np.array(mus))
    master.create_dataset("mean_list", data=np.array(means))
    master.create_dataset("error_list", data=np.array(errors))

    print("Cleaning up temporary files...")
    for fname in files:
        try:
            os.remove(fname)
        except OSError as e:
            print(f"Error deleting {fname}: {e}")
            
    print("Cleanup complete.")
    
print("Successfully merged files into final_results.h5")