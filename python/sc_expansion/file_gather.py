import h5py
import glob
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('order', type=int, help='Expansion order')
parser.add_argument('U', type=float, help='Interaction strength U')
parser.add_argument('beta', type=float, help='Inverse temperature beta')
parser.add_argument('alpha', type=float, help='Alpha parameter for mcmc')

args = parser.parse_args()

order = args.order
U = args.U
beta = args.beta
alpha = args.alpha

filename = f"./results/full_lattice_Omega_order_{order}_scan_mu_U={U}_beta={beta}_alpha={alpha}.h5"

# Creates the folder path if it's missing
os.makedirs(os.path.dirname(filename), exist_ok=True)

# 1. Find all the tiny files
files = glob.glob("results/full_lattice_data_mu_*.h5")
files.sort() 

if not files:
    print("No subfiles found to merge. Exiting.")
    exit()

# 2. Process and Merge
# Mode "a" opens for reading and writing, creating it if it doesn't exist.
with h5py.File(filename, "a") as master:
    
    # Storage for the current batch of data
    new_mus = []
    new_means = []
    new_errors = []

    print(f"Reading {len(files)} subfiles...")
    for fname in files:
        with h5py.File(fname, "r") as f:
            new_mus.append(f["mu"][()])
            new_means.append(f["mean"][()])
            new_errors.append(f["error"][()])

    # Helper function to handle the HDF5 logic
    def update_dataset(name, data_list):
        new_data = np.array(data_list)
        
        if name in master:
            # APPEND LOGIC
            dset = master[name]
            curr_size = dset.shape[0]
            # Expand the dataset to fit new data
            dset.resize((curr_size + len(new_data),))
            # Write new data starting at the old end-point
            dset[curr_size:] = new_data
            print(f"Appended {len(new_data)} items to {name}.")
        else:
            # CREATE LOGIC
            # maxshape=(None,) makes the dimension resizable
            master.create_dataset(
                name, 
                data=new_data, 
                maxshape=(None,), 
                chunks=True
            )
            print(f"Created new dataset {name} with {len(new_data)} items.")

    # Apply to all three datasets
    update_dataset("mu_list", new_mus)
    update_dataset("mean_list", new_means)
    update_dataset("error_list", new_errors)

    # 3. Cleanup
    print("Cleaning up temporary files...")
    for fname in files:
        try:
            os.remove(fname)
        except OSError as e:
            print(f"Error deleting {fname}: {e}")
            
    print("Cleanup complete.")

print(f"Successfully updated master file: {filename}")