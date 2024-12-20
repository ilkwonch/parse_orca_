import os
import h5py
import argparse
import numpy as np
from parse_orca_output import parse_orca_output, process_parsed_data


def generate_orca_input(charge, save_filename):
    """
    Generate an ORCA input file.
    """
    input_file = f"{save_filename.replace('.xyz', '.inp')}"
    with open(input_file, "w") as f:
        f.write("!b97-3c tightscf engrad SCFConvForced slowconv notrah\n")
        f.write("%SCF\n")
        f.write("    maxiter 500\n")
        f.write("END\n")
        f.write("%elprop dipole true quadrupole true end\n")
        f.write("%output Print[P_DFTD_GRAD] 1 Print[P_Hirshfeld] 1 end\n")
        f.write(f"* XYZFILE {charge} 1 {save_filename}\n")
    return input_file

def save_xyz_file(coord, number, save_filename):
    """
    Save atomic coordinates and numbers into an XYZ file.
    """
    atomic_number_to_symbol = {1: 'H', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
                               14: 'Si', 15: 'P', 16: 'S', 17: 'Cl',
                               33: 'As', 34: 'Se', 35: 'Br', 53: 'I'}

    num_atoms = len(coord)
    base_name = os.path.splitext(os.path.basename(save_filename))[0]

    xyz_content = f"{num_atoms}\n"
    xyz_content += f"{base_name}\n"

    for num, atom_coord in zip(number, coord):
        element_symbol = atomic_number_to_symbol[num]
        xyz_content += f"{element_symbol} {atom_coord[0]:.12f} {atom_coord[1]:.12f} {atom_coord[2]:.12f}\n"

    with open(save_filename, "w") as xyz_file:
        xyz_file.write(xyz_content)

def orca_to_hdf5(hdf5_file, group_name, start, end, orca_path, running_dir):
    """
    Ensure the groups `energy`, `forces`, `dipole`, and `charges` exist in the HDF5 file.
    Parse ORCA outputs and populate these groups.
    """
    with h5py.File(hdf5_file, "a") as f:
        if group_name not in f:
            raise ValueError(f"Group {group_name} not found in HDF5 file.")

        group = f[group_name]
        coords = group["coord"][:]  # (N_files, atomic_num, 3)
        numbers = group["numbers"][:]  # (atomic_num,)
        charges = group["charge"][:]  # (N_files,)
        

        if "energy" not in group:
            group.create_dataset("energy", (len(coords),), dtype="float64")
        if "forces" not in group:
            group.create_dataset("forces", (len(coords), numbers.shape[0], 3), dtype="float64")
        if "dipole" not in group:
            group.create_dataset("dipole", (len(coords), 3), dtype="float64")
        if "charges" not in group:
            group.create_dataset("charges", (len(coords), numbers.shape[0]), dtype="float64")
        
        end = end if end is not None else len(coords)
        #for i in range(start, min(end, len(coords))):
        for i in range(start, end):

            if group["energy"][i] != 0.0:
                continue
            
            
            xyz_filename = os.path.join(running_dir, f"orca_{i}.xyz")
            input_file = generate_orca_input(charges[i], xyz_filename)
            output_file = os.path.join(running_dir, f"orca_{i}.out")

            save_xyz_file(coords[i], numbers, xyz_filename)

            success = True
            try:
                os.system(f"{orca_path} {input_file} > {output_file}")

                if not os.path.exists(output_file):
                    raise FileNotFoundError(f"ORCA output file {output_file} was not generated.")

                parsed_data = parse_orca_output(output_file)
                # If the parsed_data indicates no energy, assume no convergence
                if parsed_data["total_energy"] is None:
                    success = False

                processed_data = process_parsed_data(parsed_data)

                if processed_data["energy"] is not None:
                    group["energy"][i] = processed_data["energy"]
                if processed_data["forces"] is not None:
                    group["forces"][i] = processed_data["forces"]
                if processed_data["dipole"] is not None:
                    group["dipole"][i] = processed_data["dipole"]
                if processed_data["charge"] is not None:
                    group["charges"][i] = processed_data["charge"]

            except Exception as e:
                print(f"Error processing conformation {i}: {e}")
                # Mark as unsuccessful if there's an exception
                success = False

            finally:
                if success:
                    # Extract the base name of the file (e.g., orca_6)
                    base_name = os.path.splitext(os.path.basename(xyz_filename))[0]
                    
                    for file_to_remove in os.listdir(running_dir):
                        full_path = os.path.join(running_dir, file_to_remove)
                        
                        # Match files starting with base_name
                        if file_to_remove.startswith(base_name):
                            if os.path.exists(full_path):
                                os.remove(full_path)
                else:
                    for file_to_remove in os.listdir(running_dir):
                        full_path = os.path.join(running_dir, file_to_remove)
                        if not file_to_remove.endswith(('.inp', '.xyz', '.out')):
                            if os.path.exists(full_path):
                                os.remove(full_path)            


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ORCA calculations and populate HDF5 file.")
    parser.add_argument("--hdf5_file", required=True, help="Path to the HDF5 file.")
    parser.add_argument("--group_name", required=True, help="Group name in the HDF5 file.")
    parser.add_argument("--start", type=int, required=True, help="Start index for conformations.")
    #parser.add_argument("--end", type=int, required=True, help="End index for conformations.")
    parser.add_argument("--end", type=int, help="End index for conformations (optional, defaults to all).")
    parser.add_argument("--orca_path", required=True, help="Path to the ORCA executable.")
    parser.add_argument("--running_dir", required=True, help="Directory for input files.")

    args = parser.parse_args()

    os.makedirs(args.running_dir, exist_ok=True)

    print(f"Starting ORCA calculations for group: {args.group_name} from {args.start} to {args.end}...")
    orca_to_hdf5(args.hdf5_file, args.group_name, args.start, args.end, args.orca_path, args.running_dir)
    print("ORCA calculations completed.")

