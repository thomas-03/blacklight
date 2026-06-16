import h5py
import numpy as np
import os

def convert_athdf_doubles_to_floats(input_file, output_file=None):
    """
    Converts all float64 datasets in a .athdf (HDF5) file to float32 and saves to a new file.

    Parameters:
        input_file (str): Path to the original .athdf file.
        output_file (str, optional): Path to the output file. If None, appends '_float32' to the input filename.
    """

    if output_file is None:
        base, ext = os.path.splitext(input_file)
        output_file = f"{base}_float32{ext}"

    def copy_attrs(src, dst):
        for key, val in src.attrs.items():
            dst.attrs[key] = val

    def process_group(src_group, dst_group):
        for name, item in src_group.items():
            if isinstance(item, h5py.Dataset):
                data = item[()]  # Read the data
                if data.dtype == np.float64:
                    data = data.astype(np.float32)
                dst_dset = dst_group.create_dataset(name, data=data, compression=item.compression)
                copy_attrs(item, dst_dset)
            elif isinstance(item, h5py.Group):
                new_group = dst_group.create_group(name)
                copy_attrs(item, new_group)
                process_group(item, new_group)

    with h5py.File(input_file, 'r') as src_file:
        with h5py.File(output_file, 'w') as dst_file:
            copy_attrs(src_file, dst_file)
            process_group(src_file, dst_file)

    print(f"Conversion complete: {output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert float64 datasets in a .athdf file to float32.")
    parser.add_argument("input_file", help="Path to the input .athdf file.")
    parser.add_argument("--output_file", default=None,help="Path to the output file. If not provided, appends '_float32' to the input filename.")

    args = parser.parse_args()
    convert_athdf_doubles_to_floats(args.input_file, args.output_file)