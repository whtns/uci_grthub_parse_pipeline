import loompy
import numpy as np
import argparse

def subsample_loom(input_loom_path, output_loom_path, percentage=0.1):
    """
    Subsamples a loom file to a given percentage of its cells.
    
    Args:
        input_loom_path (str): The path to the input .loom file.
        output_loom_path (str): The path where the subsampled .loom file will be saved.
        percentage (float): The percentage of cells to keep (e.g., 0.1 for 10%).
    """
    try:
        # Open the input loom file in read-only mode
        with loompy.connect(input_loom_path, "r") as ds:
            # Get the total number of cells (columns)
            n_cells = ds.shape[1]
            print(f"Original file has {n_cells} cells.")

            # Determine the number of cells to select
            n_select = int(n_cells * percentage)
            if n_select == 0 and n_cells > 0:
                n_select = 1  # Ensure at least one cell is selected if the file isn't empty
            print(f"Subsampling to {n_select} cells ({percentage * 100:.0f}%).")

            # Get a random subset of column indices
            # np.random.choice is useful for sampling without replacement
            selected_col_indices = np.sort(np.random.choice(n_cells, size=n_select, replace=False))

            # Create a new loom file with the subsampled data
            # The columns (cells) are sliced using the selected indices
            loompy.create(output_loom_path, ds[:, selected_col_indices], row_attrs=ds.ra, col_attrs=ds.ca[selected_col_indices])

        print(f"Successfully subsampled and saved to {output_loom_path}")

    except FileNotFoundError:
        print(f"Error: The file '{input_loom_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Subsample a loom file to a given percentage of cells."
    )
    
    # Add arguments
    parser.add_argument(
        "input",
        type=str,
        help="Path to the input .loom file."
    )
    parser.add_argument(
        "output",
        type=str,
        help="Path to the output .loom file."
    )
    parser.add_argument(
        "--percentage",
        type=float,
        default=0.1,
        help="The percentage of cells to keep (e.g., 0.1 for 10%%). Defaults to 0.1."
    )
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the subsampling function with the parsed arguments
    subsample_loom(args.input, args.output, args.percentage)