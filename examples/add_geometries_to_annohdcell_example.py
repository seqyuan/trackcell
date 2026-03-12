#!/usr/bin/env python
"""
Example script for adding geometries to annohdcell's final cell h5ad output.

This script demonstrates how to add polygon geometries to annohdcell's final
cell-level h5ad file (b2c_cell.h5ad) using the bin-level h5ad (b2c_2um.h5ad)
that contains cell assignment labels.
"""

import trackcell.io as tcio
import trackcell.pl as tcpl

# Example usage
if __name__ == "__main__":
    # Path to annohdcell's 2μm bin h5ad file (with cell labels)
    bin_h5ad_path = "b2c_2um.h5ad"

    # Path to annohdcell's final cell h5ad file (without geometries)
    cell_h5ad_path = "b2c_cell.h5ad"

    # Output path for cell h5ad with geometries
    output_h5ad_path = "b2c_cell_with_geom.h5ad"

    # Add geometries to the final cell h5ad
    adata = tcio.add_geometries_to_annohdcell_output(
        bin_h5ad_path=bin_h5ad_path,
        cell_h5ad_path=cell_h5ad_path,
        output_h5ad_path=output_h5ad_path,
        sample="sample1",
        labels_key="labels_joint",  # Column in bin h5ad with cell labels
        bin_size_um=2.0,            # Bin size in micrometers
        buffer_polygons=True        # Buffer polygons to account for bin size
    )

    print(f"\nGeometry addition complete!")
    print(f"Output: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"Saved to: {output_h5ad_path}")

    # Now you can use trackcell visualization functions
    print("\nYou can now visualize with trackcell:")
    print("  import trackcell.pl as tcpl")
    print("  tcpl.spatial_cell(adata, sample='sample1')")

    # Optional: visualize immediately
    # tcpl.spatial_cell(adata, sample="sample1", show=True)
