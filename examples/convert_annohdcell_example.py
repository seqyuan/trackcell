#!/usr/bin/env python
"""
Example script for converting annohdcell bin2cell output to trackcell format.

This script demonstrates how to convert a 2μm bin-level h5ad file from annohdcell
(with cell assignment labels) into a trackcell-compatible cell-level h5ad file
with polygon geometries for spatial visualization.
"""

import trackcell.io as tcio
import trackcell.pl as tcpl

# Example usage
if __name__ == "__main__":
    # Path to your annohdcell 2μm bin h5ad file (e.g., b2c_2um.h5ad)
    bin_h5ad_path = "sample_2um.h5ad"

    # Output path for trackcell-compatible h5ad
    output_h5ad_path = "sample_trackcell.h5ad"

    # Convert annohdcell format to trackcell format
    adata = tcio.convert_annohdcell_to_trackcell(
        bin_h5ad_path=bin_h5ad_path,
        output_h5ad_path=output_h5ad_path,
        sample="sample1",
        labels_key="labels_joint",  # Column in .obs with cell labels
        bin_size_um=2.0,            # Bin size in micrometers
        create_polygons=True,       # Create polygon geometries
        buffer_polygons=True        # Buffer polygons to account for bin size
    )

    print(f"\nConversion complete!")
    print(f"Output: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"Saved to: {output_h5ad_path}")

    # Now you can use trackcell visualization functions
    print("\nYou can now visualize with trackcell:")
    print("  import trackcell.pl as tcpl")
    print("  tcpl.spatial_cell(adata, sample='sample1')")

    # Optional: visualize immediately
    # tcpl.spatial_cell(adata, sample="sample1", show=True)
