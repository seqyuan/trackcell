"""Example usage for square-bin visualization in TrackCell."""

import trackcell as tcl


def main():
    adata_bin = tcl.io.read_hd_bin(
        datapath="SpaceRanger4.0/Cse1/binned_outputs/square_016um",
        sample="Cse1",
        binsize=16,
    )

    # 1. Show H&E image only with coordinate range
    tcl.pl.spatial_squarebin(adata_bin, color=None)

    # 2. Plot a gene on square bins
    tcl.pl.spatial_squarebin(
        adata_bin,
        color="EPCAM",
        cmap="Reds",
        alpha=0.8,
        alpha_img=0.4,
        rasterized=True,
    )

    # 3. Crop a region of interest
    tcl.pl.spatial_bin(
        adata_bin,
        color="EPCAM",
        crop_coord=(16000, 18000, 14000, 18000),
        cmap="Reds",
    )


if __name__ == "__main__":
    main()
