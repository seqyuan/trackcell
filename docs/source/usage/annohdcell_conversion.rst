Converting annohdcell Output to TrackCell Format
================================================

This guide explains how to convert annohdcell output files into trackcell-compatible
format with polygon geometries for spatial visualization.

Two Conversion Methods
----------------------

**Method 1**: Use :func:`~trackcell.io.convert_annohdcell_to_trackcell` to create a
new cell-level h5ad from scratch (simple count summation).

**Method 2**: Use :func:`~trackcell.io.add_geometries_to_annohdcell_output` to add
geometries to annohdcell's final cell h5ad (preserves exact count aggregation).


Method 1: Convert from 2μm Bin H5AD Only
----------------------------------------

Background
~~~~~~~~~~

**annohdcell Format**

* **Input**: 2μm bin-level h5ad (e.g. ``b2c_2um.h5ad``)
* **Structure**: Each observation is a 2μm bin with cell assignment labels in
  ``.obs["labels_joint"]``
* **Missing**: No polygon geometries for cells, preventing spatial visualization

**trackcell Format**

* **Required**: Cell-level h5ad with polygon geometries
* **Structure**:

  * ``.obs["geometry"]``: WKT strings representing cell boundaries
  * ``.uns["spatial"][sample]["geometries"]``: GeoDataFrame with Shapely Polygon objects
  * ``.obsm["spatial"]``: Cell centroid coordinates

Basic Example
~~~~~~~~~~~~~

.. code-block:: python

   import trackcell as tcl

   adata = tcl.io.convert_annohdcell_to_trackcell(
       bin_h5ad_path="sample_2um.h5ad",
       output_h5ad_path="sample_trackcell.h5ad",
       sample="sample1"
   )
   tcl.pl.spatial_cell(adata, library_id="sample1")

Advanced Options
~~~~~~~~~~~~~~~~

.. code-block:: python

   adata = tcl.io.convert_annohdcell_to_trackcell(
       bin_h5ad_path="b2c_2um.h5ad",
       output_h5ad_path="trackcell_format.h5ad",
       sample="my_sample",
       labels_key="labels_joint",      # Column with cell labels (default)
       bin_size_um=2.0,                # Bin size in micrometers
       create_polygons=True,           # Create polygon geometries
       buffer_polygons=True,           # Buffer polygons by bin size
   )


Method 2: Add Geometries to Existing Cell H5AD
----------------------------------------------

When to Use
~~~~~~~~~~~

Use this method when you want to:

* **Preserve annohdcell's exact count aggregation** (which may differ from simple summation)
* **Keep all metadata** from annohdcell's final output
* **Add only geometries** without recalculating counts

Basic Example
~~~~~~~~~~~~~

.. code-block:: python

   import trackcell as tcl

   adata = tcl.io.add_geometries_to_annohdcell_output(
       bin_h5ad_path="b2c_2um.h5ad",      # 2μm bin h5ad with cell labels
       cell_h5ad_path="b2c_cell.h5ad",    # Final cell h5ad from annohdcell
       output_h5ad_path="b2c_cell_with_geom.h5ad",
       sample="sample1"
   )
   tcl.pl.spatial_cell(adata, library_id="sample1")

What This Method Does
~~~~~~~~~~~~~~~~~~~~~

1. Reads both h5ad files (2μm bin + final cell)
2. Matches cells using ``object_id`` in the cell h5ad
3. Creates convex hull polygons from bin coordinates for each cell
4. Adds polygon geometries to the cell h5ad without modifying counts
5. Preserves all original data from annohdcell's output


Comparison
----------

.. list-table::
   :header-rows: 1

   * - Feature
     - Method 1 (``convert_annohdcell_to_trackcell``)
     - Method 2 (``add_geometries_to_annohdcell_output``)
   * - Input
     - 2μm bin h5ad only
     - 2μm bin h5ad + cell h5ad
   * - Count aggregation
     - Simple sum across bins
     - Uses annohdcell's aggregation
   * - Metadata
     - Basic metadata only
     - Preserves all annohdcell metadata
   * - Speed
     - Faster (one file)
     - Slightly slower (two files)


Common Parameters
-----------------

* **bin_h5ad_path**: 2μm bin h5ad from annohdcell
* **output_h5ad_path**: Path for output h5ad (optional)
* **sample**: Sample name (auto-detected if ``None``)
* **labels_key**: ``.obs`` column with cell labels (default ``"labels_joint"``)
* **bin_size_um**: Bin size in μm (default 2.0)
* **buffer_polygons**: Buffer polygons by bin size (default ``True``)

Method-Specific Parameters
--------------------------

**Method 1** (:func:`~trackcell.io.convert_annohdcell_to_trackcell`):

* **create_polygons**: Whether to create polygon geometries (default ``True``)

**Method 2** (:func:`~trackcell.io.add_geometries_to_annohdcell_output`):

* **cell_h5ad_path**: Path to final cell h5ad from annohdcell
* **n_jobs**: Number of parallel workers for polygon creation (default ``-1`` = all cores)


Output Format
-------------

The output h5ad contains:

``.X`` — Expression Matrix
    * **Method 1**: Summed counts across all bins
    * **Method 2**: Original counts from annohdcell
    * Sparse CSR matrix format

``.obs`` — Cell Metadata
    ``cellid``, ``object_id``, ``bin_count``, ``geometry`` (WKT), and QC metrics

``.obsm["spatial"]``
    Cell centroid coordinates (mean of bin coordinates), shape ``(n_cells, 2)``

``.uns["spatial"][sample]``
    ``geometries`` GeoDataFrame, ``images``, and ``scalefactors``


Example Workflows
-----------------

Quick Conversion (Method 1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # 1. Run annohdcell (outside Python)
   # $ annohdcell bin2cell --input spaceranger_output --output results

   # 2. Convert to trackcell
   import trackcell as tcl
   adata = tcl.io.convert_annohdcell_to_trackcell(
       bin_h5ad_path="results/b2c_2um.h5ad",
       output_h5ad_path="results/trackcell_format.h5ad",
       sample="sample1"
   )

   # 3. Visualize
   tcl.pl.spatial_cell(adata, library_id="sample1", show=True)

Preserve Exact annohdcell Output (Method 2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # 1. Run annohdcell
   # $ annohdcell bin2cell --input spaceranger_output --output results
   # Produces: b2c_2um.h5ad, b2c_cell.h5ad

   # 2. Add geometries to the final cell h5ad
   import trackcell as tcl
   adata = tcl.io.add_geometries_to_annohdcell_output(
       bin_h5ad_path="results/b2c_2um.h5ad",
       cell_h5ad_path="results/b2c_cell.h5ad",
       output_h5ad_path="results/b2c_cell_with_geom.h5ad",
       sample="sample1",
   )

   # 3. Visualize and analyse
   tcl.pl.spatial_cell(adata, library_id="sample1", show=True)


Restoring Geometries After Reading
----------------------------------

When an h5ad with geometries is saved, the GeoDataFrame is converted to WKT
strings for serialization. After reading back, restore the GeoDataFrame with
:func:`~trackcell.io.restore_geometries`:

.. code-block:: python

   import scanpy as sc
   import trackcell as tcl

   adata = sc.read_h5ad("b2c_cell_with_geom.h5ad")
   adata = tcl.io.restore_geometries(adata)

   # Now spatial visualization works
   tcl.pl.spatial_cell(adata, library_id="sample1")

What it does: converts WKT strings in ``uns["spatial"][sample]["geometries"]``
back to a GeoDataFrame.

.. code-block:: python

   adata = sc.read_h5ad("cell_with_geom.h5ad")
   # type(adata.uns["spatial"]["sample1"]["geometries"])
   # → pandas.core.frame.DataFrame

   adata = tcl.io.restore_geometries(adata, sample="sample1")
   # type(adata.uns["spatial"]["sample1"]["geometries"])
   # → geopandas.geodataframe.GeoDataFrame


Troubleshooting
---------------

**Column 'labels_joint' not found**
    Check that the input h5ad has cell assignment labels, or specify a different
    ``labels_key``.

**'spatial' not found in .obsm**
    Ensure the input h5ad has spatial coordinates and is from annohdcell.

**Column 'object_id' not found (Method 2)**
    Verify the cell h5ad comes from annohdcell's bin_to_cell output and has the
    ``object_id`` column.

**Failed to create polygon for cell X**
    Some cells have problematic bin coordinates; a fallback point geometry is used.

**Cell X has no bins assigned**
    Cells in the cell h5ad without corresponding bins in the bin h5ad
    (e.g. from different runs) will not have geometries added.


See Also
--------

* :doc:`reading` — loading Visium HD and Xenium data
* :doc:`visualization` — ``tcl.pl.spatial_cell`` plotting options
* :doc:`spatial_clustering` — YardCluster expression + spatial domain clustering
