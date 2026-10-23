Installation
============

TrackCell can be installed via pip from PyPI.

PyPI Installation
------------------

.. code-block:: bash

   pip install trackcell -i https://pypi.org/simple

Upgrade to Latest Version
--------------------------

.. code-block:: bash

   pip install --upgrade trackcell -i https://pypi.org/simple

Dependencies
------------

TrackCell requires Python ``>=3.10,<3.13`` and the following packages
(installed automatically via pip):

* scanpy >= 1.9.0
* anndata >= 0.9.0
* numpy >= 1.24.0, < 2.0
* pandas >= 2.1.0
* scipy >= 1.9.0
* scikit-learn >= 1.2.0
* statsmodels >= 0.13.0
* patsy >= 0.5.0
* joblib >= 1.2.0
* matplotlib >= 3.6.0
* numba >= 0.58.0
* geopandas >= 1.1.1
* shapely >= 2.0.0
* imageio >= 2.31.0
* h5py >= 3.0.0
* pyarrow >= 12.0.0 (parquet tables for Xenium / SpaceRanger)
* zarr >= 2.18.0, < 3.0.0

Optional Feature Extras
-----------------------

Some features require additional packages. Install them with the matching
extra:

.. list-table::
   :header-rows: 1

   * - Extra
     - Command
     - Required for
   * - ``clustering``
     - ``pip install "trackcell[clustering]"``
     - YardCluster / scanpy Leiden clustering
       (``leidenalg``, ``python-igraph``)
   * - ``integration``
     - ``pip install "trackcell[integration]"``
     - SCT + RPCA / Harmony / BBKNN integration
       (``annoy``, ``harmonypy``, ``bbknn``)
   * - ``sct``
     - ``pip install "trackcell[sct]"``
     - Optional SCT v2 backends (``pyglmGamPoi``, ``rrng``)
   * - ``ml``
     - ``pip install "trackcell[ml]"``
     - XGBoost-enhanced SCT integration (``xgboost``)
   * - ``interactive``
     - ``pip install "trackcell[interactive]"``
     - Interactive ROI selection in Jupyter (``ipympl``)

For example, the common single-cell / spatial clustering setup is:

.. code-block:: bash

   pip install "trackcell[clustering]"

The equivalent manual installs are also supported:

.. code-block:: bash

   pip install leidenalg python-igraph   # clustering extra
   pip install annoy harmonypy bbknn     # integration extra
   pip install pyglmGamPoi rrng          # sct extra
   pip install xgboost                   # ml extra
   pip install ipympl                    # interactive extra
   pip install scikit-misc              # optional: seurat_v3 HVG on raw counts

``integrate_sct_rpca`` and ``run_sct_integration`` work without ``harmonypy`` /
``bbknn`` unless you use those downstream methods.

``spatial_cluster(preprocess=True)`` works without ``scikit-misc`` (it falls
back to ``flavor="seurat"`` on log-normalized data); install it only if you
want ``seurat_v3`` HVG selection on raw counts.

Development Installation
------------------------

To install TrackCell in development mode:

.. code-block:: bash

   git clone https://github.com/seqyuan/trackcell.git
   cd trackcell
   poetry install --with dev --extras clustering

Run the fast test suite (tests that need external R/Seurat benchmark data are
marked ``slow`` and skipped by default in CI):

.. code-block:: bash

   poetry run pytest -m "not slow"

Benchmark data for the parity tests is read from ``TRACKCELL_BENCHMARK_ROOT``
(default: ``~/trackcell_benchmark``).
