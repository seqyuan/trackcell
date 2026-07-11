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

TrackCell requires the following packages:

* Python >= 3.10
* scanpy >= 1.9.0
* geopandas >= 1.1.1
* pandas >= 2.1.0
* shapely >= 2.0.0
* imageio >= 2.31.0
* numpy >= 1.24.0
* scipy (for spatial distance calculations)

These dependencies are automatically installed when installing TrackCell via pip.

Optional Dependencies (YardCluster)
------------------------------------

Spatial clustering with :doc:`usage/spatial_clustering` additionally requires:

.. code-block:: bash

   pip install trackcell[clustering]

Or manually:

.. code-block:: bash

   pip install leidenalg python-igraph

For multi-sample Harmony batch integration (``integrate='joint'`` in
:doc:`usage/spatial_clustering`):

.. code-block:: bash

   pip install harmonypy

Batch expression integration (SCT + RPCA / Harmony / BBKNN) — see
:doc:`usage/batch_integration`:

.. code-block:: bash

   pip install annoy              # RPCA anchors (required for integrate_rpca)
   pip install harmonypy          # SCT + Harmony
   pip install bbknn                # SCT + BBKNN
   pip install pyglmGamPoi          # native SCT v2 (or R glmGamPoi via method='r')

``integrate_sct_rpca`` and ``run_sct_integration`` work without ``harmonypy`` /
``bbknn`` unless you use those downstream methods.

Development Installation
------------------------

To install TrackCell in development mode:

.. code-block:: bash

   git clone https://github.com/seqyuan/trackcell.git
   cd trackcell
   pip install -e .

