.. list-table:: SCT stepwise parity (ref seed, native Python vs R)
   :header-rows: 1
   :widths: 12 8 10 10 10 10

   * - Sample
     - Native HVG J
     - R fit→Py HVG J
     - θ corr (step1)
     - RV Spearman
     - Residual corr
   * - GSM8779707
     - 0.9802
     - 1.0000
     - 0.959
     - 1.000
     - 1.000
   * - GSM8779708
     - 0.9512
     - 1.0000
     - 0.931
     - 0.997
     - 1.000
   * - GSM8779709
     - 0.9822
     - 1.0000
     - 0.921
     - 0.999
     - 1.000

.. list-table:: HVG randomness (6 R seeds × 3 samples)
   :header-rows: 1
   :widths: 30 20

   * - Metric
     - Value
   * - R inter-seed HVG Jaccard (median)
     - 0.9443
   * - R inter-seed HVG Jaccard (min – max)
     - 0.8697 – 0.9646
   * - Python vs R (ref seed) HVG Jaccard (median)
     - 0.9802
   * - Py median within R seed min–max
     - False
   * - Py median within R inter-seed IQR
     - False
