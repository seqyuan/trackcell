.. list-table:: Integration SCT stepwise (ref seed, vs R SCTransform)
   :header-rows: 1
   :widths: 14 12 10 24

   * - Sample / merged
     - Metric
     - Value
     - Notes
   * - GSM8779707
     - Native HVG J
     - 0.9182
     - native vst vs R SCTransform
   * - GSM8779707
     - R fit→Py HVG J
     - 0.9500
     - residual corr 0.994
   * - GSM8779708
     - Native HVG J
     - 0.9133
     - native vst vs R SCTransform
   * - GSM8779708
     - R fit→Py HVG J
     - 0.9582
     - residual corr 0.993
   * - GSM8779709
     - Native HVG J
     - 0.9255
     - native vst vs R SCTransform
   * - GSM8779709
     - R fit→Py HVG J
     - 0.9324
     - residual corr 0.991
   * - merged
     - Features J (native)
     - 0.9231
     - native SCT SelectIntegrationFeatures
   * - merged
     - Prep corr (R features)
     - 0.9999
     - native SCT + forced R feature list

.. list-table:: RPCA stepwise parity (ref seed, R input → Python step)
   :header-rows: 1
   :widths: 14 10 10 10 10

   * - Step
     - Metric
     - Value
     - Notes
     -
   * - 03 features
     - Jaccard
     - 1.0000
     - R SCT → Py SelectIntegrationFeatures
     - 
   * - 04 prep
     - Gene corr (median)
     - 1.0000
     - R SCT + R features → Py Prep
     - 
   * - 05 anchors
     - Local pair Jaccard
     - 0.8679
     - count ratio 1.000
     - 
   * - 06 integrate
     - Gene corr (median)
     - 0.9595
     - R prep + R anchors → Py IntegrateData
     - 
   * - 06 integrate
     - Gene corr (Py anchors)
     - 0.9594
     - R prep + Py anchors → Py IntegrateData
     - 

.. list-table:: Python vs Seurat quick comparison (ref seed)
   :header-rows: 1
   :widths: 28 12 12

   * - Metric
     - Stepwise
     - E2E
   * - Features Jaccard
     - 1.0000
     - 0.9286
   * - Prep gene corr
     - 1.0000
     - 1.0000
   * - Anchor local-pair Jaccard
     - 0.8679
     - 0.7196
   * - Integrate gene corr
     - 0.9595
     - 0.9589

.. list-table:: RPCA end-to-end native Python vs R (ref seed)
   :header-rows: 1
   :widths: 22 12

   * - Metric
     - Value
   * - SCT HVG Jaccard (median per batch)
     - 0.9175
   * - Integration features Jaccard
     - 0.9286
   * - Shared integration features
     - 1926
   * - Prep gene corr (aligned columns)
     - 1.0000
   * - Prep gene corr (shared features only)
     - 1.0000
   * - Anchor local pair Jaccard
     - 0.7196
   * - Integrated gene corr (aligned)
     - 0.9589
   * - Integrated gene corr (shared features)
     - 0.9774
