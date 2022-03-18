# 2022modTFBSs
Public repository of additional scripts used in our expanded alphabet manuscript.

These include:
    - `run_cytomod_TFBS_pipeline.sh`: used along with the below scripts, to construct (K562) modified genomes
        - `computePseudocountForKLequality.py`, `makeRegionBaseFreqPlots.sh`, `output5mC5hmCModBaseBEDsAtThreshold.py`, and `run_hypothesis_tests.sh`.
    - `processDataAndPlot.py`: used to create our main scored data and creates swarm plots
    - `run_matrix-clustering-specific_TFs.sh`: used to run RSAT `matrix-clustering`
    - `gen_treemaps.py`: used to create treemaps, of scored and clustered hypothesis pairs
    - `call_peaks-MACS-uncalibrated-short_reads.sh`: used to call peaks for our CUT&RUN datasets
