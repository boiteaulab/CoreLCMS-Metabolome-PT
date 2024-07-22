# CoreMS LC Metabolome pipeline
This data processing pipeline annotates a liquid chromatography mass spectrometry (LCMS) metabolome feature list (e.g. MSDial, MZmine) with a molecular formula library based on the analysis of a pooled sample by liquid chromatography with ultrahigh resolution mass spectrometry (LC-21T FT-ICR MS).

CoreMS_LC_internalstd_pt.py:
Python Script for assessing the peak area and retention time of quality control internal standards across a data set. 

CoreMS_LC_assignments_pt.py:
Python script for generating a molecular formula library (.csv) from an LCMS data set (Thermo .RAW file)
Molecular formula library generation requires the CoreMS package.
https://github.com/EMSL-Computing/CoreMS

CoreMS_LC_metabolome_pt.py:
Python script for annotating a metabolome feature list with a molecular formula library (.csv)

CoreMS_LC_plots_pt.py:
Python script for generating figures

siloxanes_pos.ref:
reference files for internal mass calibration peaks. 

Example MS data is available in the MassIVE repository under accession #MSV000094385. 
[https://massive.ucsd.edu/](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000094385)
