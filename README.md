# SpectralFFS
Matlab Code associated to "Multi-color fluorescence fluctuation spectroscopy in living cells via spectral detection", https://doi.org/10.1101/2020.12.18.423407, written together with Salvatore Chiantia

Multi-color Fluorescence Fluctuation Spectroscopy Code (for MATLAB) – Instructions

System requirements: MATLAB 2015b or newer, at least 16 GB RAM system (loading spectral .tif files
requires sufficient RAM)


Calibration of spectra (needed for both SFSCS and RSICS) 

Extract and pool spectra from single-species images

• Run script spectra_one_component_memb_czi_batch.m (for membrane-localized FPs) or spectra_one_component_cyto_czi_batch.m (for cytosolic FPs) to load all files from a single directory
• Pool spectra from multiple species and generate spectral_patterns.txt file by running scripts spectral_patterns_2ch.m, spectral_patterns_3ch, or spectral_patterns_4ch depending on the number of FP species


SFSCS analysis pipeline 

Export

• Export acquired .czi line scan files as .tif in Zen Software. Save one .tif file per three spectral bins (saved as RGB triple), e.g. 8 files for 23 bins with the last file containing only 2 bins, all to one directory (hence one directory for each measurement). Save files with index at the beginning of filename (e.g. 01_measurement1, 02_measurement1 etc.)

Analysis

• Set acquisition parameters in header of SFSCS_analysis_2species.m/ ...3species.m/ ..._4species.m and run file
• Select directory where exported files are located
• Select directory where output files shall be saved
• Do rectangular selection of membrane section in kymograph. Make sure to start and close the
rectangle outside the plot window
• Wait for import to be finished (filenames are printed in command window). This may take a
bit of time
• Do polygonal selection around membrane (again close outside the plot window). Exclude
bright intracellular structures (e.g. vesicles)
• If a background correction is needed (backgroundcorrection=1), select a polygonal ROI on
enclosing the background signal in the kymograph (on one side of the membrane, usually the
intracellular side)
• The lateral alignment of the kymograph proceeds and the selected region is displayed
• Load spectral_patterns.txt file. The spectral decomposition will now proceed. If you selected
binning, time binning will be performed
• The average spectral fractions are determined and printed in the command window and dialog
box. You can now choose to continue the analysis or terminate, for example if not all species
are present in this measurement
• The bleaching correction is performed. As part of it, a GUI pops up to removed intensity
segments (e.g. bright peaks) that may distort the double-exponential fit. Remove those segments if occurring and press Done afterwards. The fit will be determined for each species and the bleaching correction applied
• Raw correlation functions are plotted
• The correlation functions are now calculated in segments and plotted as overlay
• A GUI pops up allowing you to remove segments from the analysis. Evaluate the average CFs (solid lines) and scatter of points (CFs of the segment) as well as intensity trace (fluctuating around average or long-term oscillations/ deviations?) to keep or remove segments with associated buttons. Adjust y axis limits to re-scale the axes (which may be helpful if different species are present in different amounts, thus showing different amplitudes). Iteratively go through the segments and confirm with Done button at the end of the inspection. For stable measurements, very few segments need to be removed. For an illustration of the procedure, watch the video associated to Dunsing&Chiantia, Jove 2018
• The final fitting is now performed
• Final graphs and output files (e.g. fit parameters) are saved in specified directory


Pooling of SFSCS results from multiple analyzed measurements

• Run scripts SFSCS_parametertxtfiles_pooled_2species.m/ ...3species.m/ ...4species.m to pool fit parameters from multiple analyzed files. Load the files from the same directory
• The script generates a table with header saved as a .txt file, containing Ns, rel.cc.’s, tau’s, Brightness values, and bleaching fractions
RSICS and TRICS analysis pipeline Decomposition
• Run spectral_decomposition_RSICS_czi_batch.m to load all files and decompose three- or four- species RSICS .czi image stacks, using the spectral_patterns.txt file containing the spectral patterns (load this file once at the beginning after specifying the directory where decomposed image stacks are saved)
Analysis
• Run script Image_selection.m to select one (or several) sub-regions of acquired images for later analysis. Imaging parameters can be adjusted in this script. The output is one or more “batch” files that will be read by the following script. Image selection is performed via the roipoly MATLAB function.
• Run ARICCS_SpectralDecomposition.m (ARICCS_SpectralDecomposition_4Channels for 4- species analysis) and select the directory where the batch files are stored. This script will automatically analyze the files and save the results as .txt files and figures.
• For TRICS analysis, the same Image_selection.m script can be used before TRICS.m. The calculation of the 4d matrix containing the correlation data is performed at the moment via a slow algorithm by the function manualCCC.m. The second part of TRICS.m is not strictly needed and presents only several alternatives for data plotting.
