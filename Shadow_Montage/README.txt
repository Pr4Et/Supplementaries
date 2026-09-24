The Matlab scripts were used to process 4D-STEM data for the paper “Shadow montage and cone-beam reconstruction in 4D-STEM tomography” (https://doi.org/10.1093/mam/ozaf126). They are intended either for hdf5 datasets acquired by ARINA pixelated detector using the Savvyscan software or for raw binary files from other labs. The scripts with the name ELA refer to another pixelated detector and are intended to process prz files acquired using CEOS Panta Rhei software.    

The python script, shadow_montage_tiltseries.py, demonstrates the fastest and most straightforward generation of shadow montage tilt series based on metadata JSON file created with https://github.com/elbaum-lab/Session-JSON. The tiltseries can be reordered according to the list of angles in the JSON file and then aligned and transformed to a tomogram using IMOD, for example.
Using the matlab script for unsupported cases requires first determining the direction of the scan engine in comparison to the pixelated detector by running ShadowMontage_checkDirections.m. 
For generating a stack of different synchronizations (for example when the defocus is unknown) we run ShadowMontage_Parallel_Sync_imageCTF.m.
The first verion (more experimental), ShadowMontage_tiltseries_imageCTF.m, generates shadow montages both in 2D (*_Tiltview.mrc) and in 3D (*_3dmontage_S*.mrc) for every tilt view. The choice of sych number may be either based on previous step or automatically (-1) based on image shifts in DPC analysis.
The *Tiltview.mrc files are then used to find the tilt series alignment. ForAretomo.m adjusts the boundaries and Aretomo (versions 1 or 2) is used to generate the alignment table (*.aln) and a 3D SART reconstruction if needed. 
Based on the alignment data in *aln file and 3D shadow montage (*_3dmontage_S*.mrc)  a 3D backprojection from the series is generated in Backproject_align3d_afterAretomo.m. This long process in most cases may be replaced by the more efficient version 2 or the python script.
Note that the 3D backprojection from the longer processing pipeline should be further processed by 3D deconvolution (see details in Supplementaries/Analysis_4DSTEM_tiltseries at main · Pr4Et/Supplementaries · GitHub).
A 3D shadow montage can be also generated as a 3D cone-beam backprojection (using Astra Toolbox), which  is demonstrated in cone_rec3d.m. For processing prz file we convert them to hdf5.  
The scripts with CTF in the name are versions with CTF correction (only sign flipping of spectral components, without inversion for negative defocus). These separate scripts require more information on the settings.

Installation in Matlab: Requires a modified MatTomo @MRCImage library (https://github.com/Pr4Et/Supplementaries/tree/main/Analysis_4DSTEM_tiltseries/MatTomo) (based on PEET project: https://bio3d.colorado.edu/imod/matlab.html).
Requires hdf5 libraries from https://www.hdfgroup.org/download-hdf5/.
For versions of Matlab above 2024 use the supplied uigetfile.m.
Running exe files in Windows: The executable software requires only Matlab runtime of 2024a version installed on Windows (see executables folder).

Updates

ShadowMontage_tiltseries_ver2.m  - Use modified CTF correction in overfocus to prevent mixed contrast contributions between the shadow overall illumination and the shadow patterns.