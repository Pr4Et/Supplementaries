These Matlab scripts were used to process 4D-STEM data for the paper “Shadow montage and cone-beam reconstruction in 4D-STEM tomography”. They are intended either for hdf5 datasets acquired by ARINA pixelated detector using the Savvyscan software or for raw binary files from other labs. The scripts with the name ELA refer to another pixelated detector and are intended to process prz files acquired using CEOS Panta Rhei software.    
Processing generally follows several stages.
Once, for determining the direction of the scan engine in comparison to the pixelated detector, we run ShadowMontage_checkDirections.m.
For getting a sense of the data we run ShadowMontage_Parallel_Sync. The result is a stack of shadow montage for different synchronization numbers.
Then ShadowMontage_tiltseries.m generates shadow montages both in 2D (*_Tiltview.mrc) and in 3D (*_3dmontage_S*.mrc) for every tilt view. The choice of sych number may be either based on previous step or automatically (-1) based on image shifts in DPC analysis.
The *Tiltview.mrc files are used to find the tilt series alignment. We further relax the boundaries using ForAretomo.m and run in Aretomo (versions 1 or 2) to acquire the alignment table (*.aln) and a 3D SART reconstruction if needed. 
Based on the alignment data in *aln file and 3D shadow montage (*_3dmontage_S*.mrc)  a 3D backprojection from the series is generated in Backproject_align3d_afterAretomo.m.
The 3D backprojection is processed further by 3D deconvolution (see details in Supplementaries/Analysis_4DSTEM_tiltseries at main · Pr4Et/Supplementaries · GitHub).
A 3D shadow montage can be also generated as a 3D cone-beam backprojection (using Astra Toolbox), which  is demonstrated in cone_rec3d.m. For processing prz file we convert them to hdf5.  
The scripts with CTF in the name are versions with CTF correction (only sign flipping of spectral components, without inversion for negative defocus). These separate scripts require more information on the settings.

Installation in Matlab: Requires @MRCImage library from MatTomo, PEET project: https://bio3d.colorado.edu/imod/matlab.html and hdf5 libraries from https://www.hdfgroup.org/download-hdf5/.
Installation on Windows: The executable software requires only Matlab runtime of 2024a version installed on Windows (see executables folder).
