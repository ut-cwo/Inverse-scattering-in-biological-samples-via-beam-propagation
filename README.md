# Inverse-scattering-in-biological-samples-via-beam-propagation
This repository contains the MATLAB reconstruction code for MSBP-based inverse-scattering to recover the sample's three-dimensional refractive index, as used in the following publication:

**Kim, Jeongsoo, et al. "Inverse-scattering in biological samples via beam-propagation." bioRxiv (2025).**

**Abstract**:
Multiple scattering limits optical imaging in thick biological samples by scrambling sample-specific information. Physics-based inverse-scattering methods aim to computationally recover this information, often using non-convex optimization to reconstruct the scatter-corrected sample. However, this nonconvexity can lead to inaccurate reconstructions, especially in highly scattering samples. Here, we show that various implementation strategies for even the same inverse-scattering method significantly affect reconstruction quality. We demonstrate this using multi-slice beam propagation (MSBP), a relatively simple nonconvex inverse-scattering method that reconstructs a scattering sample’s 3D refractive-index (RI). By systematically conducting MSBP-based inverse-scattering on both phantoms and biological samples, we showed that an amplitude-only cost function in the inverse-solver, combined with angular and defocus diversity in the scattering measurements, enabled high-quality, fully-volumetric RI imaging. This approach achieved subcellular resolution and label-free 3D contrast across diverse, multiple scattering samples. These results lay the groundwork for robust use of inverse-scattering techniques to achieve biologically interpretable 3D imaging in increasingly thick, multicellular samples, introducing a new paradigm for deep-tissue computational imaging.

# Experimental dataset

The experimental dataset for running the code can be downloaded from the link below. This dataset contains angular scattering measurements (complex-field measurements) of a scattering phantom, C.elegans, intestinal organoids and zebrafish embryo. Parameters of the angle-scanning imaging system and parameters used for MSBP-based inverse scattering are also included.

**Download Link:** https://dataverse.tdl.org/dataset.xhtml?persistentId=doi:10.18738/T8/YR1ONC

# Running the code

**To reproduce the MSBP-based inverse-scattering results for scattering samples, please follow the steps below.**

1. Please download the Recon_corefunction folder, which contains the core functions required to run the reconstruction code.

2. Download the data from the link provided in the Experimental Dataset section and save it in the same directory as the Recon_corefunction folder.

3. The scripts Recon_main.m is used to reconstruct the experimental datasets. Place that in the same directory as the Recon_corefunction folder, then run the scripts in MATLAB to reproduce the reconstruction results.

4. Since the code contains comments explaining each section, it is recommended to execute the code step by step and read through the related comments when running the code for the first time.

<p align="center">
<img src="ImageFolder/Celegans.jpg" width="700"/>
</p>

<p align="center">
<strong> Figure 1: </strong> C.elegans inverse-scattering results 
</p>

<p align="center">
<img src="ImageFolder/Organoids.jpg" width="700"/>
</p>

<p align="center">
<strong> Figure 2: </strong> Intestinal organoids inverse-scattering results 
</p>

<p align="center">
<img src="ImageFolder/Zebrafish.jpg" width="700"/>
</p>

<p align="center">
<strong> Figure 2: </strong> Zebrafish embryo inverse-scattering results 
</p>
