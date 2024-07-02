# A Radial Cell Model for Classification of Protein Localization

This repository aims to perform classification of signal in fluorescent microscopy images by assuming a spherical shape approximation of cells and then basing the classification on the localization of each of the proteins of interest in relation to the sphere that represents a cell and its membrane boundary. The spherical model is derived from the fluorescent signal given by nuclei. Once nuclei are segmented via a difference of gaussians and threshold, centroids are computed for resulting connected components (CCs) and these centroids represent the center of the nuclei detected per slice. Although detections are performed in 2 dimensions as slices are sequentially processed with the same parameters, the detections on each slice can be assigned a z coordinate that equates to the z slice that they were found on. This allows us to construct a 3D point cloud that then has the z dimension scaled to ensure it represents an isotropic point cloud. k nearest neighbors with k=1 is computed to find the average distance between adjacent nuclei. This is approximately twice the radius of any given cell nucleus. Then, Density-based spatial clustering of applications with noise (DBSCAN) will use the implied radius as the search radius to cluster the point detections in 3D.


While the algorithm described above can be generalized to any similar fluorescent microscopy experiments in tissue culture ranging from monolayers to organoids, there are some system specific corrections that need to be accounted for.
The Stage Scanning confocal microscope developed by Dr. Matthew Parent in Dr. Abhishek Kumar's Lab is prone to producing a minor artifact in images along the x dimension (columns) that results from the illumination profile of the Powell lens used in the illumination pathway. Additionally, the system multiplexes two channels at a time. Given that a user might be interested in imaging multiple channels, experiments often involve using a single channel as a reference channel while the channels of interest are then multiplexed with the reference channel. This ensures that as image volumes are taken, we can ensure that there is information to use as reference to register the channels of interest. Here, we use phase correlation with registration code written by Dr. Matthew Parent to find the shifts needed in the y and z dimensions .


The folders herein contain codes that will perform:
1. Preprocessing of the 3D image volumes
  * Normalization: Reduction of powell lens illumination profile artifact, i.e. beam correction
  * Registration: Channels that were multiplexed with a reference channel will be registered by shifting the second channel to have the same coordinates


2. Analysis of the images
  * Nuclei detection: Difference of Gaussians, threshold, centroids of CCs
  * Sphere Model: find *r* from 2D nearest neighbors, scale isotropically, 3D clustering of pointcloud with *r*, average cluster points to find 3D localization
  * Protein localization classification: resample images to be isotropic, segment protein channels, determine if signal is within distance *r* to classify cells as expressing the protein inside or not


3. Visualization of Analysis
  * Stacked Bar plots
  * Number of voxels per nucleus
  * Intensity of signal
    - Size invariant: Average voxel intensity per nucleus - Normalization by number of protein voxels
    - Size dependent: Integrated voxel intensity per nucleus - does not normalize by total number of protein channel voxels inside the radius


Future directions of this algorithm might consider the radial distribution function (RDF). We could do this in an efficient manner by leveraging some HPC with sufficient RAM. We would need to compute the distance matrix for all pixels in the binary image volume in relation to the 3D nuclei centroids. From this, we can then do a voxel distance based RDF to show the probability of finding a protein within a given distance away from the cell. Given that we have ~1400+ cells detected in our images, this would enable statistically sound arguments with additional discretization as the framework for a more continuous classification regime. Additionally, this could also scale according to computed voxel intensities to then also know what the intensities are across different distances.

The above improvement on the classification model is a matter of moving *r* from $$r_{1} to r_{2}$$ by some small step *dr* where *dr* is simply the smallest possible Euclidian distance, *d* in the data. In this case, our images will have voxel size [dx dy dz] so *d* would be sqrt(dx^2 + dy^2 + dz^2).

Lastly, one might consider bootstrapping the data. One might also want to consider where the sample is in focus to guide the boundaries of the bootstrapping. This bootstrapping of a specific subset may introduce some bias towards detection of information in planes that are in focus.
