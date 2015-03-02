# Ski_Form

Project written in MATLAB to identify location of Xs in scans of ski forms used at Mt. Hood. Uses a combination of image analysis and machine learning to turn paper records into computer records. 

process_allfiles.m
	This script loads images from scans/ and uses image processing to extract subimages. These are then saved to files images_scans.mat & images_stef_scans.mat.

run_machine_learning loads .mat files from earlier and runs 10-fold cross-validation on three machine learning methods: ridge regression, optimal separating hyperplanes (SVM), and feed forward neural networks.

See write_up.pdf for a more in-depth look at the mathematics.