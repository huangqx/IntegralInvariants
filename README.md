# This software provides a simple implementation of the volume descriptor described in the following paper:
Integral Invariants for Robust Geometry Processing. Helmut Pottmann, Johannes Wallner, Qixing Huang, and Yongliang Yang. CAGD' 2009.

Usage:
Para.dessDim = 6; // number of descriptors
Para.rMin = 0.05; // minimum radius of the ball
Para.rMax = 0.1;  // maximum radius of the ball 
Para.gridRes = 0.01; // grid resolution for discretzing the   
                        volume integral

// Load the mesh

Mesh = read_from_obj('Data/bunny.obj');

// Compute the volume descriptors

vertexDess = volume_invariant_3d(Mesh, Para);

// Save the descriptors

save_descriptors(vertexDess, 'Data/bunny_dess.txt');
