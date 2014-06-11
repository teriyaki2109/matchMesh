*** PLEASE REFRAIN FROM DISTRIBUTING THIS BETA VERSION ***
*** USE IT FOR RESEARCH PURPOSE ***

matchMesh by Tony TUNG (tonytung.org)- 2014 (c)

SUMMARY:
Surface Alignment by Geodesic Mapping [Tung et al., PAMI14] & [Tung et al., CVPR10]
Use provided Makefile for compilation and don't forget to create an ./output/ folder folder before running the executable.
Example: run with ./matchMesh mesh_0080_0.off mesh_0081_0.off or ./matchMesh Horse/16.off Cat/009.off 

DESCRIPTION:
This software can be used to align two surfaces represented by 3D meshes. For example, Mesh1 and Mesh2.
Each vertex of Mesh2 is matched to a vertex of Mesh1 using geodesic mapping.
* The resulting correspondances are given in ./output/matchVec.txt
* The resulting mesh after alignment ./output/reMesh.off has vertices of Mesh1 with mesh connectivity of Mesh2.
* Intermediate results outMesh*.off and outUnique*.off show reference points and ambiguity level respectively (see papers).

DETAILS:
The code calls the library libamrg.a to compute surface geodesics and extremal points.
The code has been tested on 32-bit Debian Linux.
The method implemented here is slightly different from the papers (for the sake of code compactness), but follows the general scheme.
(Differences are in the initial node matching and final mapping processes, which are here simplified and somehow less robust!)

This code works best between consecutive frames of 3D video sequences (i.e., not optimized for wideframe alignment).
Better performance can still be achieved by fine tuning several parameters that were determined empirically (e.g., thresholds ana variables in geomap.h).

For testing, we recommend the public mesh sequences of "breakdance" (e.g., lock, free, etc.) made by Univ. of Surrey,
or the ones processed by INRIA Grenoble which are topologically consistent over time (see [Cagniart et al., ECCV10]).
As we are using geodesics, surface meshes should be manifold (1 connected component) and with sufficient resolution.
All results are ouput in the folder ./output/ , so don't forget to create it before running the executable.


Have fun!

Tony.


--
Tony TUNG <www.tonytung.org>

