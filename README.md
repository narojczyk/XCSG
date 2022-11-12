# (XCSG) X-mer Crystal Structure Generator

Program generates an f.c.c. monomer crystal structure at close packing.
Additional inclusions can be introduced into the structure. The inclusions 
can form channels or layers of arbitraty orientation and diameter or 
thickness respectively. The diameters of particles that form inclusions 
can be controlled. Multiple inclusions in arbitrary orientations are
possible. Monomers (spheres) inside or outside inclusions can be paired into
dimers, forming a *Degenerated Crystal* (DC) of dimers (without inclusions)
or a mixed system of monomer and dimers.

###### Acknowledgements
The project build upon older code by M. Kowalik


## Manifest

Project directory should contain:

* **src**
	Subdirectory with source files, a makefile and a shell script it requires.

## Project Roadmap

- [x] Change sphere types ID's with respect to monomers they form (1 - 8) or
  inclusions (100+).
- [x] Inclusions may be formed by multimers.
- [x] Finite size inclusions (clusters)
- [ ] Add hcp and bcc lattices.
- [ ] Add 3+ multimer particles.

