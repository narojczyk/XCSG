# (XCSG) X-mer Crystal Structure Generator

Program generates an f.c.c. monomer crystal structure at close packing.
Additional inclusions can be introduced into the structure. The inclusions 
can form channels or layers of arbitraty orientation and diameter or 
thickness respectively. The diameters of particles that form inclusions 
can be controlled. Multiple inclusions in arbitrary orientations are
possible. Monomers (spheres) outside inclusions can be paired into dimers, 
forming a *Degenerated Crystal* (DC) of dimers, in no inclusions are set, or
a mixed system of monomer and dimers. 

###### Acknowledgements
The project build upon older code by M. Kowalik


## Manifest

Project directory should contain:

* **src**
	Subdirectory with source files, a makefile and a shell script it requires.

## Project Roadmap

- [ ] Change sphere types ID's with respect to monomers they form (1 - 8) or
  inclusions (100+). *in progress*
- [ ] Inclusions may be formed by multimers. *in progress*
- [ ] Finite size inclusions. *in progress*
- [ ] Add hcp and bcc lattices.
- [ ] Add 3+ multimer particles.

