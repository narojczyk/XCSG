# MultiMer Structures (MMS)

Program generates a structure at close packing, composed from monomers,
multimers or their mixtures. Particles can form arbitrary lattices. The
structure may include additional inclusions composed by yet another types of
particles. The inclusions may be take the form of channels, layers (at any
orientation) or any of their combinations.

The project is an extension of the former fccDCgen project


## Manifest

Project directory should contain:

*	**src**

	Subdirectory with source files, a makefile and a shell script it requires.

## Roadmap

* Add hcp and bcc lattice types.
* Add >2 multimer particles.
* Change sphere types ID's with respect to monomers they form (1 - 8) or
  inclusions (100+).
* Inclusions may be formed by multimers.
* Finite size inclusions
* Rename to Multimer Crystal Structures (MCS)
