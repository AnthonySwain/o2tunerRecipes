# Brief Overview of the optimisations for the O2 build

I would reccomend using the "ParticleSpecificMaps" above, which uses the O2pdgspecific build, for the following reasons: <br>
a: Significantly increased functionality with PDG specific maps <br>
b: Can assign multiple voxel maps at the same time through one command <br>
c: I developed the more useful framework for the optimisations with the O2pdgspecific build <br>
d: This is the full extent of the readme that I am going to write... <br>

Apart from the genetic algorithm...

# Genetic Algorithm
This is something that is not implemented in the O2pdgspecific build. This is a ML algorithm but is very barebones and O2Tuner provides a much more sophisticated genetic algorithm. However, this has the one advantage that it is significantly faster than O2Tuner... currently. Either way, I wouldn't reccomend its use... but this is just to know I made this once. 