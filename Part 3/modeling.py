import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology

import os
import sys
import csv





#----------------------------
# Define Input Files
#----------------------------
datadirectory = "./inputs/"
topology_file = datadirectory+"topology.txt"

#----------------------------
# Set MC Sampling Parameters
#----------------------------

num_frames = 5000
if '--test' in sys.argv: num_frames=10
num_mc_steps = 10

#----------------------------
# Create movers
#----------------------------

# rigid body movement params
rb_max_trans = 2
rb_max_rot = 0.1

# flexible bead movement
bead_max_trans = 1.00

# Each protein model to move separately
rigid_bodies = [["HerA_A"],["HerA_B"],["HerA_C"],["HerA_D"],["HerA_E"],["HerA_F"]]

################################################

#---------------------------------
# Build the Model Representation
#---------------------------------

# Initialize model
m = IMP.Model()

# Create list of components from topology file
topology = IMP.pmi.topology.TopologyReader(topology_file)
domains = topology.component_list

print('#'*10,domains)


bm = IMP.pmi.macros.BuildModel(m,
                    component_topologies=domains,
                    list_of_rigid_bodies=rigid_bodies)

representation = bm.get_representation()

# add colors to the components
for nc,component in enumerate(domains):
    name = component.name
    sel = IMP.atom.Selection(representation.prot,molecule=name)
    ps = sel.get_selected_particles()
    clr = IMP.display.get_rgb_color(float(nc)/len(domains))
    for p in ps:
        if not IMP.display.Colored.get_is_setup(p):
            IMP.display.Colored.setup_particle(p,clr)
        else:
            IMP.display.Colored(p).set_color(clr)


# Randomize the initial configuration before sampling
# representation.shuffle_configuration(50)






#----------------------------
#  Define Degrees of Freedom
#----------------------------

# Add default mover parameters to simulation
representation.set_rigid_bodies_max_rot(rb_max_rot)
representation.set_floppy_bodies_max_trans(bead_max_trans)
representation.set_rigid_bodies_max_trans(rb_max_trans)

outputobjects = [] # reporter objects (for stat files)
sampleobjects = [] # sampling objects

# Add the movers to the sample and output object lists
outputobjects.append(representation)
sampleobjects.append(representation)







#-----------------------------------
# Define Scoring Function Components
#-----------------------------------

# Here we are defining a number of restraints on our system.
# For all of them we call add_to_model() so they are incorporated into scoring
# We also add them to the outputobjects list, so they are reported in stat files


# Excluded Volume Restraint
# To speed up this expensive restraint, we operate it at resolution 20
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         representation, resolution=20)
ev.add_to_model()
outputobjects.append(ev)


# Crosslinks - dataset 1
# To use this restraint we have to first define the data format
# Here assuming that it's a CSV file with column names that may need to change
# Other options include the linker length and the slope (for nudging components together)

columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None

xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'xlinks.txt',
                                   length=35.0,
                                   slope=0.04,
                                   columnmapping=columnmap,
                                   resolution=1.0,
                                   label="Crosslinks",
                                   csvfile=True)

xl1.add_to_model()             	# crosslink must be added to the model

sampleobjects.append(xl1) 		#crosslink restraint is storing a sampled particle
outputobjects.append(xl1)







#--------------------------
# Monte-Carlo Sampling
#--------------------------
# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    representation,
                                    monte_carlo_sample_objects=sampleobjects,
                                    output_objects=outputobjects,
                                    crosslink_restraints=[xl1],    # allows XLs to be drawn in the RMF files
                                    monte_carlo_temperature=1,
                                    simulated_annealing=False,
                                    simulated_annealing_minimum_temperature=1,
                                    simulated_annealing_maximum_temperature=2.5,
                                    simulated_annealing_minimum_temperature_nframes=200,
                                    simulated_annealing_maximum_temperature_nframes=20,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    number_of_best_scoring_models=100,
                                    monte_carlo_steps=num_mc_steps,
                                    number_of_frames=num_frames,
                                    global_output_directory="output",
                                    atomistic=True)

# Start Sampling
mc1.execute_macro()
