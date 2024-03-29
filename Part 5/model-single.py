from modeller import *
from modeller.automodel import *
#from modeller import soap_protein_od

env = environ()
a = automodel(env, alnfile='HerA-NurA_align.ali',
              knowns='HerA-NurA', sequence='HerA-NurA_UniProt',
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
a.ending_model = 5

# MD optimization:
# a.md_level = refine.fast

a.make()
