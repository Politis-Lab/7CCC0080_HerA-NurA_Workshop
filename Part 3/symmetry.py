import os, csv
import numpy as np
from pymol import cmd, stored

w_dir = "./output_atomistic/"        	# path to where the files are

# Get all models
models = os.listdir(w_dir)
print "\n%s files found in %s" % (len(models), w_dir)
print "Calculating RMSD symmetry score.. "

symmetryList = []
for m in models:
    print "  Loaded file "+m

    cmd.load(w_dir+m,m[:-4])

    # Find number of chains in the file
    stored.chains = []
    cmd.iterate(str(m[:-4]), 'stored.chains.append(chain)')

    stored.chains = sorted(list(set(stored.chains)))
    rev_chains = sorted(stored.chains,reverse=True)

    # Duplicate loaded object m and call n
    cmd.create("dummy",str(m[:-4]))
    print "create dummy, %s;" % str(m[:-4])

    # Change first chain to zz
    # cmd.alter("dummy & chain %s" % stored.chains[0], "chain='Z%s'" % stored.chains[0])
    # print "alter dummy & chain %s, chain='Z%s'" % (stored.chains[0],stored.chains[0])

    # Change each dummy chain to Zx first then change to reverse chain list
    for index, ch in enumerate(stored.chains):
        cmd.alter("dummy & chain %s" % ch, "chain='Z%s'" % ch)
        print "alter dummy & chain %s, chain = 'Z%s';" % (ch,ch)

    # Loop over remaining chains and change
    for index, ch in enumerate(stored.chains):
        cmd.alter("dummy & chain Z%s" % ch, "chain='%s'" % rev_chains[index])
        print "alter dummy & chain Z%s, chain='%s';" % (ch,rev_chains[index])

    # Now change chain zz to the first chain in the reverse list
    # cmd.alter("dummy & chain ZZ", "chain='%s'" % rev_chains[0])
    # print "alter dummy & chain ZZ, chain='%s'" % rev_chains[0]

    # Now align chain A with chain A
    rmsd_perChain = []
    for ind, pc in enumerate(stored.chains):
        cmd.align("%s & chain %s" % (m[:-4],pc), "dummy & chain %s" % pc)
        print "align %s & chain %s, dummy & chain %s;" % (m[:-4],pc,pc)

        rmsd = cmd.rms_cur("dummy",str(m[:-4]))
        print "rms_cur dummy, %s" % str(m[:-4])
        print "rmsd: %s angstroms" % rmsd

        rmsd_perChain.append(rmsd)

    print "rmsds for each chain comparison: %s" % rmsd_perChain

    sumRMSD = np.sum(rmsd_perChain)
    print "Sum per-chain rmsd: %s" % sumRMSD

    cmd.delete(str(m[:-4]))
    print "delete %s" % str(m[:-4])

    cmd.delete("dummy")
    print "delete dummy"

    symmetryList.append([m[:-4], sumRMSD])

# with open('RMSD_symmetry.csv', 'a') as f:
# 		writer = csv.writer(f, dialect='excel')
# 		writer.writerows(zip(symmetryList))

with open('score_symmetry.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(symmetryList)
