import os, csv
from pymol import cmd, stored
import pandas as pd

w_dir = "./output_atomistic/"        	# path to where the files are
crosslinking = pd.read_csv("./HerA-HerA_crosslinks.csv")

# Get all models
models = os.listdir(w_dir)
print "\n%s files found in %s" % (len(models), w_dir)
print "Calculating cross-linking distance score.. "

distList = []
for m in models:
    print "  Loaded file "+m

    cmd.load(w_dir+m,m[:-4])

    distances = []
    for id, row in crosslinking.iterrows():
        dst = cmd.distance('dst_%s' % id, '/%s//%s/LYS`%s/NZ' % (m[:-4],row['Chain1'],row['Resi1']), '/%s//%s/LYS`%s/NZ' % (m[:-4],row['Chain2'],row['Resi2']))
        print "distance dst_%s, '/%s//%s/LYS`%s/NZ', '/%s//%s/LYS`%s/NZ'" % (id, m[:-4],row['Chain1'],row['Resi1'], m[:-4],row['Chain2'],row['Resi2'])

        distances.append(dst)

    crosslinking['Distance'] = distances

    SXL = []
    for id, row in crosslinking.iterrows():
        if row['Distance'] < 35:
            SXL.append(0)
        else:
            score = (35-row['Distance'])**2
            # print score
            SXL.append(score)

    sumScore = sum(SXL)

    distList.append([m[:-4],sumScore])

    cmd.delete(str(m[:-4]))
    print "delete %s" % str(m[:-4])

with open('score_crosslinks.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(distList)
