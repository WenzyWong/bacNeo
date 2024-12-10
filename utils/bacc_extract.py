import pandas as pd
import os
import sys

args = sys.argv[:]
path = args[1]
sample = args[2]

dir = os.listdir( path )
os.chdir( path )
mpa_path = path + '/' + sample + '.mpa'
mpa = pd.read_table(mpa_path, header = None, names = ['Taxonomy', 'Abundance'])
# Extract bacteria taxonomy
mpa_sub = mpa[ mpa.Taxonomy.str.contains("Bacteria") ]

fo_f = open(sample + "_family.txt", "w")
fo_g = open(sample + "_genus.txt", "w")
fo_s = open(sample + "_species.txt", "w")

for i in mpa_sub.index:
    taxon = mpa_sub['Taxonomy'][i]
    taxon = str(taxon)
    count = mpa_sub['Abundance'][i]
    count = str(count)
    if ("f__" in taxon) and ("g__" not in taxon) and ("s__" not in taxon):
        fo_f.write(taxon + "\t" + count + "\n")
    if ("g__" in taxon) and ("s__" not in taxon):
        fo_g.write(taxon + "\t" + count + "\n")
    if "s__" in taxon:
        fo_s.write(taxon + "\t" + count + "\n")

fo_f.close()
fo_g.close()
fo_s.close()

