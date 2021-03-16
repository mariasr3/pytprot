import os
import sys
import re

file_re = re.compile("(\w+)\.(\w+)_1(\w+)\..*\.pdb")

for filename in os.listdir("."):
    m = file_re.match(filename)
    if not m:
        print ("Skipping %s" %filename)
        continue

    protein_name = m.group(1)
    pdb_name = m.group(2)
    chain = m.group(3)

    chains = set([ line[21] for line in open(filename) if line.startswith("ATOM") and line[21]!=" "])
    chains.discard(chain)
    dnachains = "".join(sorted(chains))
    command = f"mv {filename} {protein_name}.DNA.{pdb_name}_{chain}_{dnachains}.pdb"
    os.system(command)


