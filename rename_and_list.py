import sys
import glob
import re
import os

path = sys.argv[1]

for filename in glob.glob(path):
	try:
		a = re.split("_", filename)[0]
		b = re.split("_", filename)[3]
		new_filename = a + "_" + b + ".fastq.gz"
		os.rename(filename, new_filename)
	except:
		pass

matches = set()
for filename in glob.glob(path):
    sep = filename.split("_")
    if len(sep) > 1:
        matches.add(sep[0])

for RIL in matches:
	with open("Seq_list.txt", "a") as file:
		file.writelines(RIL + "\n")
