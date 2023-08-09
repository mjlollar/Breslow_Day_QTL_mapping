#### Author: Matthew Lollar
#### QTL scan framework using the Breslow-Day test statisitic
#### UNIDIRECTIONAL
#### Last update: Aug 8, 2023
#### mjlollar1@gmail.com

### Required Python libraries
import argparse
import numpy as np
import pandas as pd

### Input arguments
parser = argparse.ArgumentParser(description='Breslow-Day based QTL mapping script (get bd cell values), last update: 07/2023')
parser.add_argument('--i', help='Input File (use full path if not in cwd)', required=True, type=str)
parser.add_argument('--o', help='Output File Prefix', required=True, type=str)
parser.add_argument('--s', help='Sterile list File Name (use full path if not in cwd)', type=str, required=True)
parser.add_argument('--f', help='Sterile list File Name (use full path if not in cwd)', type=str, required=True)
parser.add_argument('-ff', help
parser.add_argument('--uset', help='Specify Y or Mito scan (options: "M" or "Y")', type=str, required=True, choices=('M', 'Y'))
parser.add_argument('--focal', help='Run a subset of unidirectional scans (options: 0 (Fr) or 2 (Zi))', type=int, choices=(0,2), required=True)
args = parser.parse_args()

### Chromosome window boundaries (zero-indexed, unidirectional chromosomes last two windows, mito then Y)
## Adjust as needed, current based on Matt's 50kb window format
X_end = 545
Chr2_end = 1524
Chr3_end = 2579

#Chromosome window ranges
len_X = list(range(0, X_end))
len_2 = list(range(X_end, Chr2_end))
len_3 = list(range(Chr2_end, Chr3_end))

mito_win = [2579] #Adjust as needed
y_win = [2580] #Adjust as needed

### Read in data, components are specific to current input for MJL as of 7/31/23
df = pd.read_csv(args.i, sep=',') #MJL pipeline format input uses tab-delim
df.drop(df.columns[[0,1,2]], axis=1, inplace=True) # not necessary if your input contains only genotypes
with open(args.s) as list_steriles:
	sterile_ids = [line.rstrip() for line in list_steriles]
with open(args.f) as list_fertiles:
	fertile_ids = [line.rstrip() for line in list_fertiles]
list_steriles.close() #sanity
list_fertiles.close()

# sanity remove empty characters in sterile/fertile list (may occur if txt files end in newline)
if '' in sterile_ids:
	while('' in sterile_ids):
		sterile_ids.remove('')
if '' in fertile_ids:
	while('' in fertile_ids):
		fertile_ids.remove('')

# index list for loop
index_list = df.columns.values.tolist()

### Initialize Breslow-Day cell count lists
###                   W2F                      W2NF
###              W1F   W1NF                 W1F    W1NF
###        S    bd1     bd3              S   bd5     bd7
###        F    bd2     bd4              F   bd6     bd8

### Unidirectional scan
### Initialize Breslow-Day cell count lists
#forward scans
bd_1 = []
bd_2 = []
bd_3 = []
bd_4 = []
bd_5 = []
bd_6 = []
bd_7 = []
bd_8 = []

### Scan function
def BD_scan(chr_1, chr_2, focal):
	if focal == 0:
		focal_1 = 0
		focal_2 = 2
	else:
		focal_1 = 2
		focal_2 = 0

	for w1 in chr_1:
		for w2 in chr_2:
			b1=0
			b2=0
			b3=0
			b4=0
			b5=0
			b6=0
			b7=0
			b8=0
			for index in index_list:
				if df.at[w1, index] == -999 or df.at[w2, index] == -999:
					pass #skip comparison if no call at either window
					print("skipping no call....")
				elif df.at[w1, index] == focal_1: #W1F
					if df.at[w2, index] == focal_2: #W2F
						if index in sterile_ids:
							b1 += 1
						elif index in fertile_ids:
							b2 += 1
						else:
							raise Exception("Error in calculations; Sterile/Fertile indices are incorrect") #sanity error catching
					else: #W2NF
						if index in sterile_ids:
							b5 += 1
						elif index in fertile_ids:
							b6 += 1
						else:
							raise Exception("Error in calculations; Sterile/Fertile indices are incorrect")
				else: #W1NF
					if df.at[w2, index] == focal_2: #W2F
						if index in sterile_ids:
							b3 += 1
						elif index in fertile_ids:
							b4 += 1
						else:
							raise Exception("Error in calculations; Sterile/Fertile indeces are incorrect")
					else: #W2NF
						if index in sterile_ids:
							b7 += 1
						elif index in fertile_ids:
							b8 += 1
						else:
							raise Exception("Error in calculations; Sterile/Fertile indices are incorrect")
			bd_1.append(b1)
			bd_2.append(b2)
			bd_3.append(b3)
			bd_4.append(b4)
			bd_5.append(b5)
			bd_6.append(b6)
			bd_7.append(b7)
			bd_8.append(b8)

if args.uset == 'M':
	if args.focal == 0:
		print("Running unidirectional scan, 0 Mito focal genotype")
		BD_scan(mito_win, len_X, 0)
		BD_scan(mito_win, len_2, 0)
		BD_scan(mito_win, len_3, 0)
		out_name = args.o + '_mito_FR_scan_bd_cells.csv'
	else:
		print("Running unidirectional scan, 2 Mito focal genotype")
		BD_scan(mito_win, len_X, 2)
		BD_scan(mito_win, len_2, 2)
		BD_scan(mito_win, len_3, 2)
		out_name = args.o + '_mito_ZI_scan_bd_cells.csv'
if args.uset == 'Y':
	if args.focal == 0:
		print("Running unidirectional scan, 0 Y focal genotype")
		BD_scan(y_win, len_X, 0)
		BD_scan(y_win, len_2, 0)
		BD_scan(y_win, len_3, 0)
		out_name = args.o + '_Y_FR_scan_bd_cells.csv'
	else:
		print("Running unidirectional scan, 2 Y focal genotype")
		BD_scan(y_win, len_X, 2)
		BD_scan(y_win, len_2, 2)
		BD_scan(y_win, len_3, 2)
		out_name = args.o + '_Y_ZI_scan_bd_cells.csv'

##output
cell_df = pd.DataFrame([bd_1, bd_2, bd_3, bd_4, bd_5, bd_6, bd_7, bd_8])
cell_df = cell_df.transpose()
cell_df.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
cell_df.to_csv(out_name, header=True, index=False)
