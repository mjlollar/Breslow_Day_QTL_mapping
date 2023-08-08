#### Author: Matthew Lollar
#### QTL scan framework using the Breslow-Day test statisitic
#### Last update: Aug 8th, 2023
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
#forward scans (0focal)
bd_1f0 = []
bd_2f0 = []
bd_3f0 = []
bd_4f0 = []
bd_5f0 = []
bd_6f0 = []
bd_7f0 = []
bd_8f0 = []
#reverse scans (2focal)
bd_1r2 = []
bd_2r2 = []
bd_3r2 = []
bd_4r2 = []
bd_5r2 = []
bd_6r2 = []
bd_7r2 = []
bd_8r2 = []
#forward scans (2focal)
bd_1f2 = []
bd_2f2 = []
bd_3f2 = []
bd_4f2 = []
bd_5f2 = []
bd_6f2 = []
bd_7f2 = []
bd_8f2 = []
#reverse scans (0focal)
bd_1r0 = []
bd_2r0 = []
bd_3r0 = []
bd_4r0 = []
bd_5r0 = []
bd_6r0 = []
bd_7r0 = []
bd_8r0 = []

### Scan function
def BD_scan(chr_1, chr_2, focal, scantype, direction):
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
			if scantype == 'bi':
				if focal == 0:
					if direction == 'forward':
						bd_1f0.append(b1)
						bd_2f0.append(b2)
						bd_3f0.append(b3)
						bd_4f0.append(b4)
						bd_5f0.append(b5)
						bd_6f0.append(b6)
						bd_7f0.append(b7)
						bd_8f0.append(b8)
					else:
						bd_1r0.append(b1)
						bd_2r0.append(b2)
						bd_3r0.append(b3)
						bd_4r0.append(b4)
						bd_5r0.append(b5)
						bd_6r0.append(b6)
						bd_7r0.append(b7)
						bd_8r0.append(b8)
				else:
					if direction == 'forward':
						bd_1f2.append(b1)
						bd_2f2.append(b2)
						bd_3f2.append(b3)
						bd_4f2.append(b4)
						bd_5f2.append(b5)
						bd_6f2.append(b6)
						bd_7f2.append(b7)
						bd_8f2.append(b8)
					else:
						bd_1r2.append(b1)
						bd_2r2.append(b2)
						bd_3r2.append(b3)
						bd_4r2.append(b4)
						bd_5r2.append(b5)
						bd_6r2.append(b6)
						bd_7r2.append(b7)
						bd_8r2.append(b8)
			else:
				print('sanity check, problem with "BD_scan" function parsing')

### Run BD scans, forward
print("getting X by 2, 0 focal")
BD_scan(len_X, len_2, 0, 'bi', 'forward')
print("getting X by 3, 0 focal")
BD_scan(len_X, len_3, 0, 'bi', 'forward')
print("getting 2 by 3, 0 focal")
BD_scan(len_2, len_3, 0, 'bi', 'forward')
print("getting X by 2, 2 focal")
BD_scan(len_X, len_2, 2, 'bi', 'forward')
print("getting X by 3, 2 focal")
BD_scan(len_X, len_3, 2, 'bi', 'forward')
print("getting 2 by 3, 2 focal")
BD_scan(len_2, len_3, 2, 'bi','forward')

### Run BD scans, reverse
print("getting 2 by X, 0 focal")
BD_scan(len_2, len_X, 0, 'bi', 'r')
print("getting X by 3, 0 focal")
BD_scan(len_3, len_X, 0, 'bi', 'r')
print("getting 3 by 2, 0 focal")
BD_scan(len_3, len_2, 0, 'bi', 'r')
print("getting 2 by X, 2 focal")
BD_scan(len_2, len_X, 2, 'bi', 'r')
print("getting 3 by X, 2 focal")
BD_scan(len_3, len_X, 2, 'bi', 'r')
print("getting 3 by 2, 2 focal")
BD_scan(len_3, len_2, 2, 'bi', 'r')


assert len(bd_1f0) == len(bd_2f0) == len(bd_3f0) == len(bd_4f0) == len(bd_5f0) == len(bd_6f0) == len(bd_7f0) == len(bd_8f0) # Sanity
assert len(bd_1r0) == len(bd_2r0) == len(bd_3r0) == len(bd_4r0) == len(bd_5r0) == len(bd_6r0) == len(bd_7r0) == len(bd_8r0) # Sanity
assert len(bd_1f2) == len(bd_2f2) == len(bd_3f2) == len(bd_4f2) == len(bd_5f2) == len(bd_6f2) == len(bd_7f2) == len(bd_8f2) # Sanity
assert len(bd_1r2) == len(bd_2r2) == len(bd_3r2) == len(bd_4r2) == len(bd_5r2) == len(bd_6r2) == len(bd_7r2) == len(bd_8r2) # Sanity

#Output
cell_df_f0 = pd.DataFrame([bd_1f0, bd_2f0, bd_3f0, bd_4f0, bd_5f0, bd_6f0, bd_7f0, bd_8f0])
cell_df_r0 = pd.DataFrame([bd_1r0, bd_2r0, bd_3r0, bd_4r0, bd_5r0, bd_6r0, bd_7r0, bd_8r0])
cell_df_f2 = pd.DataFrame([bd_1f2, bd_2f2, bd_3f2, bd_4f2, bd_5f2, bd_6f2, bd_7f2, bd_8f2])
cell_df_r2 = pd.DataFrame([bd_1r2, bd_2r2, bd_3r2, bd_4r2, bd_5r2, bd_6r2, bd_7r2, bd_8r2])
cell_df_f0 = cell_df_f0.transpose()
cell_df_r0 = cell_df_r0.transpose()
cell_df_f2 = cell_df_f2.transpose()
cell_df_r2 = cell_df_r2.transpose()
cell_df_f0.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
cell_df_r0.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
cell_df_f2.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
cell_df_r2.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
out_name_f0 = args.o + '_0focal_forward_scan_bd_cells.csv'
out_name_r0 = args.o + '_0focal_reverse_scan_bd_cells.csv'
out_name_f2 = args.o + '_2focal_forward_scan_bd_cells.csv'
out_name_r2 = args.o + '_2_focal_reverse_scan_bd_cells.csv'
cell_df_f0.to_csv(out_name_f0, header=True, index=False)
cell_df_r0.to_csv(out_name_r0, header=True, index=False)
cell_df_f2.to_csv(out_name_f2, header=True, index=False)
cell_df_r2.to_csv(out_name_r2, header=True, index=False)
