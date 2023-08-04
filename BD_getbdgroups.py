#### Author: Matthew Lollar
#### QTL scan framework using the Breslow-Day test statisitic
#### Last update: July 5th, 2023
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
parser.add_argument('--u', help='Run unidirectional scans (optional)', action='store_true', required=False)
parser.add_argument('--uset', help='Run a subset of unidirectional scans (0=FRmito/ZIy, 1=both, 2=ZImito/FRy)', type=int, choices=range(0, 3), required=False)
args = parser.parse_args()
if args.u and (args.uset is None):
	parser.error("Flag '--u' requires additional flag '--uset.'")

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

### Unidirectional scan
if args.u == True:
	mito_win = [2579] #Adjust as needed
	y_win = [2580] #Adjust as needed
	### Initialize Breslow-Day cell count lists
	#forward scans
	bd_1f_m = []
	bd_2f_m = []
	bd_3f_m = []
	bd_4f_m = []
	bd_5f_m = []
	bd_6f_m = []
	bd_7f_m = []
	bd_8f_m = []
	bd_1f_y = []
	bd_2f_y = []
	bd_3f_y = []
	bd_4f_y = []
	bd_5f_y = []
	bd_6f_y = []
	bd_7f_y = []
	bd_8f_y = []
	#reverse scans
	bd_1r_m = []
	bd_2r_m = []
	bd_3r_m = []
	bd_4r_m = []
	bd_5r_m = []
	bd_6r_m = []
	bd_7r_m = []
	bd_8r_m = []
	bd_1r_y = []
	bd_2r_y = []
	bd_3r_y = []
	bd_4r_y = []
	bd_5r_y = []
	bd_6r_y = []
	bd_7r_y = []
	bd_8r_y = []

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
					else:
						bd_8r2.append(b8)
						bd_1r2.append(b1)
						bd_2r2.append(b2)
						bd_3r2.append(b3)
						bd_4r2.append(b4)
						bd_5r2.append(b5)
						bd_6r2.append(b6)
						bd_7r2.append(b7)
						bd_8r2.append(b8)
			elif scantype == 'mito':
				if focal == 0: #add to either forward or reverse lists
					bd_1f_m.append(b1)
					bd_2f_m.append(b2)
					bd_3f_m.append(b3)
					bd_4f_m.append(b4)
					bd_5f_m.append(b5)
					bd_6f_m.append(b6)
					bd_7f_m.append(b7)
					bd_8f_m.append(b8)
				else:
					bd_1r_m.append(b1)
					bd_2r_m.append(b2)
					bd_3r_m.append(b3)
					bd_4r_m.append(b4)
					bd_5r_m.append(b5)
					bd_6r_m.append(b6)
					bd_7r_m.append(b7)
					bd_8r_m.append(b8)
			elif scantype == 'y':
				if focal == 0: #add to either forward or reverse lists
					bd_1f_y.append(b1)
					bd_2f_y.append(b2)
					bd_3f_y.append(b3)
					bd_4f_y.append(b4)
					bd_5f_y.append(b5)
					bd_6f_y.append(b6)
					bd_7f_y.append(b7)
					bd_8f_y.append(b8)
				else:
					bd_1r_y.append(b1)
					bd_2r_y.append(b2)
					bd_3r_y.append(b3)
					bd_4r_y.append(b4)
					bd_5r_y.append(b5)
					bd_6r_y.append(b6)
					bd_7r_y.append(b7)
					bd_8r_y.append(b8)
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
cell_df_f0 = cell_df_f.transpose()
cell_df_r0 = cell_df_r.transpose()
cell_df_f2 = cell_df_f.transpose()
cell_df_r2 = cell_df_r.transpose()
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

### Run BD scans
### Warning: No sanity check for equal bd groups in uni subset
if args.u == True:
	if args.uset == 1:
		#forward mito
		print("getting focal 0 mito")
		BD_scan(mito_win, len_X, 0, 'mito')
		BD_scan(mito_win, len_2, 0, 'mito')
		BD_scan(mito_win, len_3, 0, 'mito')

		#reverse mito
		print("getting focal 2 mito")
		BD_scan(mito_win, len_X, 2, 'mito')
		BD_scan(mito_win, len_2, 2, 'mito')
		BD_scan(mito_win, len_3, 2, 'mito')

		#forward y
		print("getting focal 0 Y")
		BD_scan(y_win, len_X, 0, 'y')
		BD_scan(y_win, len_2, 0, 'y')
		BD_scan(y_win, len_3, 0, 'y')

		#reverse y
		print("getting focal 2 Y")
		BD_scan(y_win, len_X, 2, 'y')
		BD_scan(y_win, len_2, 2, 'y')
		BD_scan(y_win, len_3, 2, 'y')

		##output
		cell_df_mito_f = pd.DataFrame([bd_1f_m, bd_2f_m, bd_3f_m, bd_4f_m, bd_5f_m, bd_6f_m, bd_7f_m, bd_8f_m])
		cell_df_mito_r = pd.DataFrame([bd_1r_m, bd_2r_m, bd_3r_m, bd_4r_m, bd_5r_m, bd_6r_m, bd_7r_m, bd_8r_m])
		cell_df_y_f = pd.DataFrame([bd_1f_y, bd_2f_y, bd_3f_y, bd_4f_y, bd_5f_y, bd_6f_y, bd_7f_y, bd_8f_y])
		cell_df_y_r = pd.DataFrame([bd_1r_y, bd_2r_y, bd_3r_y, bd_4r_y, bd_5r_y, bd_6r_y, bd_7r_y, bd_8r_y])
		cell_df_mito_f = cell_df_mito_f.transpose()
		cell_df_mito_r = cell_df_mito_r.transpose()
		cell_df_y_f = cell_df_y_f.transpose()
		cell_df_y_r = cell_df_y_r.transpose()
		cell_df_mito_f.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
		cell_df_mito_r.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
		cell_df_y_f.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
		cell_df_y_r.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
		out_name_mito_f = args.o + '_mito_uni_FR_scan_bd_cells.csv'
		out_name_mito_r = args.o + '_mito_uni_ZI_scan_bd_cells.csv'
		out_name_y_f = args.o + '_y_uni_FR_scan_bd_cells.csv'
		out_name_y_r = args.o + '_y_uni_ZI_scan_bd_cells.csv'
		cell_df_mito_f.to_csv(out_name_mito_f, header=True, index=False)
		cell_df_mito_r.to_csv(out_name_mito_r, header=True, index=False)
		cell_df_y_f.to_csv(out_name_y_f, header=True, index=False)
		cell_df_y_r.to_csv(out_name_y_r, header=True, index=False)

	if args.uset == 0: #FRmito, ZIy
		#forward mito
		print("getting 0 mito, 2 Y scans")
		BD_scan(mito_win, len_X, 0, 'mito')
		BD_scan(mito_win, len_2, 0, 'mito')
		BD_scan(mito_win, len_3, 0, 'mito')

		#reverse y
		BD_scan(y_win, len_X, 2, 'y')
		BD_scan(y_win, len_2, 2, 'y')
		BD_scan(y_win, len_3, 2, 'y')

		cell_df_mito_f = pd.DataFrame([bd_1f_m, bd_2f_m, bd_3f_m, bd_4f_m, bd_5f_m, bd_6f_m, bd_7f_m, bd_8f_m])
		cell_df_y_r = pd.DataFrame([bd_1r_y, bd_2r_y, bd_3r_y, bd_4r_y, bd_5r_y, bd_6r_y, bd_7r_y, bd_8r_y])
		cell_df_mito_f = cell_df_mito_f.transpose()
		cell_df_y_r = cell_df_y_r.transpose()
		cell_df_mito_f.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
		cell_df_y_r.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
		out_name_mito_f = args.o + '_mito_uni_FR_scan_bd_cells.csv'
		out_name_y_r = args.o + '_y_uni_ZI_scan_bd_cells.csv'
		cell_df_mito_f.to_csv(out_name_mito_f, header=True, index=False)
		cell_df_y_r.to_csv(out_name_y_r, header=True, index=False)

	if args.uset == 2: #ZImito, FRy
		#reverse mito
		print("getting 2 mito, 0 y scans")
		BD_scan(mito_win, len_X, 2, 'mito')
		BD_scan(mito_win, len_2, 2, 'mito')
		BD_scan(mito_win, len_3, 2, 'mito')

		#forward y
		BD_scan(y_win, len_X, 0, 'y')
		BD_scan(y_win, len_2, 0, 'y')
		BD_scan(y_win, len_3, 0, 'y')

		cell_df_mito_r = pd.DataFrame([bd_1r_m, bd_2r_m, bd_3r_m, bd_4r_m, bd_5r_m, bd_6r_m, bd_7r_m, bd_8r_m])
		cell_df_y_f = pd.DataFrame([bd_1f_y, bd_2f_y, bd_3f_y, bd_4f_y, bd_5f_y, bd_6f_y, bd_7f_y, bd_8f_y])
		cell_df_mito_r = cell_df_mito_r.transpose()
		cell_df_y_f = cell_df_y_f.transpose()
		cell_df_mito_r.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
		cell_df_y_f.columns = ['bd1','bd2','bd3','bd4','bd5','bd6','bd7','bd8']
		out_name_mito_r = args.o + '_mito_uni_ZI_scan_bd_cells.csv'
		out_name_y_f = args.o + '_y_uni_FR_scan_bd_cells.csv'
		cell_df_mito_r.to_csv(out_name_mito_r, header=True, index=False)
		cell_df_y_f.to_csv(out_name_y_f, header=True, index=False)
