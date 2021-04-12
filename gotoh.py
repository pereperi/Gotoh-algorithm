#!/usr/bin/env python

import sys
import re
from Bio import SeqIO

def read_matrix(filename):	#function to read substitution matrix
	matrix = {}	#initialize matrix
	f = open(filename)	#open file
	chars = f.readline().rstrip().split('\t')	#read first line (separated by tabs)
	aa_alphabet = []	#initialize a list that will contain all aa
	for x in chars[1:]:
		aa_alphabet.append(x)	#add the amino acids to aa_alphabet
	lines = f.readlines()	#read the rest of the lines
	row = 0	#to keep track of which row we are reading
	for line in lines:	#for every line (row)
		matrix[aa_alphabet[row]] = {}	#initialize a dictionary inside the main dictionary for the aa we are reading (one different every row)
		values = line.rstrip().split('\t')	#get the line splitted by tabs and without the endlines
		for col in range(len(aa_alphabet)):	#every element in values contains a value for one aa, with the same order as aa_alphabet
			matrix[aa_alphabet[row]][aa_alphabet[col]] = int(values[col+1])	#add that value to the dictionary of the aa of the current row
		row += 1	#when done with the whole line, increase row by one and repeat process
	return matrix	#return the matrix (as a dictionary of dictionaries), and the list of all amino acids

		
record = list(SeqIO.parse(sys.argv[1], "fasta"))

seqI = re.sub(r"[\n-]", "", "".join(record[0].seq))
seqJ = re.sub(r"[\n-]", "", "".join(record[1].seq))

matrix = read_matrix(sys.argv[2])
lenI = len(seqI)
lenJ = len(seqJ)

gop = -10
gep = -0.5

m = [[0 for x in range(lenJ + 1)] for y in range(lenI + 1)]
ix = [[0 for x in range(lenJ + 1)] for y in range(lenI + 1)]
iy = [[0 for x in range(lenJ + 1)] for y in range(lenI + 1)]
tbm = [[0 for x in range(lenJ + 1)] for y in range(lenI + 1)]
tbx = [[0 for x in range(lenJ + 1)] for y in range(lenI + 1)]
tby = [[0 for x in range(lenJ + 1)] for y in range(lenI + 1)]

for i in range(1, lenI + 1):
	iy[i][0] = gop + (i-1) * gep
	ix[i][0] = float('-inf')
	m[i][0] = float('-inf')
	tbm[i][0] = 1
	tbx[i][0] = 1
	tby[i][0] = 1
	

for j in range(1, lenJ + 1):
	ix[0][j] = gop + (j-1) * gep
	iy[0][j] = float('-inf')
	m[0][j] = float('-inf')
	tbm[0][j] = -1
	tbx[0][j] = -1
	tby[0][j] = -1

for i in range(1, lenI + 1):
	for j in range(1, lenJ + 1):
		s = matrix[seqI[i-1]][seqJ[j-1]]
		#M
		sub = m[i-1][j-1] + s
		x = ix[i-1][j-1] + s
		y = iy[i-1][j-1] + s
		if sub >= x and sub >= y:
			m[i][j] = sub
			tbm[i][j] = 0
		elif x > y:
			m[i][j] = x
			tbm[i][j] = -1
		else:
			m[i][j] = y
			tbm[i][j] = 1
		#Ix	
		sub = m[i][j-1] + gop
		x = ix[i][j-1] + gep
		if sub >= x:
			ix[i][j] = sub
			tbx[i][j] = 0
		else:
			ix[i][j] = x
			tbx[i][j] = -1
		#Iy
		sub = m[i-1][j] + gop
		y = iy[i-1][j] + gep
		if sub >= y:
			iy[i][j] = sub
			tby[i][j] = 0
		else:
			iy[i][j] = y
			tby[i][j] = 1


print('Optimal score: ', max(m[i][j], ix[i][j], iy[i][j]))
i = lenI
j = lenJ
alnI = []
alnJ = []

if m[i][j] >= ix[i][j] and m[i][j] >= iy[i][j]:
	state = 0
elif ix[i][j] > iy[i][j]:
	state = -1
else:
	state = 1
	
while i != 0 or j != 0:
	if state == 0:
		state = tbm[i][j]
		i += -1
		j += -1
		alnI.append(seqI[i])
		alnJ.append(seqJ[j])
	elif state == -1:
		state = tbx[i][j]
		j += -1
		alnI.append("-")
		alnJ.append(seqJ[j])
	else:
		state = tby[i][j]
		i += -1
		alnI.append(seqI[i])
		alnJ.append("-")

seqI_aln = "".join(reversed(alnI))
seqJ_aln = "".join(reversed(alnJ))

symbols = ''
identity = 0
similarity = 0
total = 0
posI = 0
posJ = 0
pos = {}
for x in range(len(seqI_aln)):
	if seqI_aln[x] == seqJ_aln[x]:
		symbols += '*'
		identity += 1
		similarity += 1
		total += 1
		posI += 1
		posJ += 1
	elif seqI_aln[x] != '-' and seqJ_aln[x] != '-':	
		total += 1
		posI += 1
		posJ += 1
		if matrix[seqI_aln[x]][seqJ_aln[x]] in (0,1):
			symbols += '.'
			similarity += 1
		elif matrix[seqI_aln[x]][seqJ_aln[x]] >= 2:
			similarity += 1
			symbols += ':'
		else:
			symbols += ' '
	else:
		symbols += ' '
		if seqI_aln[x] != '-':
			posI += 1
		elif seqJ_aln[x] != '-':
			posJ += 1
	if (x+1)%100 == 0:
		pos[x+1] = [posI,posJ]	#for each slice of 100 that will be printed, at which position of each sequence am I?

print('Identity:', str((identity/total)*100)+'%', f"({identity}/{total})")	
print('Similarity:', str((similarity/total)*100)+'%', f"({similarity}/{total})")
print()

l = len(seqI_aln)
part = l//100
ind = 0
for x in range(part):
	print(seqI_aln[ind:ind+100], pos[ind+100][0])
	print(symbols[ind:ind+100])
	print(seqJ_aln[ind:ind+100], pos[ind+100][1])
	print()
	ind += 100

print(seqI_aln[ind:], lenI)
print(symbols[ind:])
print(seqJ_aln[ind:], lenJ)
	

