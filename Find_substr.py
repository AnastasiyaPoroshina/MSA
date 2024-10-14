import time
import sys
import random
import argparse
import os
import itertools
from itertools import chain
from functools import reduce
from threading import Timer
import functools
print = functools.partial(print, flush=True)
import re



#Parser
def createParser ():
	parser = argparse.ArgumentParser()
	parser.add_argument ('-n', '--name', nargs='?', default='GCF_000004095.1_Hydra_RP_1.0_genomic.fna')
	parser.add_argument ('-o', '--out_name', nargs='?', default='')
	return parser

#Получаем геномные строки из строк
def get_strain_sequence(clist):
	i = 0
	strain_list = []
	global strain_names
	while i < len(clist):
		if clist[i][0] == '>':
			strain_names.append(clist[i])
			strain_clist = []
			end_flag = True
			try:
				while clist[i+1][0] != '>':
					strain_clist.append(clist[i+1]) 
					i+=1
			except:
				end_flag = False
			if end_flag == True:
				strain_list.append(strain_clist)
				#print("get sequence", len(strain_list) ,", number of strings is", len(strain_clist))
			i+=1
		else:
			i+=1
	return strain_list


#Объединяем геномные строки в большие строки
def union_strain_sequence(clist):
	i = 0
	strain_string_list = []
	for y in clist:
		strain_string = ""
		for x in y:
			strain_string += x
		strain_string_list.append(strain_string)
	return strain_string_list

def find_subsequence(str1,str2,f):
	list_of_occurences = [m.start() for m in re.finditer(str1, str2)]
	print(str1, file = f)
	print(list_of_occurences, file = f)


def capitalize_list_of_strings(clist):
	new_clist = []
	for string in clist:
		new_string = string.upper()
		new_clist.append(new_string)
	return new_clist


#================================================= MAIN ===============================================================


parser = createParser()
namespace = parser.parse_args (sys.argv[1:])

#Имя входной задачи
input_name = namespace.name
input_name_path = "C:\\Users\\Vsk\\Desktop\\" + input_name
out_file = "C:\\Users\\Vsk\\Desktop\\out_" + input_name



try:
	input_file = open(input_name_path)
except:
	print("file is not opened")
	exit()


input_file_list = input_file.readlines()
if input_file_list[-1][0] != '>':
	input_file_list.append('>')


strain_names = []
strain_lists_divided = get_strain_sequence(input_file_list)
strain_names.pop()
strain_strings_list = union_strain_sequence(strain_lists_divided)



nucles = ['A','T','G','C']

subsequences1 = itertools.product(nucles, repeat=1)
string_subsequence1 = []
for subsequence in subsequences1:
	new_subsequence = ""
	for nucle in subsequence:
		new_subsequence += nucle
	string_subsequence1.append(new_subsequence*10)

subsequences2 = itertools.product(nucles, repeat=2)
string_subsequence2 = []
for subsequence in subsequences2:
	new_subsequence = ""
	for nucle in subsequence:
		new_subsequence += nucle
	string_subsequence2.append(new_subsequence*10)

subsequences3 = itertools.product(nucles, repeat=3)
string_subsequence3 = []
for subsequence in subsequences3:
	new_subsequence = ""
	for nucle in subsequence:
		new_subsequence += nucle
	string_subsequence3.append(new_subsequence*10)


subsequences = []
subsequences.extend(string_subsequence1)
subsequences.extend(string_subsequence2)
subsequences.extend(string_subsequence3)

#print(*subsequences,sep = '\n')

strain_strings_list = capitalize_list_of_strings(strain_strings_list)



with open(out_file, 'w') as f:
	print("====================================================================================================================================", file=f)
	print("Input file name", input_name, file=f)
	print("Absolute input file name", input_name_path, file=f)
	print("Number of strings in input file", len(input_file_list), file=f)
	print("Number of strains", len(strain_strings_list), file=f)
	print("Number of strain names", len(strain_names), file=f)
	print("====================================================================================================================================", file=f)
	print("====================================================================================================================================", file=f)
	for i in range(1):
	#for i in range(len(strain_names)):
		print(strain_names[i], file=f)
		print("Len of strain sequence is", len(strain_strings_list[i]), file=f)
		for subsequence in subsequences:
			find_subsequence(subsequence,strain_strings_list[i], f)
		print("", file = f)
		print("====================================================================================================================================", file=f)
		print("", file = f)

