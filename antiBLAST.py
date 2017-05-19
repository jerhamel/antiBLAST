__author__ = 'Jeremie Hamel'

import os, subprocess
#################################################
# antiBLAST creates fasta files from antiSMASH  #
# cluster results, runs BLAST queries and sorts #
# the results into .csv files.					#
# It is composed of 3 main parts: Creator,		#
# Executor and Analyzor							#
#################################################

#################################################
#ANALYZOR
#################################################
def read_result(infile):
	dict_lines = {}             #Initiate linecount dictionary
	line_num = 0
	output_type = "-"
	for line in infile:
		line = line.split("\n")
		dict_lines[line_num] = line  # Build line number dictionnary
		line_num += 1

		if "Query=" in str(line):
			full_line = str(line).split(" ")        #Isolate ID components
			full_id = (full_line[1])[:-2]
			full_id_list = full_id.split("_")
			bacteria = full_id_list[0]+"_"+full_id_list[1]
			scaffold = full_id_list[2]
			pathway = full_id_list[3]
			orf = full_id_list[4]+"_"+full_id_list[5]

	if infile.endswith("nr.out"):
		output_type = "nr"
	if infile.endswith("pdb.out"):
		output_type = "pdb"

	if output_type == "nr":
		top_hit_accession, description, e_value, identity_percent, positive_percent, gap_percent = process_nr(dict_lines)

	if output_type == "pdb":
		top_hit_accession, description, e_value, identity_percent, positive_percent, gap_percent = process_pdb(dict_lines)

	return bacteria, scaffold, pathway, orf, top_hit_accession, description, e_value, identity_percent, positive_percent, gap_percent, output_type


def process_nr(dict_lines):
	for line_num in dict_lines:
		if "Sequences producing" in str(dict_lines[line_num]):      #Find top hit accession number
			top_hit_line = dict_lines[line_num+2]
			top_hit_accession = (str(top_hit_line).split(' '))[0]
			top_hit_accession = top_hit_accession[2:]

			for line_num in dict_lines:
				if str(">"+top_hit_accession) in str(dict_lines[line_num]):     #Find top hit full result
					description = str(dict_lines[line_num]).split(' ')
					description = description[1:]
					description = ' '.join(description)
					description = description[:-6]
					check_line = line_num+1
					while "Score =" not in str(dict_lines[check_line]):         #Find E-value
						check_line += 1
					e_value = (str(dict_lines[check_line]).split(","))[1]
					e_value = e_value.split(" ")[4]

					following_line = str(dict_lines[check_line+1]).split(",")   #Find identities, positives and gaps
					following_line = str(following_line).split(" ")

					identity_percent = (following_line[4])[1:-3]
					positive_percent = (following_line[9])[1:-3]
					gap_percent = (following_line[14])[1:-4]


	return top_hit_accession, description, e_value, identity_percent, positive_percent, gap_percent

def process_pdb(dict_lines):                    #Almost identical to process_nr so far, in case I wanna change stuff
	for line_num in dict_lines:
		if "Sequences producing" in str(dict_lines[line_num]):      #Find top hit accession number
			top_hit_line = dict_lines[line_num+2]
			top_hit_accession = (str(top_hit_line).split(' '))[0]
			top_hit_accession = top_hit_accession[2:]

			for line_num in dict_lines:
				if str(">"+top_hit_accession) in str(dict_lines[line_num]):     #Find top hit full result
					description = str(dict_lines[line_num]).split(' ')
					description = description[3:]                           #Only difference w/ other function
					description = ' '.join(description)
					description = description[:-6]
					check_line = line_num+1
					while "Score =" not in str(dict_lines[check_line]):         #Find E-value
						check_line += 1
					e_value = (str(dict_lines[check_line]).split(","))[1]
					e_value = e_value.split(" ")[4]

					following_line = str(dict_lines[check_line+1]).split(",")   #Find identities, positives and gaps
					following_line = str(following_line).split(" ")

					identity_percent = (following_line[4])[1:-3]
					positive_percent = (following_line[9])[1:-3]
					gap_percent = (following_line[14])[1:-4]


	return top_hit_accession, description, e_value, identity_percent, positive_percent, gap_percent

def write_result(result):
	result_split = result.split(",")
	file_name = result_split[0]+".cvs"
	header = "Bacteria,Scaffold,Pathway,ORF,Database,Top Hit Acc. Num.,Description,E-Value,Identity %,Positive %,Gap %"
	if file_name not in os.listdir():           #Find if file already present
		print("Creating", file_name)
		file = open(file_name, "w")             #Create file
		file.write(header+"\n"+result)          #Write
		file.close()                            #Close
	else:
		print("Appending result to", file_name)
		file =open(file_name, "a")
		file.write("\n"+result)
		file.close()


def analyzor(path):
	path_out = path+"out/"
	results_list = []
	for file in os.listdir(path):
		print(file)
		file_path = path_out + file
		infile = open(file_path, 'r')
		bacteria, scaffold, orf, pathway, top_hit_accession, description, e_value, identity_percent, positive_percent, gap_percent, output_type = read_result(
			infile)
		infile.close()
		results = [bacteria, scaffold, orf, output_type, top_hit_accession, description, e_value, identity_percent,
				   positive_percent, gap_percent]
		results = ",".join(results)  # Join for listing
		results_list.append(results)

	results_list_sorted = sorted(results_list)  # Sort for easier file writing and appending
	for result_sorted in results_list_sorted:
		write_result(result_sorted)

	print("RESULT SORTING COMPLETED")
	print("Job complete. Thank you for choosing antiBLAST 1.0 suite")
	print("The CPU will now autodestruct")
	print("(Just kidding)")

#################################################
#CREATOR
#################################################

def find_subdirectories(path):
	dir_list2=[]             # Initialize list of directories
	dir_number = 0
	sub = os.listdir(path)  # Lists all files and directories in path

	for p in sub:
		dir_number += 1
		pDir = os.path.join(path, p)  # Joins path and files/directory name for correct file path
		if os.path.isdir(pDir):  # Enquires if file or directory
			dir_list2.append(pDir)
			#print(pDir)

	return dir_list2, dir_number

def find_cluster_files(sub):
	cluster_file_list = []
	cluster_number = 0
	for file in os.listdir(sub):        #Search each subdirectory
		if file.endswith('geneclusters.txt'):
			cluster_file = os.path.join(sub, file)
			cluster_file_list.append(cluster_file)
			cluster_number += 1

	return cluster_file_list, cluster_number


def read_sequence(orf, sequence_file):
	sequence = "-"
	seqfile = open(sequence_file, "r")
	line_number = 0                 #Line count
	dict_lines = {}
	for line1 in seqfile:

		line = line1.split('\n')
		dict_lines[line_number] = line      #Build line number dictionnary
		line_number += 1

	for line_num in dict_lines:
		if "/locus_tag=" in str(dict_lines[line_num]):      #Find ORF
			if orf in str(dict_lines[line_num]):
				target_line = line_num + 2          #Find translation line

				if "/translation=" in str(dict_lines[target_line]):
					full_line_list  = str(dict_lines[target_line]).split('"')
					local_line_list = str(full_line_list[1])
					seq_line1 = local_line_list.split("'")
					sequence = seq_line1[0]
					target_line += 1
					while "." not in str(dict_lines[target_line]):     #Filter line
						seq_line2 = str(dict_lines[target_line]).split()
						seq_line2 = str(seq_line2[1])
						seq_line2 = seq_line2.split("'")
						seq_line2 = seq_line2[0]
						if seq_line2.endswith('"'):
							seq_line2 = seq_line2[:-1]      #Remove last caracter
						sequence = sequence + seq_line2
						target_line += 1


	seqfile.close()

	return sequence

def read_cluster(file, ID, path):
	infile = open(file, 'r')            #Open geneclusters.txt
	scaffold_prog = 0
	for line in infile:
		line2 = line.strip('\n')         #Separate by returns
		scaffold_prog += 1
		scaffold_elements = line2.split('\t')        #Separate by tabs
		orfs_list = scaffold_elements[3].split(';')     #Separate ORFs into a list
		description = ID+"_"+scaffold_elements[0]+"_"+scaffold_elements[2]
		sequence_file = path+"/scaffold1.1.final.gbk"       #Sequence file name and path
		print(len(orfs_list), "ORFs found for", ID)
		orf_prog = 0
		for orf in orfs_list:
			orf_prog += 1               #Progression
			sequence = read_sequence(orf, sequence_file)        #Read and return sequence
			header = ">"+description+"_"+orf
			filename = "fasta/"+ID+"/"+ID+"_"+orf+".fasta"
			#print(orf_prog, "of", len(orfs_list), " ORFs in scaffold", scaffold_prog, "of", len(line))

			if not sequence == "-":
				if not os.path.exists("fasta/"+ID):
					os.makedirs("fasta/"+ID)
				print("Writing", filename)
				fasta_file = open(filename, "w+")           #Create fasta file
				fasta_file.write(header+"\n"+sequence)      #Write
				fasta_file.close()                          #Close

	infile.close()

def creator(path):

	run_name_index = 0						#Fix name selection when path changes
	for character in str(path):
		if character == "/":
			run_name_index += 1
	dir_list, dir_number = find_subdirectories(path)
	for sub_path in dir_list:
		cluster_list_path, cluster_number = find_cluster_files(sub_path)

		for cluster_file in cluster_list_path:
			full_name = sub_path.split('/')  # Find bacteria ID
			run_name = full_name[run_name_index]
			split_run_name = run_name.split('.')
			bacteria_ID = split_run_name[0]
			read_cluster(cluster_file, bacteria_ID, sub_path)

	print("FASTA FILE CREATION COMPLETED")

###################################################
#EXECUTOR
###################################################

def executor(blastp_path):
	sub = os.listdir("fasta/")          #List directories in fasta/
	#dir_number = len(sub)
	subprocess.call("mkdir out", shell=True)

	for subdir in sub:                  #Find fasta files in each Run directory
		fasta_list = []
		sub_fasta = os.listdir("fasta/"+subdir)

		for file in sub_fasta:          #To make sure not to create a list of lists
			file_path = os.path.abspath("fasta/"+subdir+"/"+file)            #Get absolute path
			fasta_list.append(file_path)

		for file1 in fasta_list:
			filename_core = (file1.split("/"))[-1]
			filename_core = (filename_core.split("."))[0]

			outname_nr = "out/"+filename_core+"_nr.out"
			outname_pdb = "out/"+filename_core+"_pdb.out"

			blast_command_nr = blastp_path + " -query " + file1 + " -out " + outname_nr + " -db nr -remote"
			blast_command_pdb = blastp_path + " -query " + file1 + " -out " + outname_pdb + " -db pdb -remote"

			print("Running blastp NR database query on", filename_core)         #Run NR query
			subprocess.call(blast_command_nr, shell=True)

			print("Running blastp PDB database query on", filename_core)        #Run PDB query
			subprocess.call(blast_command_pdb, shell=True)

	print("BLAST QUERIES COMPLETED")

if __name__ == '__main__':
	working_dir_path = (os.getcwd())+"/"				#Get working directory
	blastp_path = "/home/jhamel/local/ncbi-blast-2.6.0+/bin/blastp"

	print("Creating fasta files")
	creator(working_dir_path)
	print("Conducting BLAST queries")
	executor(blastp_path)
	print("Sorting results")
	analyzor(working_dir_path)