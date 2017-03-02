###################################################################

# Developer: Cosi Jackson
# Date: February 10, 2017
# Professor: David Knox
# Class: Algorithms for Molecular Biology 
# Assignment 2

# The purpose of this assigment is to determine the top 5 most common k-mers that 
# appear in the file and the top 5 most common k-mers that appear in the most sequences.
# The time complexity  and memory limitations of this assignement will also be determined.


# referencing: http://stackoverflow.com/questions/11011756/is-there-any-pythonic-way-to-combine-two-dicts-adding-values-for-keys-that-appe
#			   http://stackoverflow.com/questions/29240807/python-collections-counter-most-common-complexity
#			   https://docs.python.org/2/library/collections.html


###################################################################

# Big O Notation - Time Complexity of the program

# When the parse_fasta function is called it runs in: O(n) - n being the file of sequences and nucleotides 
# it goes through the whole file when its called
# For the rest of the for loop where it is called, it runs in O(n^2)

# When the count_kmers function and the count_kmers2 function get called, they run in: O(n+k) because it runs through the sequences in the entire file and then for each k-mer it has to look at the length we inputted
# the function does a slice based on what we input 
# 'k' is either the value of a parameter or the number of elements in the parameter.

# When the count_kmers function gets called it is within a for loop, making the time complexity: O(n(n+k)) = O(n^2+nk)

# When the count_kmers2 function gets called it is also whithin a for loop, making the time complexity also: O(n^2+nk)

# The smaller for loops that print the top 5 could have a run time complexity of O(n) - for indexing through the dictionary 
# and in worst case would have to loop through all of the k-mers in the file

# most_common() time complexity:
# "We use most_common to return k > 1 elements, then we use nlargest method of heapq. 
# This is an O(k) + O((n - k) log k) + O(k log k) algorithm, which is very good for a small constant k, since it's essentialy linear. 
# The O(k) part comes from heapifying the initial k counts, the second part from n - k calls to heappushpop method 
# and the third part from sorting the final heap of k elements. Since k <= n we can conclude that the complexity is:
# O(n log k)" - Stack Overflow 


# Final run time would be approximately O(n^2+nk)+ O(nlog(k)) 
# (ignoring constants)


###################################################################

# Memory Cost of the program

# The memory limitations of this program are associated with the data structures used.
# We use a dictionary to hold k-mers (key) and number of times they appear (value)
# The max size of this dictionary is going to be (4^kmer_Length) with kmer_Length restricted to sizes between 3 and 8
# An array is also used to append the kmer_Count results from each sequence
# The memory limitation comes with the array and not the dictionary. This is because the dictionary will concatonate the same k-mers that appear. 
# The array will keep them all separate. For example, if we were looking for 1-mers, the dictionary would only hold keys of {A, C, T, G}.
# However, the array would hold each occurrence of [A, C, T, G], which would be the entire sequence.
# Therefore, if we input a sequence that is so long, we could run out of memory when appending to the array.

# There are other memory limitations within this program, however, they are insignificant compared to the memory limitations
# used with the data structures. These smaller limitations are those such as: arithmetic operations (+=1 takes less memory than counter = counter + 1)
# and declaring local variables saves more memory than declaring global variables.  

###################################################################


import argparse 				  # argparse = parser for command line
from collections import Counter   # going to using collections and Counter built in python methods to output the top 5 


###################################################################


def usage():      # print the usage of this application if they do not enter the proper arguments in the command line
  
	print ("How to use this program")
	print ("Python CosiJackson_HW2.py -F <file_name> -L <kmer_Length>")
	print ("Where -F is the FASTA file ending in .fasta or .txt to be read and analyzed.")
	print ("Where -L is the length of a k-mer in order to find the most common of that length.")
	print ("The k-mer length must be between numbers 3 and 8.")


###################################################################

# the purpose of this function is to parse the fasta file based on sequence name and sequence (the nucleotides).


def parse_Fasta(file):														# define the function with file as only parameter
	sequence_Name, sequence = None, []										# the two parameters are the sequence name and the sequence. we will be storing the sequence into an array.
	for line in file:														# iterating through the file line by line
		line = line.rstrip()												# removes the new line character
		if line.startswith(">"):											# condition - used for parsing the file into separate sequences, separated by sequence name + sequence
			if sequence_Name: yield (sequence_Name, ''.join(sequence))		# we use yield because we are iterating through a file. yield is a generator which only iterates through the file once in order to parse it.
			sequence_Name, sequence = line, []								# using yield, this function will only return the generator object - sequence name and sequence 
		else:
			sequence.append(line)											# adding the sequence line to the array which holds the sequence string 
	if sequence_Name: yield (sequence_Name, ''.join(sequence))				# this function will only return the generator object - sequence name and sequence 


# How this function and yield works:
# We start by iterating through the file and every line starts the same 
# First, we strip off the new line character and come to our first condition which is to check if there is a '>' 
# If there isn't, we have a sequence. When there is '>', we go into the if statement and sequence_Name = None so we skip the yield. 
# We then get into the else statement where we set sequence_Name = line, so name = >sequence_Name1 and sequence = []
# We do this through out the whole file in order to parse it. Every time we come across a sequence string we append it to the array
# that contains the rest of the sequence associated with it's sequence_Name. 
# yield means return so it returns a tuple of (sequence_Name, sequence)
# Also, the ''.join(sequence) command turned the array into a string, so we have a sequence in the format of a string.
# Once the yield is done, we continously move throught the file and check the conditions and append to an array. 

###################################################################


# the purpose of this function is to find the most common k-mers appearing in all the sequences. 
# we use this function after we run the parse_fasta function.


def count_Kmers(sequences, kmer_Length):				# pass in an array of sequence from the output of the parse function and the length of the kmer we are looking for 
	counts = {}											# Start with an empty dictionary
	num_kmers = len(sequences) - kmer_Length + 1	    # Calculate how many kmers of length k there are
	for i in range(num_kmers):							# Loop over the kmer start positions
		kmer=sequences[i:i+kmer_Length]					# Slice the string to get the kmer
		if kmer not in counts:							# Add the kmer to the dictionary if it's not there
			counts[kmer] = 0							# if the k-mer is not in the dictionary, add it to the dictionary
		counts[kmer] += 1								# Increment the count for this kmer because it is in the dictionary 
	return counts							    		# Return the final counts of the kmer and how many times it appears in it's sequence 

###################################################################

# the purpose of this function is to find the most common k-mers appearing in the most sequences. 
# we use this function after we get the results from the parse_fasta function 
# this function is very similar to the above function, but instead of needing to count each k-mer that appears in the sequence, we only want to know if it appears at all 

def count_Kmers2(sequences, kmer_Length):				# pass in an array of sequence from the output of the parse function and the length of the kmer we are looking for 
	counts2 = {}										# Start with an empty dictionary
	num_Kmers = len(sequences) - kmer_Length + 1	    # Calculate how many kmers of length k there are
	for i in range(num_Kmers):							# Loop over the kmer start positions
		kmer=sequences[i:i+kmer_Length]				    # Slice the string to get the kmer
		if kmer not in counts2:							# Add the kmer to the dictionary if it's not there
			counts2[kmer] = 1							# set the count equal to 1 - because we only care if it appears at least once within the individual sequences 
	return counts2					  		   		    # Return the k-mers and whether or not they appear in a sequence or not.
														# if it appears in a sequence, the dictionary will hold the k-mer name as key and a value of 1
														# when we call this function, we will iterate through to find which k-mers appear in the most sequences  

###################################################################

# Main application function :
# The purpose of the main function is to determine the top 5 most common k-mers that appear in the entire sequence, and the top 5 most common k-mers that appear in the most sequences. 

# The main function first sets the flags used in the command line arguments. 
# It sets each argument to a variable so we can use the input and analyze it. 
# If the user does not correctly input the proper command line argument, the usage function will be called. 
# Next, we parse the file using the parse_fasta function in order to separate every sequence. 
# With every sequence separated, we call the count_kmer function so that each kmer returns a count of how many times it appears in its seqeunce. this is stored in a dictionary. 
# With the individual dictionaries made, we append all of the dictionaries to an array so that all the kmers and # of times they appear are grouped together, rather than individually. 
# We then go through the array and add it to a dictionary and call Counter() on it, which will sort the dictionary in ascending order based on most common k-mer.
# Next, we call the most_common() function, which returns the top 5 most common k-mers. If there are ties for the 5th most common, we print all the ties. 
# For problem 2, the only thing different is the count_kmers2 function, because we only want to know if the k-mer appears in a sequence at least once.
# This program is meant to help a molecular biology student search for the top k-mers that appear in large sequences of nucleotides.


###################################################################


def main():
	parse = argparse.ArgumentParser()            		 # this allows for flags to be in the command line
	parse.add_argument('-F', '-f', '--filename')  		 # the fasta file
	parse.add_argument('-L', '-l', '--length', type=int) # the k-mer length, which will be restricted to lengths 3 thru 8

	arg = parse.parse_args()

	file = arg.filename 		# setting command line arguments to variables to be used within file 
	kmer_Length = arg.length 	

	if(file == None or kmer_Length == None):			 # if the user does not input arguments into the command line correctly 
		print(usage())									 # usage function gets called so we can tell the user how to use the program
		exit(2) 										 # exit 2 so we don't keep going through the file and trying to run something 

	if ((kmer_Length < 3) or (kmer_Length > 8)): 		 # if the user choses a k-mer size less than 3 or greater than 8 we run the print statement because k-mer length is supposed to be between 3 and 8
		print("The length of the most common k-mer is not in the appropriate range. Please choose a number between 3 and 8.")
		exit(2)									 		 # exit 2 so we don't keep going through the file and trying to run something 


	print("\nFile: " + file)                   			 # printing the name of the .fasta or .txt file that it is being searched through        
	print("K-mer Length: " + str(kmer_Length) + '\n')    # printing the kmer length from user input


###################################################################


	# parse by sequence name and store in array 

	sequence_Nucleotides = []							  # make an array of the nucleotides in a sequence based on the output from the parse_fasta (seq) function so that the count_kmer function can go through them 
	with open(file) as f:								  # with open is more efficient because it closes the file on its own python is awesome.
		for sequence_Name, sequence in parse_Fasta(f):	  # call parse_fasta function to separate each sequence in a file 
			sequence_Nucleotides.append(sequence)		  # append the sequence and store into an array 


###################################################################


	# Problem 1
	
	print("\nTop k-mers occurrences in all the sequences:\n")


	all_Sequences = []	
	final_Dict1 = {}				
	final_Dict1 = Counter(final_Dict1)	

	for sequences in sequence_Nucleotides:	  					# sequences is an arbitrary index that will loop through the individual sequences 
		each_Sequence = count_Kmers(sequences, kmer_Length)		# each_sequence will hold the indiviual sequence as they are returned from the count_kmers function 
		all_Sequences.append(each_Sequence)						# we will then append all of the sequences to an array so they will no longer be individual 

	counter = 0				
	output_Kmers = 0		
	
	for index in all_Sequences:							# i is an arbitrary variable used to index through all of the sequences
		final_Dict1 = final_Dict1 + Counter(index)		# we then add the sequences to a dictionary (final_Dict1) and call Counter on the dictionary so that we can find top 5 k-mers 
														# we iterate through each sequence using i from the array, adds it to a dictionary, and then calls counter so that it is sorted is ascending order and the top 5 can be printed 
	for key,value in final_Dict1.most_common():	 		# we will iterate through the keys and values of the dictionary in order to find the top 5 k-mers. most_common() is a built in function that allows us to do this
		if (counter >= 5 and output_Kmers == value):	# if we have printed 5 times or more, but our previous value is still equal to our current value, print the current value 
			print ('%s: %i' % (key,value))				# print the key and value for the k-mers that occur an equal # of times as the 5th k-mer in the top 5
			output_kmers = value 						# continue through the most_common and print if it fulfills the conditions
		elif (counter < 5):								# if we have printed less than 5 k-mers, continue through the most_common and print them 
			print ('%s: %i' % (key,value))				# printing the k-mer name and the value - key and value from dictionary
			counter = counter + 1						# continuing through the most_common 
			output_Kmers = value


###################################################################
	

	# Problem 2 - this problem is very similar to the one above, so a lot of variable names have a 2 added to the end (because we are doing almost the same thing and a 2 for problem 2)
	
	print("\nTop k-mers in the most sequences:\n")


	all_Sequences_Part2 = []
	final_Dict2 = {}						
	final_Dict2 = Counter(final_Dict2)			
	
	for s in sequence_Nucleotides:	  						# s is an arbitrary index used to separate the k-mers
		each_Sequence_Part2 = count_Kmers2(s, kmer_Length)	# for each sequence, we call the second count kmer function (count_kmers2) so that we can know if a k-mer appears within that sequence 
		all_Sequences_Part2.append(each_Sequence_Part2)		# we add each sequence into an array

	counter2 = 0	
	output_Kmers2 = 0
	
	for index in all_Sequences_Part2:						 # indexing through all the sequences
		final_Dict2 = final_Dict2 + Counter(index)			 # calling Counter on the sequences to have them ordered in ascending order of top 5 k-mers
	for key_Part2,value_Part2 in final_Dict2.most_common(5): # print the top 5 k-mers
		print ('%s: %i' % (key_Part2,value_Part2))  		 # I do not use the conditional statement like I did above because in a large file test it printed out too many, so I stop right at 5. 


###################################################################


if __name__ == '__main__':
	main()

###################################################################
