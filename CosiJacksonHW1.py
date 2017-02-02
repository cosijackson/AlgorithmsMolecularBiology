###################################################################

# File : CosiJacksonHW1.py

# The purpose of this homework is to : 

# use command line parameters to control the behavior of the application
# read a sequence file in FASTA format
# find a specific sequence to be processed within the file
# calculate the %GC content across an entire sequence / specified sequence
# search the sequence for exact matches to a given k-mer / specific k-mer (use sequence ATG as sample)
# report those matches using an annotation file format (GFF)

# Developer : Cosi Jackson
#             CSCI 4314, Spring 2017, HW1


# Referencing - Okeson_Sample.py // Developer : Alex Okeson // Spring 2016 - Formatting of the header, main function, and comments
# Referencing - https://www.tutorialspoint.com/python/ - Python built in functions and notation 


###################################################################


#   Sample command line arguments to run program: 
#   python CosiJacksonHW1.py -F file.fasta.txt -S chrI -K ACA
#   -F specifies a .fasta  or .txt input file to be read in and analyzed
#   -S specifies the sequence to be analyzed
#   -K specifies the kmer to be found within a sequence and number of occurrences within a sequence


###################################################################


import sys, getopt, math, os.path, argparse

# sys = system specific parameters and functions - useful for command line operations
# getopt = parser for command line
# math = mathematical functions
# os.path = common pathname manipulations 
# argparse = parser for command line

###################################################################


def usage():      # print the usage of this application
  
  print "How to use this program"
  print "\n python CosiJacksonHW1.py -F <file_name> -S <seq_name> -K <kmer>" 
  print "\n where -F is the FASTA file ending in .fasta or .txt to be read and analyzed."
  print "\n where -S is the sequence name in the file to be analyzed."
  print "\n where -K is the k-mer (ex. ATA) to be matched within the sequence."


###################################################################

# the purpose of this function is to calculate the GC percent within a specified sequence

def gc_Count(sequence):
  gc_counter = 0                      # set gc_counter equal to zero to initialize
  for x in sequence:                  # going through every nucleotide within the sequence
    if(x == "G" or x == "C"):         # if the nucleotide is a G or a C, add one to the gc_counter and continue through sequence 
      gc_counter = gc_counter + 1     # += ... saves a step in assembly language
  gc_percent = int((float(gc_counter) / float(len(sequence))) * 100) # simple calculation of gc_percent using floats for the decimal numbers but returning the value as an integer so it is not super long 
  return(gc_percent)                  # return gc_percent


###################################################################

# the purpose of this function is to search for a specific k-mer within a specified sequence and store the positions in which the kmer occurs

def kmer_Search(sequence, kmer):
  end = len(sequence)                         # initializing variables
  beg = 0
  kmer_matches = []                           # storing the k-mer match positions within a list 
  while sequence.find(kmer, beg, end) != -1:  # python built in method of .find - looking within the sequence to find the k-mer, starting at beginning and ending at end of seq. 
    x = sequence.find(kmer, beg , end)        # storing every position where k-mer is found in arbitrary varaible (x) so I can append each position into list
    kmer_matches.append(x+1)                  # appending positions from variable (x) into list # we do not want to start at 0 so we add 1 (class rule)
    beg = x+1                                 # (x+1) because x is the position where we find the kmer and we do not want to repeat / find it again in that same posititon 
  return(kmer_matches)                        # return the list of k-mer matches

# find method works by scanning the list until it finds something so we have to increment 
# with x+1 so that we start at position where we last found it and avoid repeats

# add size of kmer to position they were found to get where the kmer position ends


###################################################################

# the purpose of this function is to determine if a kmer has been found within the reverse complement of the sequence

def reverse_Comp(kmer):
  # kmer = kmer.strip('\n')                            # strip the kmer - not sure if this is necessary
  comp_Dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}     # dictionary uses keys and values (every A becomes a T, etc.)
  result = []                                          # creating an array called result

  for nucleotide in kmer:                              # iterating through each nucleotide in the specified kmer
    result.append(comp_Dict[nucleotide])               # append to results the complement of the nucleotide

  result = reversed(result)                            # built in python function to create the reverse complement string
  result = ''.join(result)                             # taking the array and making it into a string
  return result

###################################################################


# Main application function :
# The main function first sets the flags used in the command line arguments. 
# It sets each argument to a variable so we can use the input from the user and analyze it. 
# If the user does not correctly input the proper command line argument, the usage function will be called. 
# Next, we obtain the specified sequence from the specified file that the user gave in the command line.
# The user will also input a specified k-mer to be analyzed within the sequence.
# We analyze the specified sequence by printing the length of the sequence and the %GC content within the sequence. 
# We also analyze the specified k-mer within the sequence by printing the position(s) in which the k-mer occurs in GFF format.
# We count positions starting at 1 rather than 0. 
# This program is meant to help a molecular biology student search for specific k-mers in large sequences of nucleotides.
# This program will output the data in GFF format.


###################################################################


def main():

  parse = argparse.ArgumentParser()             # this allows for flags to be in the command line
  parse.add_argument('-F', '-f', '--filename')  # the fasta file
  parse.add_argument('-S', '-s', '--sequence')  # specified sequence input
  parse.add_argument('-K', '-k', '--kmer')      # specific k-mer to match sequence 

  arg = parse.parse_args()
  # arg, unknown = parse.parse_known_args(['--filename', '--sequence', '--kmer'])
  

  file = arg.filename
  sequence_input = arg.sequence
  kmer = arg.kmer.upper()         # added the upper function in case kmer input is typed with lower case letters 
  program_name = parse.prog



  #if (len(sys.argv) < 4):   # (0 .. CosiJacksonHW1.py) + (1 .. file) + (2 .. sequence) + (3 .. k-mer) = 4 arguments
    # usage()
    # exit (2)                # we have to exit the usage function
  

  # the following matches the sequence name given from the arguments to the sequence 

  nucleotides=[]
  #flag = False
  with open(file,'r') as f:                       # opening the .fasta or .txt file to be analyzed
    for line in f:                                # beginning loop
      if line.startswith('>' + sequence_input):
        #flag = True
        #if flag:
        f.readlines                               # reading the lines of the sequence
        for line in f:
          condition = line.startswith('>')        # going through the sequence and using the condition of the '>' as a stopping point
          if line.startswith('') != condition:    # if there is no '>' it will append the nucleotides of the sequence to be analyzed
            nucleotides.append(line.strip('\n'))  # we use .strip so the '\n' character is not appended to the nucleotide list 
          else:                                   # if there is a '>' it will stop appending the nucelotides of the sequence because the sequence has ended 
              #flag = False   
            break                                 # break to end because we have found our specified sequence

  sequence = ''.join(nucleotides).upper()           # we input our list of nucelotides into a string and call the list of nucelotides a new variable - sequence 
                                                    # also using upper function in case the kmer input is lower case letters

# output print statements

  print("\nFile: " + file)                               # printing the name of the .fasta or .txt file that it is being searched through        
  print("Sequence: " + sequence_input)                   # printing the sequence name that the user has inputted
  print("K-mer: " + kmer + '\n')                         # printing the kmer than is being searched for from user input
  print("Seq Length: " + str(len(sequence)))             # printing the sequence length
  print("GC content: " + str(gc_Count(sequence))) + "%"  # printing the GC% of that sequence 


# output GFF format

# output GFF for non-reverse complement sequences


  positions_list = kmer_Search(sequence, kmer)           # putting the output of kmer_Search function into variable positons_list, which holds each position the specified kmer has been found 
  
  if (len(positions_list) != 0):                         # checking that list of positions is not empty 
    for k in range(len(positions_list)):                 # this loop will go through each element (position) in the list 
      print(str(sequence_input) + "\t" + str(program_name) + "\tMatch\t" + str(positions_list[k]) + "\t" + str(positions_list[k]+len(kmer) - 1) + "\t+\t.")    # for every position in the list, it will print out the sequence, program name, match, starting and ending position of kmer found, and a "+" which signifies that it was found in forward sequence
      k = k + 1       # iterate to the next position within the list
  else: print("\nThis k-mer was not found in the forward sequence. Please try a different k-mer or change the specified sequence.\n")       # if positions list is empty, this print statment will run


# output GFF format for reverse complement sequences

  reverse = reverse_Comp(kmer)    # putting the output of reverse_Comp function into variable reverse, which will print out all reverse complements of specified kmers in the sequence


  positions_list2 = kmer_Search(sequence, reverse) # putting the output of kmer_Search function (which now takes in reverse as an argument) into variable positons_list

  if (len(positions_list2) != 0):                  # checking that list of positions for reverse complement kmer matches is not empty
    for k in range(len(positions_list2)):          # this loop will go through each element (position) in the list
      print(str(sequence_input) + "\t" + str(program_name) + "\tMatch\t" + str(positions_list2[k]) + "\t" + str(positions_list2[k] + len(kmer) -1) + "\t-\t.\n")   # for every position in the list, it will print out the sequence, program name, match, starting and ending position of kmer found, and a "-" which signifies that it was found in reverse complement of sequence
      k = k + 1         # iterate to the next position within the list
  else: print("\nThis k-mer was not found in the reverse complement sequence. Please try a different k-mer or change the specified sequence.\n")       # if positions list is empty, this print statment will run




#  calls the main() function which will run the program

if __name__ == '__main__': 
  main()




