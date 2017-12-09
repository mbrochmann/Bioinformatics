#!/usr/bin/python

import sys
import time

nucleotide_list = ["A", "C", "G", "T"]

nucleotide_dict = {"A":0, "C":1, "G":2, "T":3}

nucleotide_dict_rev = {v: k for k, v in nucleotide_dict.items()}

nucleotide_dict_bin = {"A":0b00, "C":0b01, "G":0b10, "T":0b11}

nucleotide_dict_bin_rev = {v: k for k, v in nucleotide_dict_bin.items()}

reverse_complement_dict = {"A":"T", "C":"G", "G":"C", "T":"A"}

def PatternToNumber(pattern):
    pattern_length = len(pattern)
    number = 0
    for i in range(0, pattern_length):
        number += nucleotide_dict[pattern[i]] * 4**(pattern_length-i-1)
    return number

def PatternToNumberBin(pattern):
    pattern_length = len(pattern)
    number = 0b00
    for i in range(0, pattern_length):
        number = number << 2
        number += nucleotide_dict_bin[pattern[i]]
    return int(number)

def NumberToPattern(number, length):
    pattern = ""
    #while (number > 0):
    for i in range(0, length):
        pattern = nucleotide_dict_rev[number % 4] + pattern
        number = number // 4
    return pattern

def NumberToPatternBin(number, length):
    pattern = ""
    #while (number > 0):
    for i in range(0, length):
        pattern = nucleotide_dict_rev[number & 0b11] + pattern
        number = number >> 2
    return pattern

# A hashed dictionary is faster
# than calling PatternToNumber each time
def MakePatternToNumberDict(k):
    PatternToNumberDict = {}
    for i in range(0, 4**k):
        PatternToNumberDict[NumberToPatternBin(i,k)] = i
    return PatternToNumberDict

# A hashed dictionary is faster
# than calling NumberToPattern each time
def MakeNumberToPatternDict(k):
    NumberToPatternDict = {}
    for i in range(0, 4**k):
        NumberToPatternDict[i] = NumberToPatternBin(i,k)
    return NumberToPatternDict

def GenerateFrequencyArray(length):
    frequency_array = {}
    for i in range(0, 4**length):
        frequency_array[NumberToPatternBin(i,length)] = 0
    return frequency_array

def ReverseComplement(pattern):
    pattern_rc = ""
    pattern_length = len(pattern)
    for i in range(0, pattern_length):
        pattern_rc = reverse_complement_dict[pattern[i]] + pattern_rc
    return pattern_rc


def PatternCount(input_file, search_string, start_char = 0):

#input_file = str(sys.argv[1])
#search_string = str(sys.argv[2])
    search_string_length = len(search_string)
    
    start_char = 0
    
    if(len(sys.argv)==4):
        start_char = int(sys.argv[3])
        
    match_count = 0;
        
    # IMPORTANT:
    #
    # This algorithm does not count appearances of search_string
    # that are broken by a newline
    #
    # start_char allows the search to start at some offset 
    # from the beginning of the file
    # useful if it is known that the search_string does not 
    # appear before start_char
    
    with open(input_file) as f:
        for line in f.readlines():
            line_length = len(line)
            if(start_char >= line_length):
                start_char -= line_length
                continue
            for i in range(start_char, line_length):
                for j in range(0, search_string_length):
                    if (line[i+j] == search_string[j]):
                        if (j+1 == len(search_string)):
                            match_count += 1
                            sys.stdout.write('%i ' % i)
                    else: 
                        break
    sys.stdout.write('\n')
    #print "match_count = %i" % match_count
    return match_count


def FrequentWordsFromFile(input_file, kmer_length):

    # key = pattern, value = count
    #FrequentPatterns = {}
    FrequentPatterns = GenerateFrequencyArray(kmer_length)
    
    with open(input_file) as f:
        for line in f.readlines():
            line_length = len(line)
            for i in range(0, line_length-kmer_length):
                kmer = line[i:i+kmer_length]
                if kmer in FrequentPatterns:
                    FrequentPatterns[kmer] += 1
                else:
                    FrequentPatterns[kmer] = 1

    last_value = -1
    for key, value in sorted(FrequentPatterns.items(), key = lambda kv:(kv[1],kv[0])):
        print("%s: %i" % (key, value))
        last_value = value

    sys.stdout.write('\n')
    sys.stdout.write('\n')

    print("The following are the most frequent patterns, each appearing %i times" % last_value)
    for key, value in FrequentPatterns.items():
        if(value == last_value):
            sys.stdout.write('%s ' % key)
            #print(key)
    sys.stdout.write('\n')
    sys.stdout.write('\n')

    for key, value in FrequentPatterns.items():
        sys.stdout.write('%i ' % value)
        
    sys.stdout.write('\n')
    sys.stdout.write('\n')


# a version that operates on a string and 
# and returns the frequency_array as an array, not a dictionary
# it is the first step in the ClumpFinding function
def FrequentWords(line, kmer_length, PatternToNumberDict):

    # key = pattern, value = count
    #FrequentPatterns = {}
    #FrequentPatterns = GenerateFrequencyArray(kmer_length)
    FrequentPatterns = [0]*4**kmer_length
    
    line_length = len(line)
    for i in range(0, line_length-kmer_length):
        kmer = line[i:i+kmer_length]
        FrequentPatterns[PatternToNumberDict[kmer]] += 1

    return FrequentPatterns


def ClumpFinding(Genome, k, L, t):
    """ 
    Function to find "clumps", or presence of t or more occurrences of a pattern
    of length k in a window of size L in Genome. 

    clump_array has one entry for each possible pattern of length k
    at the end of the function, each entry is true if we have found
    one or more clumps of the corresponding pattern in the genome 

    Genome is the text string representing the genome
    k is the pattern length
    t is the minimum number of times we want a pattern to be repeated in the window
    L is the window size
    """

    start_time = time.time()

    PatternToNumberDict = MakePatternToNumberDict(k)
    #NumberToPatternDict = MakeNumberToPatternDict(k)

    print(time.time() - start_time)

    #clump_array = GenerateFrequencyArray(k)
    clump_array = [False]*4**k

    num_clumping_kmers = 0

    with open(Genome) as f:
        for line in f.readlines():
            line_length = len(line)
            window = line[0:L]
            frequency_array = FrequentWords(window, k, PatternToNumberDict)
            # check for clumps in frequency_array
            for i in range(0, 4**k):
                if frequency_array[i] >= t:
                    clump_array[i] = 1
            print(time.time() - start_time)
            # now move the window one character at a time
            # decrementing the count for the first pattern in the old window
            # and incrementing the count for the last pattern in the new window
            #indexL = PatternToNumberDict[line[0:k]]
            indexR = PatternToNumberDict[line[L-k:L]]
            #indexR = PatternToNumber(line[L-k:L])
            for i in range(1, line_length-L):
                indexL = PatternToNumberDict[line[i-1:i+k-1]]
                frequency_array[indexL] -= 1
                #indexL = indexL & ((0b10 << (k << 1)) - 0b01 - (0b11 << ((k<<1) - 2)))
                #indexL = indexL << 2
                #indexL += nucleotide_dict_bin[line[i+k-1]]
                #indexL += nucleotide_list.index(line[i+k-1])
                indexR = PatternToNumberDict[line[i+L-k:i+L]]
                #indexR = indexR & ((0b10 << (k << 1)) - 0b01 - (0b11 << ((k<<1) - 2)))
                #indexR = indexR << 2
                #indexR += nucleotide_dict_bin[line[i+L-1]]
                #indexR += nucleotide_list.index(line[i+L-1])
                frequency_array[indexR] +=1
                updated_frequency = frequency_array[indexR]
                # check if recently read kmer forms clump in current frequency_array
                if updated_frequency >= t:
                    if not clump_array[indexR]:
                        clump_array[indexR] = True
                        num_clumping_kmers += 1
                        #print(line[i-1:i+k-1])

    print(time.time() - start_time)

    """
    for i in range(0, 4**k):
        if clump_array[i] == 1:
            num_clumping_kmers += 1
            print(NumberToPatternDict[i])
    """

    print(num_clumping_kmers)

    return clump_array
