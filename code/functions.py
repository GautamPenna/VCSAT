## PURPOSE: Scan through the raw data and determine the types of protiens that are present as coding regions
## There are a lot of proteins that we do not need for this project specifically, examples including core protein and others.
## After getting an output, place the unique list in a text file and delete the ones that you do not want.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import re
import os
from PIL import Image
from PIL import ImageDraw, ImageFont

def get_genomes_from_alignment_file(alignment_file):
    num_rows = int(input('How many rows are in your votes file: '))
    complete_text = ''
    genomes = []
    ids = []
    with open(alignment_file, 'r') as file:
        for line in file:
            line = line.replace('\n','') #takes away any new lines
            line = line.replace(' ','') #takes away any spaces
            complete_text += line
        first_index = 0
        while first_index != -1:
            first_index = complete_text.find('___')
            id_index = complete_text.find('>')
            genome_pc = complete_text[first_index + 3:first_index + 3 + num_rows + 1] #genome post consensus. The plus one accounts for the fact that python doesnt read the last character
            id_pc = complete_text[id_index + 1:id_index+9]
            genomes.append(genome_pc)
            ids.append(id_pc)
            complete_text = remove_indexes(complete_text, [id_index,first_index, first_index + 1, first_index + 2])
                    

    genomes.pop()
    ids.pop()

    #print(genomes[1])

    return(genomes, ids)

def translate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-', 
        'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            try: protein+= table[codon]
            except: protein += '-'

    return protein 

def remove_indexes(string, indexes):
    result = ""
    for i, char in enumerate(string):
        if i not in indexes:
            result += char
    return result

def determining_protein_types(file_name, all_pros_file_path): #Our example is /Users/tompenn/Documents/KE_LAB/Hep_B_surface_protien_project/(1)coding_regions.txt'
    
    #OLD WAY TO PARSE THROUGH DATA ENTRY LIST, didn't work for HepC due to more > than Hep B.
    
    '''with open(file_name, 'r') as file_one:
        transcript = file_one.read().replace("\n", "") #remove all spaces between entries
        data_points_entry_list = [] #create empty list
        first_index = 0 
        count = 0
        while first_index != -1:
            first_index = transcript.find('>')
            new_transcript = transcript[:first_index] + transcript[first_index + 1:]
            second_index = new_transcript.find('>')
            second_index += 1
            data_points_entry_list.append(transcript[first_index + 1:second_index])
            transcript = new_transcript[second_index-1:]
            count += 1'''


    #new way to parse through and create data_points_entry_list
    # First find |, then go forward to find the next >, go back to find hte first >, then delete this entry.

    with open(file_name, 'r') as file_one:
        transcript = file_one.read().replace("\n", "") #remove all spaces between entries
        data_points_entry_list = [] #create empty list
        first_index = 0 
        while first_index != -1:
            first_index = transcript.find('|')
            #print(first_index)
            gsigns = []
            for i in range(first_index):
                if transcript[i] == '>': 
                    gsigns.append(i)
                    #print(i)
                else: pass
            try: 
                zeroth_index = min(gsigns)
                #print(zeroth_index)
            except: 
                zeroth_index = 0
                #print('Default')
            new_transcript = transcript[first_index + 1:]
            second_index = new_transcript.find('>') + 1
            data_points_entry_list.append(transcript[zeroth_index + 1:second_index])
            transcript = new_transcript[second_index-1:]
        
    file_one.close()

    file = open(all_pros_file_path, 'w')
    file.write('ID,Entry')

    file.write('\n')

    protien_type_list = []
    for entry in data_points_entry_list:
        first_index = entry.find('|')
        second_index = entry.find('[')
        protien_type = entry[first_index + 1:second_index-1]
        protien_type_list.append(protien_type)
        if 'join(' in entry: entry = entry.replace("join(", "")
        else: pass
        id_num = entry[0:10]
        file.write(str(id_num) + ',' + str(entry))
        file.write('\n')


    file.close()   #Cant go through 40K sequences every time, going to make it into a file that I can always look at.

    return list(set(protien_type_list)) # with this command and above for loop, we got the list of protein classifications for surface protiens

def key_protein_lengths(protein_type_list, output_file_name, all_pros_file_path): ### End product was a file that had lengths of proteins and their associated accension #'s labeled in a csv format.
    protien_lengths_file = open(output_file_name, 'w')
    protien_lengths_file.write('ID,Protien_Name,Sequence,Length,Genotype')
    protien_lengths_file.write('\n')

    with open(all_pros_file_path, 'r') as file:
        #print('File Opened!')
        for line in file:
            id_num = line[0:10]
            entry = line[11:]
            name_start = line.find('|')
            name_end = line.find('[')
            name = line[name_start+1:name_end-1]


            name_line = line[name_start:]
            
            display_end = min([name_line.find('['), name_line.find(',')])
            display_name = name_line[1:display_end]



            if 'genotype' in line: 
                ind_genotype = line.find('genotype') + 9
                new_line = line[ind_genotype:]
                ind_end_genotype = new_line.find(']')
                #print(ind_genotype, ind_end_genotype)
                genotype = line[ind_genotype:ind_genotype + ind_end_genotype]
            else: genotype = 'Nada'
            for type_protien in protein_type_list:
                if name==type_protien:
                    start_index = entry.find(']') + 1
                    sequence = entry[start_index:]
                    sequence = sequence.replace('\n','')
                    sequence = sequence.replace(' ','')
                    len_seq = len(sequence)
                    protien_lengths_file.write(str(id_num)+','+str(display_name)+','+str(sequence)+','+str(len_seq)+','+str(genotype))
                    protien_lengths_file.write('\n')

    protien_lengths_file.close()
    #print('Your file has been created in your requested pathway.')

    #Purpose - Look for wanted Protein types and find the associated lengths to find the longest one.
    #Replace the list below with a list that you made




    #need to get its ID again, going most are the same length
    #some are a 1 character shorter length, so I am going to delete those since they don't have the protien we are looking for anyway


    #The Protein_SeqEntry.txt file was made in order to associate an accension # with protein type name.

def finding_longest_protein(last_file_name, genotypes, folder_name):  
    sequence_lengths_df = pd.read_csv(last_file_name)
    #print(genotypes)
    #print(sequence_lengths_df.head(10))
    for element in genotypes: 
        work_df = sequence_lengths_df[sequence_lengths_df['Genotype'] == element] ##### BASED ON WHAT GENOTYPE YOU WANT, CHANGE THAT HERE !!!!!!
        print(work_df.head(10)) # - CANNOT USE N/A bc python pandas doesn't recognize this as an actual genotype but as no entry present here.
        accension_nums = work_df['ID'].tolist()
        unique_an = list(set(accension_nums))  ### Need to do this because each accension # has multiple entries for the various proteins
        for ID in unique_an:
            subset = work_df[work_df['ID'] == ID]
            ID_lens = subset['Length'].tolist()
            len_max = max(ID_lens)
            work_df.loc[(work_df['ID'] == ID) & (work_df['Length'] == len_max), 'max_status'] = 'MAX_VAL'  ### Assigns the ones that have the longest length with a max val status
        
        final_df = work_df[work_df['max_status'] == 'MAX_VAL']
        final_df.to_csv(str(folder_name) + '/' + str(element)+"_genotype_max_lengths.txt", index=False) #UPLOAD THE CSV FILE SO THE RUN TIME ISNT SO LONG

        ## OUTPUT - files of different genotypes that have the longest associated protein length. Each accension number should only be included ONCE by this stage of the code
        ### END OF 3's CODE

        ### Loops through seq_len file that has all protein types - looks for the one that is the longest length and choose that one
        ### We want the one with the longest lenght since this insures that most of the surface protein is included

def consensus_determination_file_creator(input_file, output_file):
    ### Creating the Consensus_Determination file according to the format required by the official MAAFT website, use the code below:
    read_file = open(input_file, 'r')
    final_df = pd.read_csv(read_file)
    with open(output_file, 'w') as file: ### CHANGE THE NAME OF THE FILE WHEN NEEDED
        i=0
        for index, row in final_df.iterrows():
            #file.write('>'+ row['ID'] + ' ' + row['Protien_Name'] + '.' + ' #' + '\n')
            file.write('>' + 'protein_name' + '.#' + '\n') #ONLY USED TO CREATE PHYLOGENETIC TREE
            file.write(row['Sequence'] + '\n')
            i += 1

#make a deletion function that deletes certain acesnion #'s

def delete_acension_ids(max_lengths_file, output_file, del_anum):
    df = pd.read_csv(max_lengths_file)
    for entry in del_anum:
        entry +=  '.1'
        row_index = df.index[df['ID'] == entry].tolist()[0]
        df.drop(row_index, inplace=True)
        print('Deleted Ascension Number: ' + str(entry))
    df.to_csv(output_file, index = False)

def votes_calculator (alignment_file, votes_file): 
    # Genotype specific things you need to know - length of consensus sequence, which should be the same as every sequence

    complete_text = ''
    genomes = []
    with open(alignment_file, 'r') as file: #change name of genotype here
        for line in file:
            line = line.replace('\n','')
            line = line.replace(' ','')
            complete_text += line
        remove_sign = 0
        while remove_sign != -1:
            remove_sign = complete_text.find('>')
            complete_text = complete_text[remove_sign + 1:]
            first_index = complete_text.find('___') # FOR SOME REASON HEP B HAS 4, NEED TO RE-EVAULUTE WHEN DETERMINING GEN D, E
            last_index = complete_text.find('>')
            if last_index == -1:
                genome_pc = complete_text[first_index + 3:] #genome post consensus  ### CHANGE TO 4 FOR HEP B
            else: 
                genome_pc = complete_text[first_index + 3:last_index] #genome post consensus ### CHANGE TO 4 FOR HEP B
            genomes.append(genome_pc)
            #complete_text = complete_text[:first_index] + complete_text[first_index + 4:]
            #complete_text = complete_text[last_index + 1: ]
            #complete_text = remove_indexes(complete_text, [first_index, first_index + 1, first_index + 2, first_index + 3, last_index])
    genomes.pop()

    sumi = 0
    for entry in genomes:
        sumi += len(entry)
    
    length = int(sumi/len(genomes))
    #print(length)

    with open(votes_file, 'w') as file: #Change name of genotype here
        file.write('Position,Consensus_Nucleotide,A_count,T_count,G_count,C_count,Other_count,Total_count_check,Percent_Variability')
        file.write('\n')
        for i in range(0,length): ## INPUT FOR LENGTH OF CONSENSUS SEQUENCE
            A_count = 0
            T_count = 0
            C_count = 0
            G_count = 0
            other_count = 0
            for entry in genomes:
                #print(len(entry))
                if entry[i] == 'A':
                    A_count += 1
                elif entry[i] == 'T':
                    T_count += 1
                elif entry[i] == 'G':
                    G_count += 1
                elif entry[i] == 'C':
                    C_count += 1
                else:
                    other_count += 1
            con_count = max(A_count, T_count, C_count, G_count, other_count)
            con_nuc = ''
            if con_count == A_count: 
                con_nuc = 'A'
                percent_variability = (G_count + T_count + C_count + other_count)/(A_count + G_count + T_count + C_count + other_count)
            if con_count == T_count: 
                con_nuc = 'T'
                percent_variability = (G_count + A_count + C_count + other_count)/(A_count + G_count + T_count + C_count + other_count)
            if con_count == C_count: 
                con_nuc = 'C'
                percent_variability = (G_count + T_count + A_count + other_count)/(A_count + G_count + T_count + C_count + other_count)
            if con_count == G_count: 
                con_nuc = 'G'
                percent_variability = (A_count + T_count + C_count + other_count)/(A_count + G_count + T_count + C_count + other_count)
            if con_count == other_count: 
                con_nuc = '-'
                percent_variability = (G_count + T_count + C_count + A_count)/(A_count + G_count + T_count + C_count + other_count)
            total_count_check = A_count + G_count + T_count + C_count + other_count
            file.write(str(i) + ',' + con_nuc + ',' + str(A_count) + ',' + str(T_count) + ',' + str(G_count) + ',' + str(C_count) + ',' + str(other_count) + ',' + str(total_count_check) + ',' + str(percent_variability*100))
            file.write('\n')

def create_variability_graph(votes_file, title):
    A_votes = pd.read_csv(votes_file) #change file name here to make grapsh of different genotypes

    plt.plot(A_votes['Position'], A_votes['Percent_Variability'], color='skyblue')
    
    plt.xlabel('Position')
    plt.ylabel('Percent Variation')
    plt.title(title)

    plt.show()

def orf_options(votes_file, ORF_match_file, genotype):
    df = pd.read_csv(votes_file)
    consensus_sequence = ''

    for index, row in df.iterrows():
        consensus_sequence += row['Consensus_Nucleotide']
    
    ####### ATG ########
    substring_1 = 'ATG'
    matches_1 = re.finditer(substring_1, consensus_sequence, re.IGNORECASE)
    indicies_1 = [match_1.start() for match_1 in matches_1]
    variabilities_1 = []
    for entry in indicies_1:
        current_var_1 = []
        current_var_1.append(float(df.loc[entry, 'Percent_Variability']))
        current_var_1.append(float(df.loc[entry + 1, 'Percent_Variability']))
        current_var_1.append(float(df.loc[entry + 2, 'Percent_Variability']))
        pos_avg = sum(current_var_1)/len(current_var_1)
        variabilities_1.append(pos_avg)
    mod_3_indicies_1 = []
    for index in indicies_1:
        mod_3_indicies_1.append(int(math.fmod(index,3)))


    ####### TAA ########
    substring_2 = 'TAA'
    matches_2 = re.finditer(substring_2, consensus_sequence, re.IGNORECASE)
    indicies_2 = [match_2.start() for match_2 in matches_2]
    variabilities_2 = []
    for entry in indicies_2:
        current_var_2 = []
        current_var_2.append(float(df.loc[entry, 'Percent_Variability']))
        current_var_2.append(float(df.loc[entry + 1, 'Percent_Variability']))
        current_var_2.append(float(df.loc[entry + 2, 'Percent_Variability']))
        pos_avg = sum(current_var_2)/len(current_var_2)
        variabilities_2.append(pos_avg)
    mod_3_indicies_2 = []
    for index in indicies_2:
        mod_3_indicies_2.append(int(math.fmod(index,3)))

    ####### TAG ########
    substring_3 = 'TAG'
    matches_3 = re.finditer(substring_3, consensus_sequence, re.IGNORECASE)
    indicies_3 = [match_3.start() for match_3 in matches_3]
    variabilities_3 = []
    for entry in indicies_3:
        current_var_3 = []
        current_var_3.append(float(df.loc[entry, 'Percent_Variability']))
        current_var_3.append(float(df.loc[entry + 1, 'Percent_Variability']))
        current_var_3.append(float(df.loc[entry + 2, 'Percent_Variability']))
        pos_avg = sum(current_var_3)/len(current_var_3)
        variabilities_3.append(pos_avg)
    mod_3_indicies_3 = []
    for index in indicies_3:
        mod_3_indicies_3.append(int(math.fmod(index,3)))


    ####### TGA ########
    substring_4 = 'TGA'
    matches_4 = re.finditer(substring_4, consensus_sequence, re.IGNORECASE)
    indicies_4 = [match_4.start() for match_4 in matches_4]
    variabilities_4 = []
    for entry in indicies_4:
        current_var_4 = []
        current_var_4.append(float(df.loc[entry, 'Percent_Variability']))
        current_var_4.append(float(df.loc[entry + 1, 'Percent_Variability']))
        current_var_4.append(float(df.loc[entry + 2, 'Percent_Variability']))
        pos_avg = sum(current_var_4)/len(current_var_4)
        variabilities_4.append(pos_avg)
    mod_3_indicies_4 = []
    for index in indicies_4:
        mod_3_indicies_4.append(int(math.fmod(index,3)))


    ### MATCHES FINDING
    with open(ORF_match_file, 'w') as file:
        file.write('Genotype,Start_index,Stop_index,length,var_avg')
        file.write('\n')
        for i in range(len(mod_3_indicies_1)):
            for j in range(len(mod_3_indicies_2)):
                if mod_3_indicies_1[i] == mod_3_indicies_2[j]:
                    start_index = int(indicies_1[i])
                    stop_index = int(indicies_2[j])
                    length = stop_index - start_index
                    var_avg = 0.5*(variabilities_1[i] + variabilities_2[j])
                    file.write(str(genotype)+','+str(start_index)+','+str(stop_index)+','+str(length)+','+str(var_avg)+'\n')
            for k in range(len(mod_3_indicies_3)):
                if mod_3_indicies_1[i] == mod_3_indicies_3[k]:
                    start_index = int(indicies_1[i])
                    stop_index = int(indicies_3[k])
                    length = stop_index - start_index
                    var_avg = 0.5*(variabilities_1[i] + variabilities_3[k])
                    file.write(str(genotype)+','+str(start_index)+','+str(stop_index)+','+str(length)+','+str(var_avg)+'\n')
            for m in range(len(mod_3_indicies_4)):
                if mod_3_indicies_1[i] == mod_3_indicies_4[m]:
                    start_index = int(indicies_1[i])
                    stop_index = int(indicies_4[m])
                    length = stop_index - start_index
                    var_avg = 0.5*(variabilities_1[i] + variabilities_4[m])
                    file.write(str(genotype)+','+str(start_index)+','+str(stop_index)+','+str(length)+','+str(var_avg)+'\n')

def orf_analysis(orf_file, votes_file, alignment_file, length_condition, variability_condition):
    df = pd.read_csv(orf_file) #change file name here to make grapsh of different genotypes
    df = df[df['length'] > length_condition] ### REMOVING CASES WHERE START CODON IS AFTER STOP CODON
    df = df[df['var_avg'] < variability_condition] ### REMOVING CASES WHERE the variability is too much - setting 5% as the alpha value

    print('Below is a full list of all the ORF matches based on your conditions.')

    test = df.sort_values(['length'], ascending=[0])
    print(test) # Use this print statement to figure out which one you want to explore more

    option = ''
    while option != 'done':
        row_analysis = int(input('Choose which row you want to analyze first: '))
        row = test.iloc[row_analysis]
        print('You chose the ORF displayed above. Below are the variabilities surrounding the stop and start codons')

        start_index = row['Start_index']
        stop_index = row['Stop_index']

        df_votes = pd.read_csv(votes_file)
        print(' ')
        print('START CODON CHARACTERISTICS: ')
        print(df_votes.iloc[[start_index]][['Position', 'Consensus_Nucleotide', 'A_count' , 'T_count' , 'G_count' , 'C_count' , 'Other_count' , 'Percent_Variability']])
        print(df_votes.iloc[[start_index + 1]][['Position', 'Consensus_Nucleotide', 'A_count' , 'T_count' , 'G_count' , 'C_count' , 'Other_count' , 'Percent_Variability']])
        print(df_votes.iloc[[start_index + 2]][['Position', 'Consensus_Nucleotide', 'A_count' , 'T_count' , 'G_count' , 'C_count' , 'Other_count' , 'Percent_Variability']])
        #start_codon_index = int(start_index['Position'])
        print(' ')
        print('STOP CODON CHARACTERISTICS: ')
        print(df_votes.iloc[[stop_index]][['Position', 'Consensus_Nucleotide', 'A_count' , 'T_count' , 'G_count' , 'C_count' , 'Other_count' , 'Percent_Variability']])
        print(df_votes.iloc[[stop_index + 1]][['Position', 'Consensus_Nucleotide', 'A_count' , 'T_count' , 'G_count' , 'C_count' , 'Other_count' , 'Percent_Variability']])
        print(df_votes.iloc[[stop_index + 2]][['Position', 'Consensus_Nucleotide', 'A_count' , 'T_count' , 'G_count' , 'C_count' , 'Other_count' , 'Percent_Variability']])
        #stop_codon_index = int(stop_index['Position'])

        num_rows = df_votes.shape[0]

        continue_1 = input('Would you like to continue with this analysis (Y/N): ')
        if continue_1 == 'Y':
            complete_text = ''
            genomes = []
            ids = []
            with open(alignment_file, 'r') as file:
                for line in file:
                    line = line.replace('\n','')
                    line = line.replace(' ','')
                    complete_text += line
                first_index = 0
                while first_index != -1:
                    first_index = complete_text.find('____')
                    id_index = complete_text.find('>')
                    genome_pc = complete_text[first_index + 4:first_index + 4 + num_rows] #genome post consensus
                    id_pc = complete_text[id_index + 1:id_index+9]
                    #print(id_pc)
                    #print(genome_pc)
                    genomes.append(genome_pc)
                    ids.append(id_pc)
                    complete_text = remove_indexes(complete_text, [id_index,first_index, first_index + 1, first_index + 2, first_index + 3])
                    

            genomes.pop()
            ids.pop()

            inter_file = 'ids_genomes_Type' + str(row['Genotype']) + '.txt'

            with open(inter_file, 'w') as file:
                file.write('id,genome'+'\n')
                for i in range(len(ids)):
                    file.write(str(ids[i]) + ',' + str(genomes[i]) + '\n')
            
            print('An intermediate file was created to aid in our progress. Please keep it where it is made. The file is called: ' + inter_file)

            df_ten = pd.read_csv(inter_file)
            stop_codon = input('From the output printed above, what is the stop codon present on the consensus sequence: ')

            print('Present below are all the ascension nubers that do not agree with the start and stop codon at the assigned values.')
            missing_ATG = [] #Empty list of id nubers that don't have start codon in position 33-35
            missing_stop = []
            dash_ids = []

            for index, row in df_ten.iterrows():
                genome = row['genome']
                start_interest = genome[start_index:start_index + 3]
                stop_interest = genome[stop_index:stop_index + 3]
                
                if start_interest != 'ATG':
                    missing_ATG.append(row['id'])
                    print(row['id'] + str(' :') + start_interest + ' at positions: ' + str(start_index) + ' to ' + str(start_index + 2))
                    if start_interest == '---': dash_ids.append(row['id'])

                if stop_interest != stop_codon:
                    missing_stop.append(row['id'])
                    print(row['id'] + str(' :') + stop_interest + ' at positions: ' + str(stop_index) + ' to ' + str(stop_index + 2)) 
                    if stop_interest == '---': dash_ids.append(row['id'])
            print('Present above are all the ascension numbers that you should note down in a seperate file.')
            print('Decide which ones you want to delete and which alignment sequence you want to modify to continue with the analysis.')

            if os.path.exists(inter_file):
                os.remove(inter_file)
            else:
                print("This should never happen.")
        else: pass

        del_dashes = input('Would you like to delete all start/stop codon values that are "---"? (Y/N): ')
        if del_dashes == 'Y':
            max_lengths_file = input('What max lengths file do you want to use?: ')
            output_file = input('Type the pathway of the output file you want: ')
            element = ' '
            while element!= 'done':
                element = input('What other ascension numbers do you want to delete? Type done when you are finished: ')
                dash_ids.append(element)

            dash_ids.pop()

            delete_acension_ids(max_lengths_file, output_file, dash_ids)
            print('Your new max-lengths file has been created in the requested directory. Please refer to the guide to figure out the next steps.')

        option = input('Do you want to continue to analyze another ORF or are you done? If done, type "done". If not, click enter: ')

        #by the end of this code, you should have written down all the acension #'s that you want to modify and delete for a more conserved analysis.
        #this is where I am currently at in the notes app.
        
### First, if you have ones that you need to delete, delete those and do a re-alignment
### Then, after finding the ORF again, if you have ones you need to change, then change using the output from the alignment MAAFT output files

def alignment_sequences_modifier(alignment_file):

    complete_text = ''
    genomes = []
    ids = []
    
    with open(alignment_file, 'r') as file:
        for line in file:
            line = line.replace('\n','')
            line = line.replace(' ','')
            complete_text += line
        first_index = 0
        id_index = 0
        while id_index != -1:
            id_index = complete_text.find('>')
            id_pc = complete_text[id_index + 1:id_index+9]
            ids.append(id_pc)
            remove_sign = complete_text.find('>')
            complete_text = complete_text[remove_sign + 1:]
            first_index = complete_text.find('____')
            last_index = complete_text.find('>')
            if last_index == -1:
                genome_pc = complete_text[first_index + 4:] #genome post consensus
            else: 
                genome_pc = complete_text[first_index + 4:last_index] #genome post consensus
            genomes.append(genome_pc)
                        

    genomes.pop()
    ids.pop()

    inter_file = 'ids_genomes_Type.txt'

    with open(inter_file, 'w') as file:
        file.write('id,genome'+'\n')
        for i in range(len(ids)):
            file.write(str(ids[i]) + ',' + str(genomes[i]) + '\n')

    df = pd.read_csv(inter_file)
    print(df)
    
    entry = ''
    while entry != 'done':
        ID = input('Which Ascension numbers do you want to modify: ')
        row_index = df.index[df['id'] == ID].tolist()[0]
        genome = df.loc[row_index, 'genome']
        element = ''
        while element != 'done':
            num = int(input('Which position do you want to change: '))
            new_nuc = input('What nucleotide do you want to change it to: ')
            choice = input('It is currently: ' + str(genome[num]) + '. Do you want to replace it with a ' + str(new_nuc) + ' (Y/N): ')
            if choice == 'Y':
                df.loc[row_index, 'genome'] = genome[:num] + new_nuc + genome[num + 1:]
            else: pass
            element = input('Do you want to continue with this Ascension Numer? Once you are done, type "done". To continue, click enter: ')
        entry = input('Do you have any other Ascension Numbers to modify? To continue, click enter. If not, type "done": ')
    

    with open(alignment_file, 'w') as file:
        for i in range(len(ids)):
            file.write('>' + str(df['id'].tolist()[i] + '_1____' + str(df['genome'].tolist()[i]) + '\n'))
    
    if os.path.exists(inter_file): os.remove(inter_file)
    else: print("This should never happen.")

    print('Your old Ascension Alignment file has been modified with the new nucleotides.')

def color_coded_var_graph(votes_file, name_of_image, ORF_start, ORF_end):
    df = pd.read_csv(votes_file)
    row_count = int(ORF_end) - int(ORF_start) + 1 #len(df)
    num = math.sqrt(row_count)
    num += 10
    num = int(round(num,0))
    axis = num*10

    image = Image.new("RGB", (axis, axis + 60), "white") #height and width are the number
    I1 = ImageDraw.Draw(image)
    #font = ImageFont.load_default()
    


    # Add Text to an image

    for i in range(row_count):
        nuc = df.loc[i + int(ORF_start), 'Consensus_Nucleotide']
        var = int(df.loc[i, 'Percent_Variability'])

        bbox = I1.textbbox((0, 0), nuc)#, font='Times New Roman')

        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]   

        

        fill = (255,255,255)
        if var <= 1: fill = (220,220,220)
        elif (var > 1 and var <= 10): fill = (246,189,192) #gradients received from https://www.schemecolor.com/red-gradient.php
        elif (var > 10 and var <= 20): fill = (241,149,155)
        elif (var > 20 and var <= 30): fill = (240,116,112)
        elif (var > 30 and var <= 40): fill = (234,76,70)
        elif (var > 40 and var <= 50): fill = (220,28,19)



        ### AXIS DETERMINATION
        y = int((i - (i % num))/num)*10 + 1
        x = (i % num)*10 + 1

        I1.rectangle([x, y, x + text_width, y + text_height], fill = fill)
        I1.text((x, y), str(nuc), fill=(0,0,0))

        start_color = (255, 255, 255)  # white
        end_color = (255, 0, 0)    # red

        I1.text((20,axis - 20), 'Variability Legend:', fill = (0,0,0))

        lbound = axis-40

        for i in range(lbound):
            r = int(start_color[0] + (end_color[0] - start_color[0]) * i / lbound)
            g = int(start_color[1] + (end_color[1] - start_color[1]) * i / lbound)
            b = int(start_color[2] + (end_color[2] - start_color[2]) * i / lbound)
            I1.line([(i + 20, axis), (i + 20, axis + 15)], fill=(r, g, b))

        I1.text((20,axis + 25), '0%', fill = (0,0,0))
        I1.text((lbound,axis + 25), '49%', fill = (0,0,0))

    
    # Save the edited image
    image.save(name_of_image)

def protein_votes_file (alignment_file, votes_file):

    # Genotype specific things you need to know - length of consensus sequence, which should be the same as every sequence

    complete_text = ''
    genomes = []
    with open(alignment_file, 'r') as file: #change name of genotype here
        for line in file:
            line = line.replace('\n','')
            line = line.replace(' ','')
            complete_text += line
        remove_sign = 0
        while remove_sign != -1:
            remove_sign = complete_text.find('>')
            complete_text = complete_text[remove_sign + 1:]
            first_index = complete_text.find('__')
            last_index = complete_text.find('>')
            if last_index == -1:
                genome_pc = complete_text[first_index + 3:] #genome post consensus
            else: 
                genome_pc = complete_text[first_index + 3:last_index] #genome post consensus
            genomes.append(genome_pc)
            #complete_text = complete_text[:first_index] + complete_text[first_index + 4:]
            #complete_text = complete_text[last_index + 1: ]
            #complete_text = remove_indexes(complete_text, [first_index, first_index + 1, first_index + 2, first_index + 3, last_index])
    genomes.pop()

    #print(genomes)


    
    sumi = 0
    for entry in genomes:
        sumi += len(entry)
    
    length = int(sumi/len(genomes))
    con_nuc = ''
    for i in range (0, length):
        pos_pro = []
        for genome in genomes:
            pos_pro.append(genome[i])
        loop = list(set(pos_pro))
        max_count = 0
        con_pro = ''
        for element in loop:
            count = pos_pro.count(element)
            if count > max_count:
                con_pro = element
                max_count = count
        
        con_nuc += con_pro

    ### Consensus Nucleotide should be created

    with open(votes_file, 'w') as file: #Change name of genotype here
        file.write('Position,Consensus_Nucleotide,Same_Count,Diff_Count,Total_count_check,Percent_Variability')
        file.write('\n')
        for i in range(0,length): ## INPUT FOR LENGTH OF CONSENSUS SEQUENCE
            same_count = 0
            diff_count = 0
            for entry in genomes:
                #print(len(entry))
                if entry[i] == con_nuc[i]:
                    same_count += 1
                else:
                    diff_count += 1

            percent_variability = (diff_count)/(same_count + diff_count)
            total_count_check = same_count + diff_count
            file.write(str(i) + ',' + con_nuc[i] + ',' + str(same_count) + ',' + str(diff_count) + ',' + str(total_count_check) + ',' + str(percent_variability*100))
            file.write('\n')

def protein_translator(alignment_file, ORF_start, ORF_end, new_file_path, plength):
    complete_text = ''
    genomes = []
    ids = []
    with open(alignment_file, 'r') as file:
        for line in file:
            line = line.replace('\n','')
            line = line.replace(' ','')
            complete_text += line
        first_index = 0
        while first_index != -1:
            first_index = complete_text.find('____')
            id_index = complete_text.find('>')
            genome_pc = complete_text[first_index + 4:first_index + 5 + int(plength)] #genome post consensus
            id_pc = complete_text[id_index + 1:id_index+9]
            #print(id_pc)
            #print(genome_pc)
            genomes.append(genome_pc)
            ids.append(id_pc)                    
            complete_text = remove_indexes(complete_text, [id_index,first_index, first_index + 1, first_index + 2, first_index + 3])
                    
    genomes.pop()
    ids.pop()

    orfs_for_all = []
    for entry in genomes:
        orf = entry[ORF_start:ORF_end+1]
        orfs_for_all.append(orf)

    proteins = []
    for entry in orfs_for_all: proteins.append(translate(entry))

    with open(new_file_path, 'w') as file:
        for i in range(len(ids)):
            file.write('>')
            file.write(ids[i] + str('protein_random_name___'))
            file.write(proteins[i])
            file.write('\n')

def color_coded_graph_v2(genome_votes_file, protein_votes_file, name_of_image, ORF_start, ORF_end):
    df_genome = pd.read_csv(genome_votes_file)
    df_protein = pd.read_csv(protein_votes_file)


    row_count = int(ORF_end) - int(ORF_start) + 1 #len(df)
    p_row_count = int((int(ORF_end) - int(ORF_start) + 1)/3)
    buffer = 60
    num = 90 #math.sqrt(row_count) + 10
    #num = int(round(num,0))

    axis = num*10

    image = Image.new("RGB", (2*axis + 30, 2*(axis) + 120), "white") #height and width are the number
    I1 = ImageDraw.Draw(image)

    def fill_color(pos_var,df_genome):
       max_var = 0
       cur_var = 0
       for index, row in df_genome.iterrows():
            cur_var = row['Percent_Variability']
            if max_var < cur_var: max_var = cur_var
       blue_color = int(255 - 255*(pos_var/max_var))
       green_color = int(255 - 255*(pos_var/max_var))
       fill_color = (255, green_color,blue_color)
       return fill_color, max_var

    for i in range(row_count):
        nuc = df_genome.loc[i + int(ORF_start), 'Consensus_Nucleotide']
        var = int(df_genome.loc[i + int(ORF_start), 'Percent_Variability'])

        bbox = I1.textbbox((0, 0), nuc)#, font='Times New Roman')

        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]   

        

        #fill = (255,255,255)
        #if var <= 1: fill = (255,255,255)
        #elif (var > 1 and var <= 10): fill = (246,189,192) #gradients received from https://www.schemecolor.com/red-gradient.php
        ##elif (var > 10 and var <= 20): fill = (241,149,155)
        #elif (var > 20 and var <= 30): fill = (240,116,112)
        #elif (var > 30 and var <= 40): fill = (234,76,70)
        #elif (var > 40 and var <= 50): fill = (220,28,19)

        ### AXIS DETERMINATION
        y = int((i - (i % num))/num)*60 + 1
        x = (i % num)*10 + 1 + buffer

        color, max_var = fill_color(var, df_genome)

        I1.rectangle([x, y, x + text_width-1, y + text_height], fill = color)
        I1.text((x, y), str(nuc) + ' ', fill=(0,0,0))

        if (i%3 == 0): I1.rectangle([x-2,y-2,x + 3*text_width + 7, y + text_height + 50], outline = 'black')
        else: pass

    for j in range(p_row_count):
        pro = df_protein.loc[j, 'Consensus_Nucleotide']
        var = int(df_protein.loc[j, 'Percent_Variability'])

        bbox = I1.textbbox((0, 0), pro)#, font='Times New Roman')

        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]   

        

        #fill = (255,255,255)
        #if var <= 1: fill = (200,200,200)
        #elif (var > 1 and var <= 10): fill = (246,189,192) #gradients received from https://www.schemecolor.com/red-gradient.php
        #elif (var > 10 and var <= 20): fill = (241,149,155)
        #elif (var > 20 and var <= 30): fill = (240,116,112)
        #elif (var > 30 and var <= 40): fill = (234,76,70)
        #elif (var > 40 and var <= 50): fill = (220,28,19)

        #start_fill = (255,255,255)
        

        ### AXIS DETERMINATION
        y = int((j - (j % (num/3)))/(num/3))*60 + 21
        y_nuc = y + 20
        x = (j % (num/3))*30 + 11 + buffer

        color, max_var2 = fill_color(var, df_protein)

        string = str(pro)

        I1.rectangle([x-9, y-3, x + 14, y + 14], fill = color)
        I1.text((x, y), string, fill=(0,0,0))
        if j+1 < 10: #aestheic letter number pairing.
            I1.text((x + 1, y_nuc), str(j + 1), fill=(0,0,0))
        elif j+1 < 100:
            I1.text((x-2, y_nuc), str(j + 1), fill=(0,0,0))
        else: I1.text((x-5, y_nuc), str(j + 1), fill=(0,0,0))
 
    start_color = (255, 255, 255)  # white
    end_color = (255, 0, 0)    # red

    I1.text((buffer,y_nuc + 30), 'Variability Legend:', fill = (0,0,0))

    lbound = axis-40

    for i in range(lbound):
        r = int(start_color[0] + (end_color[0] - start_color[0]) * i / lbound)
        g = int(start_color[1] + (end_color[1] - start_color[1]) * i / lbound)
        b = int(start_color[2] + (end_color[2] - start_color[2]) * i / lbound)
        I1.line([(i + buffer, y_nuc + 50), (i + buffer, y_nuc + 65)], fill=(r, g, b))

    I1.text((buffer,y_nuc + 75), '0%', fill = (0,0,0))
    I1.text((lbound + 50,y_nuc + 75), str(int(max_var)) + '%', fill = (0,0,0))

    I1.text((0,1), 'Nucleotide: ', fill = (0,0,0))
    I1.text((0,21), '    Protein: ', fill = (0,0,0))
    I1.text((0,41), '   Position: ', fill = (0,0,0))

    
    
    # Save the edited image
    image.save(name_of_image)

def get_consensus(votes_file):
    df = pd.read_csv(votes_file)

    consensus_sequence = ''
    for index, row in df.iterrows():
        consensus_sequence += row['Consensus_Nucleotide']

    return(consensus_sequence)

def gaps_in_seq(genome_votes_file, alignment_file, deletion_threshold):
    con_seq = get_consensus(genome_votes_file)
    df_votes = pd.read_csv(genome_votes_file)
    positions = []
    for pos in range(len(con_seq)): 
        if con_seq[pos] == '-': positions.append(pos) #talking about indexes
    num_rows = df_votes.shape[0]

    complete_text = ''
    genomes = []
    ids = []
    with open(alignment_file, 'r') as file:
        for line in file:
            line = line.replace('\n','')
            line = line.replace(' ','')
            complete_text += line
        first_index = 0
        while first_index != -1:
            first_index = complete_text.find('____')
            id_index = complete_text.find('>')
            genome_pc = complete_text[first_index + 4:first_index + 4 + num_rows] #genome post consensus
            id_pc = complete_text[id_index + 1:id_index+9]
            genomes.append(genome_pc)
            ids.append(id_pc)
            complete_text = remove_indexes(complete_text, [id_index,first_index, first_index + 1, first_index + 2, first_index + 3])
                    

    genomes.pop()
    ids.pop()

    adt = [] #above deletion threshold positions
    bdt = [] #below deletion threshold positions
    print('These are the positions of consensus sequence that have "_" in the sequence: ' + str(positions))
    for i in positions:
        var = float(df_votes.loc[i, 'Percent_Variability'])
        if var >= float(deletion_threshold): adt.append(i)
        else: bdt.append(i)
    
    print('The positions below the deletion threshold are: ' + str(bdt))

    bdt_acids = []
    for pos in bdt:
        for i in range(len(genomes)):
            genome = genomes[i]
            if genome[pos] != '-': 
                bdt_acids.append(ids[i])
                print('The Accession Number ' + str(ids[i]) + ' does not have a "-" in position ' + str(pos))

    if len(adt) > 0:
        print('The positions above the deletion threshold are: ' + str(adt))
        for pos in adt:
            pos_row = df_votes[df_votes['Position'] == pos]
            print('Here are the variabilities for position ' + str(pos))
            print(pos_row)
            switch_nuc = input('Which nucleotide would you like to switch the consensus sequence with? Type "N" if you do not want to switch with any: ')
            if switch_nuc == 'A':
                df_votes.loc[df_votes['Position'] == pos, "Consensus_Nucleotide"] = 'A'
            elif switch_nuc == 'G':
                df_votes.loc[df_votes['Position'] == pos, "Consensus_Nucleotide"] = 'G'
            elif switch_nuc == 'T':
                df_votes.loc[df_votes['Position'] == pos, "Consensus_Nucleotide"] = 'T'
            elif switch_nuc == 'C':
                df_votes.loc[df_votes['Position'] == pos, "Consensus_Nucleotide"] = 'C'
            
        

    
        





    else: print('There are no positions that have a variability difference higher than your threshold value.')

def seq_var(votes_file, start_pos, end_pos):
    ## This function serves to determine the variability of a select portion of the sequence by using the variabilities of the nucleotide in the bounds of the sequence

    seq_range = int(end_pos) - int(start_pos) + 1

    df = pd.read_csv(votes_file)
    var_list = []
    for index, row in df.iterrows():
        var_list.append(float(row['Percent_Variability']))
    
    selected_range = []

    for i in range(seq_range):
        selected_range.append(var_list[start_pos - 1 + i])

    seq_var = sum(selected_range)/len(selected_range)
    
    print('The variability of your selected sequence is: ' + str(seq_var) + '%.')

def gap_analysis(votes_file, post_alignment_file):
    con_seq = get_consensus(votes_file)
    df_votes = pd.read_csv(votes_file)
    positions = []
    for pos in range(len(con_seq)): 
        if con_seq[pos] == '-': positions.append(pos) #talking about indexes
    num_rows = df_votes.shape[0]

    genomes, ids = get_genomes_from_alignment_file(post_alignment_file)
    positions = []
    var_list = []

    replacement_amino_acids = []

    for index, row in df_votes.iterrows():
        nuc = row['Consensus_Nucleotide']
        if nuc == '-':
            var_list.append(row['Percent_Variability'])
            positions.append(int(row['Position']))
        else: pass
    
    for pos in positions:
        pos_dict = dict()
        pos_string = ''
        for genome in genomes:
            pos_string += genome[pos]
        unique_chars = list(set(pos_string))

        for char in unique_chars: pos_dict.update({char: int(pos_string.count(char))})
        #print(pos_dict)
        values_sorted = sorted(pos_dict.values(), reverse=True)
        max_val_2 = values_sorted[1]
        keys = [key for key, value in pos_dict.items() if value == max_val_2][0] #only thing wrong with this method is what if there are multiple substitions with same value? Need to account for that later.
        
        print('In position: ' + str(pos) + ' there is a gap present. The second best substitution is: ' + str(keys) + '. There are ' + str(pos_dict.get(keys)) + ' to support this substitution.')
        replacement_amino_acids.append(keys)
    
    continue_maybe = input('Would you like to replace all of these positions with the recommended amino acid (Y/N): ')
    if continue_maybe == 'Y':
        file_path = input('Type the pathway for the new alignment file to be stored: ')
        file_path_2 = input('Type the pathway for the new votes file to be stored: ')
        protein_name = input('What is the name of the protein you are analyzing: ')
        # You have lists of geomes, ids, positions and replacement amino acids
        for i in range(len(positions)):
            for j in range(len(genomes)):
                genomes[j] =  genomes[j][:positions[i]] + replacement_amino_acids[i] + genomes[j][positions[i] + 1:]
        

        with open(file_path, 'w') as file:
            for i in range(len(ids)): file.write('>' + str(ids[i]) + protein_name + '___' + str(genomes[i]) + '\n')
        file.close()

        print('Your new alignement file has been created in the specified directory.')
        
        protein_votes_file(file_path, file_path_2)
        
        print('Your new votes file has been created in the specified directory.')

    else: pass

def cross_protein_matching_1(all_pros_file, p1_alignment_file, output_file): #### NOT TRASH BUT REMEMBER THAT THE ACCESSION NUMBERS ARE DIFFERENT, THE CORE AND E1/E2 WILL NEVER BE THE SAME
    print('This function is used to link together 2 proteins that you want to align together.')
    element = ' '
    p2_list = []

    while element != 'done':
        element = input('What protein type do you want to include: Make sure to get the name exactly. ')
        p2_list.append(element)

    p2_list.pop()
    inter_file = 'inter_file.txt'
    key_protein_lengths(p2_list, inter_file, all_pros_file)

    genomes, ids = get_genomes_from_alignment_file(p1_alignment_file)
    p2_df = pd.read_csv(inter_file)
    p2_ids = list((p2_df['ID']))
    p2_seqs = p2_df['Sequence']

    print(ids)
    print(p2_ids)

    combined_id = []

    display_name = input('What do you want the display name to be: ')
    
    with open(output_file, 'w') as file2:
        file2.write('ID,Sequence,Length')
        for i in range(len(ids)):
            for j in range(len(p2_ids)):
                if (ids[i] + '.1') == (p2_ids[j]): 
                    print('Found one!')
                    sequence = genomes[i] + p2_seqs[j]
                    ID = ids[i]
                    length = len(sequence)
                    file2.write(str(ID) + ',' + str(sequence) + ',' + str(length))
                    file2.write('\n')

def cross_protein_addition(p1_sequence, p2_key_lengths_file, output_file):
    df = pd.read_csv(p2_key_lengths_file)
    name = input('What is the name of your combined protein (no spaces): ')
    for index, row in df.iterrows():
        df.loc[index, 'Sequence'] = p1_sequence + row['Sequence']
        df.loc[index, 'Protien_Name'] = name
        df.loc[index, 'Length'] = len(row['Sequence'] + p1_sequence)
    df.to_csv(output_file)
    print('your file was created in the specified directory')

def full_seq_analysis():
    full_seq = input('What is full sequence you want to analyse: ')
    start = int(input('What is start position you want to analyse: '))
    end = int(input('What is end position you want to analyse: '))

    print('Your request sequence is: ' + str(full_seq[start:end + 1]))



