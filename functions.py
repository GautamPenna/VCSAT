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



def remove_indexes(string, indexes):
    result = ""
    for i, char in enumerate(string):
        if i not in indexes:
            result += char
    return result

def determining_protein_types(file_name): #Our example is /Users/tompenn/Documents/KE_LAB/Hep_B_surface_protien_project/(1)coding_regions.txt'
    with open(file_name, 'r') as file:
        transcript = file.read().replace("\n", "")
        data_points_entry_list = []
        first_index = 0
        count = 0
        while first_index != -1:
            first_index = transcript.find('>')
            new_transcript = transcript[:first_index] + transcript[first_index + 1:]
            second_index = new_transcript.find('>')
            second_index += 1
            data_points_entry_list.append(transcript[first_index + 1:second_index])
            transcript = new_transcript[second_index-1:]
            count += 1

    protien_type_list = []
    count = 0 #Just for there to be a constant output in processing the code
    for entry in data_points_entry_list:
        first_index = entry.find('|')
        second_index = entry.find('[')
        protien_type = entry[first_index + 1:second_index]
        protien_type_list.append(protien_type)
        count += 1
        #print(count) #constant output reading
    
    file = open('All_protein_entry.txt', 'w')
    file.write('ID,Entry')
    file.write('\n')

    for entry in data_points_entry_list:
        id_num = entry[0:10]
        file.write(str(id_num) + ',' + str(entry))
        file.write('\n')

    file.close()   #Cant go through 40K sequences every time, going to make it into a file that I can always look at.

    return(list(set(protien_type_list))) # with this command and above for loop, we got the list of protein classifications for surface protiens

def key_protein_lengths(protein_type_list, output_file_name): ### End product was a file that had lengths of proteins and their associated accension #'s labeled in a csv format.
    protien_lengths_file = open(output_file_name, 'w')
    protien_lengths_file.write('ID,Protien_Name,Sequence,Length')
    protien_lengths_file.write('\n')

    with open('All_protein_entry.txt', 'r') as file:
        for line in file:
            id_num = line[0:10]
            entry = line[11:]
            name_start = line.find('|')
            name_end = line.find('[')
            name = line[name_start+1:name_end]
            for type_protien in protein_type_list:
                if name==type_protien:
                    start_index = entry.find(']') + 1
                    sequence = entry[start_index:]
                    sequence = sequence.replace('\n','')
                    sequence = sequence.replace(' ','')
                    len_seq = len(sequence)
                    protien_lengths_file.write(str(id_num)+','+str(name)+','+str(sequence)+','+str(len_seq))
                    protien_lengths_file.write('\n')

    protien_lengths_file.close()
    print('Your file has been created in your requested pathway.')

    #Purpose - Look for wanted Protein types and find the associated lengths to find the longest one.
    #Replace the list below with a list that you made




    #need to get its ID again, going most are the same length
    #some are a 1 character shorter length, so I am going to delete those since they don't have the protien we are looking for anyway


    #The Protein_SeqEntry.txt file was made in order to associate an accension # with protein type name.

def finding_longest_protein(last_file_name, genotypes, folder_name):
    sequence_lengths_df = pd.read_csv(last_file_name)
    for element in genotypes: 
        work_df = sequence_lengths_df[sequence_lengths_df['Genotype'] == element] ##### BASED ON WHAT GENOTYPE YOU WANT, CHANGE THAT HERE !!!!!!
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
        for index, row in final_df.iterrows():
            file.write('>'+ row['ID'] + ' ' + row['Protien_Name'] + '.' + ' #' + '\n')
            file.write(row['Sequence'] + '\n')

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
        first_index = 0
        while first_index != -1:
            first_index = complete_text.find('____')
            genome_pc = complete_text[first_index + 4:first_index + 1207] #genome post consensus
            genomes.append(genome_pc)
            complete_text = remove_indexes(complete_text, [first_index, first_index + 1, first_index + 2, first_index + 3])
    genomes.pop()

    with open(votes_file, 'w') as file: #Change name of genotype here
        file.write('Position,Consensus_Nucleotide,A_count,T_count,G_count,C_count,Other_count,Total_count_check,Percent_Variability')
        file.write('\n')
        for i in range(0,1203): ## INPUT FOR LENGTH OF CONSENSUS SEQUENCE
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
                con_nuc = '_'
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

            for index, row in df_ten.iterrows():
                genome = row['genome']
                start_interest = genome[start_index:start_index + 3]
                stop_interest = genome[stop_index:stop_index + 3]
                
                if start_interest != 'ATG':
                    missing_ATG.append(row['id'])
                    print(row['id'] + str(' :') + start_interest + ' at positions: ' + str(start_index) + ' to ' + str(start_index + 2))

                if stop_interest != stop_codon:
                    missing_stop.append(row['id'])
                    print(row['id'] + str(' :') + stop_interest + ' at positions: ' + str(stop_index) + ' to ' + str(stop_index + 2)) 
            print('Present above are all the ascension numbers that you should note down in a seperate file.')
            print('Decide which ones you want to delete and which alignment sequence you want to modify to continue with the analysis.')

            if os.path.exists(inter_file):
                os.remove(inter_file)
            else:
                print("This should never happen.")
        else: pass

        option = input('Do you want to continue to analyze another ORF or are you done? If done, type "done". If not, click enter: ')

        #by the end of this code, you should have written down all the acension #'s that you want to modify and delete for a more conserved analysis.
        #this is where I am currently at in the notes app.
        
### First, if you have ones that you need to delete, delete those and do a re-alignment
### Then, after finding the ORF again, if you have ones you need to change, then change using the output from the alignment MAAFT output files

def alignment_sequences_modifier(alignment_file, length):

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
            genome_pc = complete_text[first_index + 4:first_index + 4 + int(length)] #genome post consensus
            id_pc = complete_text[id_index + 1:id_index+9]
            genomes.append(genome_pc)
            ids.append(id_pc)
            complete_text = remove_indexes(complete_text, [id_index,first_index, first_index + 1, first_index + 2, first_index + 3])
                        

    genomes.pop()
    ids.pop()

    inter_file = 'ids_genomes_Type.txt'

    with open(inter_file, 'w') as file:
        file.write('id,genome'+'\n')
        for i in range(len(ids)):
            file.write(str(ids[i]) + ',' + str(genomes[i]) + '\n')

    df = pd.read_csv(inter_file)
    
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

def color_coded_var_graph(votes_file, name_of_image):
    df = pd.read_csv(votes_file)
    row_count = len(df)
    num = math.sqrt(row_count)
    num += 10
    num = int(round(num,0))
    axis = num*10

    image = Image.new("RGB", (axis, axis), "white") #height and width are the number
    
    I1 = ImageDraw.Draw(image)
    font = ImageFont.load_default()
    


    # Add Text to an image

    for i in range(row_count):
        nuc = df.loc[i, 'Consensus_Nucleotide']
        var = int(df.loc[i, 'Percent_Variability'])

        bbox = I1.textbbox((0, 0), nuc, font=font)

        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]   

        

        fill = (255,255,255)
        if var <= 1: fill = (255,255,255)
        elif (var > 1 and var <= 10): fill = (255,0,0)
        elif (var > 10 and var <= 20): fill = (0,255,0)
        elif (var > 20 and var <= 30): fill = (0,0,255)
        elif (var > 30 and var <= 40): fill = (255,255,100)
        elif (var > 40 and var <= 50): fill = (255,0,255)



        ### AXIS DETERMINATION
        y = int((i - (i % num))/num)*10 + 1
        x = (i % num)*10 + 1

        I1.rectangle([x, y, x + text_width, y + text_height], fill = fill)
        I1.text((x, y), str(nuc), fill=(0,0,0))
        

    
    # Save the edited image
    image.save(name_of_image)





    
    