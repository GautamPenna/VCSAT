import os

print('Welcome to Viral Consensus Sequence Analysis Tool (VCSAT)!')
print('Here, you will be able to determine the consensus sequence of any number of data points for a specfic virus')

option = ' '
while option != 'done':
    print('Our options are: ')
    print('    (01) Determining Protein Types in Coding Region Output File from NCBI')
    print('    (02) Picking Protein Types to Analyze')
    print('    (03) Determining the Longest Protein from Genome')
    print('    (04) Creating Files for Consensus Sequence Determination')
    print('    (05) Deleting Certain Ascension IDs')
    print('    (06) Votes Calculation')
    print('    (07) Creating an Alignment Graph')
    print('    (08) ORF Match File Creator')
    print('    (09) ORF Analyzer')
    print('    (10) Genome Modifier')
    print('    (11) Color Coded Graph Creator')
    print('    (12) Protein Votes Calculator')
    print('Type "done" if you are finished with using VCSAT.')



    option = input('Please choose which option you want to proceed with: ')

    if option == '01':
        os.system('python code/determining_protein_types_1.py')
    elif option == '02':
        os.system('python code/picking_key_protein_types.py')
    elif option == '03':
        os.system('python code/longest_protein_finder.py')
    elif option == '04':
        os.system('python code/creating_consensus_determination_input_files.py')
    elif option == '05':
        os.system('python code/deleting_ascension_nums.py')  
    elif option == '06':
        os.system('python code/votescalc.py')  
    elif option == '07':
        os.system('python code/alignemnt_graphs.py') 
    elif option == '08':
        os.system('python code/orf_options.py') 
    elif option == '09':
        os.system('python code/orf_analysis.py') 
    elif option == '10':
        os.system('python code/genome_modifier.py') 
    elif option == '11':
        os.system('python code/color_coded_graph.py')
    elif option == '12':
        os.system('python code/protein_votes_files.py')

print('Thank you for using our code. Come back soon!')
