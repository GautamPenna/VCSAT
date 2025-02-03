import functions as func

alignment_file = input('What is the alignment file you want to use: ') #/Users/tompenn/Desktop/EXAMPLE_Analysis_genotype_A/(3)post-deletion:modification/ADelCON1.txt

ORF_start = int(input('What is the position of the nucleotide that starts the ORF of interest: ')) 
ORF_end = int(input('What is the position of the nucleotide that ends the ORF of interest: '))
plength = input('What is the length of the entire protein - last number present on votes file: ')

new_file = input('What is the pathway + name you want to use for your protein alignment file: ') #/Users/tompenn/Desktop/EXAMPLE_Analysis_genotype_A/PCon1.txt

func.protein_translator(alignment_file, ORF_start, ORF_end, new_file, plength)

print('Your translated file has been created in the described directory.')