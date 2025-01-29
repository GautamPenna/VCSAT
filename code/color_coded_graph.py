import functions as func

votes_file = input('What is the votes file you want to use: ')
name = input('What is the name of your graph. Make sure to not include spaces and end the image with a .png: ')

ORF_start = input('What is the position of the nucleotide that starts the ORF of interest: ')
ORF_end = input('What is the position of the nucleotide that ends the ORF of interest: ')

func.color_coded_var_graph(votes_file, name, ORF_start, ORF_end)

print('your image has been created along with the legend image.')


# /Users/tompenn/Desktop/txtfiles/E_votes_new.txt
# /Users/tompenn/Documents/KE_LAB/Hep_B_surface_protien_project/first_run_txt_files/(6)A_genotype_votes.txt
# /Users/tompenn/Documents/KE_LAB/Hep_B_surface_protien_project/first_run_txt_files/(6)D_genotype_votes.txt