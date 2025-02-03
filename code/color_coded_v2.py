import functions as func

genome_votes_file = input('What is the genome votes file you want to use: ')
protein_votes_file = input('What is the protein votes file you want to use: ')
name_of_image = input('What is the name of your graph. Make sure to not include spaces and end the image with a .png: ')

ORF_start = input('What is the position of the nucleotide that starts the ORF of interest: ')
ORF_end = input('What is the position of the nucleotide that ends the ORF of interest: ')

func.color_coded_graph_v2(genome_votes_file, protein_votes_file, name_of_image, ORF_start, ORF_end)

print('your image has been created along with the legend image.')


#G_votes: /Users/tompenn/Desktop/EXAMPLE_Analysis_genotype_A/(3)post-deletion:modification/ADelModvotes.txt
#pvotes: 
#name: /Users/tompenn/Desktop/EXAMPLE_Analysis_genotype_A/pic2.png
