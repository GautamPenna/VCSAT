import functions as func

votes_file = input('What do you want the votes file to be called with pathname: ')
orf_file = input('What is the pathway for the orf file you want: ')
genotype = input('What is the genotype of this analysis: ')

func.orf_options(votes_file, orf_file, genotype)

print('your file has been created')