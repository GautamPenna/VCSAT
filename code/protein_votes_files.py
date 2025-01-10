import functions as func

alignment_file = input('What is the pathway for the alignment file: ') #/Users/tompenn/Desktop/PEX/Acon1.txt
votes_file = input('What do you want the votes file to be called with pathname: ')

func.protein_votes_file(alignment_file, votes_file)
print('your file has been created')