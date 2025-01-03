import functions as func

alignment_file = input('Type the path of the alignment file: ')
length = input('What is the length of the consensus sequence. You can determine this by looking at the votes files: ')


func.alignment_sequences_modifier(alignment_file, length)

