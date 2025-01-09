import functions as func

file_name = input(str('What is the file path for the coding regions file: '))
print(func.determining_protein_types(file_name))
print('By looking at the protein types above, decide which protein types you want to explore more.')
print('Furthermore, a file named All_protein_entry.txt was made if you wish to go further with this code. In order to proceed, this file must be created.')