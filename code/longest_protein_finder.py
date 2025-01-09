import functions as func

last_file_name = input('What is the pathway of the last file you made in Code 2: ')
genotypes = []
element = ' '
while element!= 'done':
    element = input('What genotypes do you want to include? Type done when you are finished: ')
    genotypes.append(element)

folder_name = input('What is the pathway for the folder you want these files to be stored in: ')

genotypes.pop()

func.finding_longest_protein(last_file_name, genotypes, folder_name)

print('Your files have been created in your desired folder')

