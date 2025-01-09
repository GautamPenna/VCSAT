import functions as func

del_anum = []
element = ' '
while element!= 'done':
    element = input('What ascension numbers do you want to delete? Type done when you are finished: ')
    del_anum.append(element)

del_anum.pop()

max_lengths_file = input('What max lengths file do you want to use?: ')
output_file = input('Type the pathway of the output file you want: ')

func.delete_acension_ids(max_lengths_file, output_file, del_anum)
print('your requested file has been created.')


#/Users/tompenn/Documents/KE_LAB/first_run_txt_files/(3)A_genotype_max_lengths.txt