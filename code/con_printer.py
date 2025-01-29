import functions as func

file_name = input(str('What is the file path for the votes file: '))
print('The following is the consensus sequenece of the file.')
print(func.get_consensus(file_name))