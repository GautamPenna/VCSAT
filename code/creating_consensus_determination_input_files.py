import functions as func

input_file = input('What text file from code 3 do you want to make a consensus determination file from? Provide the entire pathway: ')
output_file = input('What is the pathway for the output_file: ')

func.consensus_determination_file_creator(input_file, output_file)

print('your files have been created in their requested directory')

