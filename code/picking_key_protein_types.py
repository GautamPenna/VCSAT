import functions as func

'''Our Key protein types in the example are listed below: 

protien_type_list = ['L protein ', 'middle protein ', 'preS1/preS2/S envelope protein ', 'MAG: large S protein ', 
                     'truncated large surface protein ', 'small envelope protein ', 'envelope protein ', 'middle surface protein ', 
                     'Middle S protein ', 'truncated large S protein ', 'truncated small S protein ',  'M-HBs ',  'large S protein ', 'preS2/S ',  
                     'Large S protein ', 'MHBs ',  'large HBsAg protein ', 'PreS2 ', 'LHBs ', 'L-HBs ', 'M protein ',  'Envelope preotein ',  
                     'small surface protein ', 'truncated middle S protein ', 'small S protein; HBsAg ', 'surface antigen mRNA preS1/preS2/S ', 
                     'envelop protein ', 'surface antigen precursor ', 'Pre-S1/Pre-S2/S ', 'M-HBsAg ', 'pre-S1/pre-S2/S protein ', 'preS2 protein ', 
                     'pre-S1/pre-S2/S ', 'truncated small surface protein ', 'S protein ', 'Pre-S1/S2/S ',  'L-HBsAg ', 'surface protein ', 
                     'MAG: small S protein ', 'S ',  'preS1 protein ', 'PreS1 protein ', 'env protein ', 
                     'surface antigen HBsAg ', 'large protein ', 'PreS1 ', 'middle S protein ',
                     'elongated middle surface protein ', 'PreS ', 'preS1/preS2/S protein ', 'elongated small surface protein ', 
                     'PreS1/preS2/S ', 'large envelope protein ', 'large S antigen ', 'large S Protein ', 'medium S protein ', 'S-HBs ', 'S-HBsAg ', 
                     'HBV preS1 ', 'pre-S/S protein ', 'pre-S2/S ', 'truncated S protein ', 
                     'small S protein ', 'elongated large surface protein ',  'HBsAg ', 'preS1/preS2/S envelope ', 'middle envelope protein ', 
                     'pre-S2/S protein ', 'PreS1/PreS2/S ', 'MAG: middle S protein ', 'LHBs/MHBs/HBsAg ', 
                     'large surface protein ', 'large S proetin ', 'surface antigen ']
                     
                     '''

## In order to get the above input from the users, we will have to ask them manually until they say they are done.

element = ' '
protein_type_list = []

print('Before picking protein types, make sure that those picked do not contain commas in the name.')

while element != 'done':
    element = input('What protein type do you want to include: ')
    protein_type_list.append(element)

protein_type_list.pop()

output_file = input('Type out the file path you want for the output file from this code, including the name of the file: ')

func.key_protein_lengths(protein_type_list, output_file)