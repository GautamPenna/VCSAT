import functions as func

orf_file = input('What is the orf_matching file that you created in step 8: ')
votes_file = input('What is the votes file that you created in step 6: ')
alignment_file = input('Please input the path file you received from the MAAFT alignment server: ')
length_condition = int(input('What the minimum length of your genome you need: '))
variability_condition = int(input('What is the maximum average variability you can consider: '))

func.orf_analysis(orf_file, votes_file, alignment_file, length_condition, variability_condition)

#/Users/tompenn/Documents/KE_LAB/first_run_txt_files/(8)ORF_Matching.txt

#/Users/tompenn/Documents/KE_LAB/first_run_txt_files/(6)A_genotype_votes.txt

#/Users/tompenn/Documents/KE_LAB/first_run_txt_files/(5)A_genotype_alignment.txt