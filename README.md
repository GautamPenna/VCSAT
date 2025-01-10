# VCSAT - Viral Consensus Sequence Analysis Tool
*Python-based tool for analyzing viral genomes*

This respository contains code files and a step-by-step documentation on how to use VCSAT for Viral Genome Sequence Analysis. Genome downloads from NCBI can be downloaded in a variety of formats, but to parse through the data and determine a *consensus sequence* is difficult. This tool guides its users through a series of steps, starting from downloading the coding regions of requested proteins on NCBI to graphical analysis when compared against the consensus or reference sequences. 

> *Consensus Sequence: each nucleotide position is determined by the most common nucleotide across all data points*

### ✍️ Author
Gautam Penna, Ke Lab, University of Texas at Austin

**Methodology Overview**:
1.	Downloading sequences from NCBI
2.	Determining Protein Types present in Data
3.	Picking Protein Types
4.	Determining Longest Protein for every Accession Number (ID)
5.	Creating Consensus Determination Files
6.	MAFFT Alignment
7.	Variability Calculation
8.	Graphical representation of Variability and Position
9.	ORF Match File Creation
10.	ORF Analysis
11.	Ascension Number Deletion
12.	Genome Modifier
13.	Color Coded Variability Graph

## General Outline of Sequence Determination

The following steps determine the variability among the consensus sequence so it can be later used in AlphaFold and other 3d protein modeling softwares to determine structural variability of a virus. The first couple steps in this process involve downloading the required genomes from NCBI and determining the proteins required for analysis. The longest protein of each ID number (Accession Number) is used for further analysis since using incomplete DNA sequences will lead to multiple ‘-’ during the MAFFT alignment process. After creating the files required by the MAFFT software, variability between the determined consensus sequence and each individual is calculated. 

From here, it is important to analyze the ORF (open reading frame) of the start and stop codons. Finding ones that are in the same reading frame that also encompass the largest protein is crucial for further analysis. After determining the Accession Numbers that do not align with the ORF determined, deletion of these ID’s from the dataset and a run-through of steps 6-13 is required for accurate depiction of the protein.

**Example Files**

To 


## Downloading Code for Usage

In order to use the code present /code file, download the respository to your personal computer. After doing so, open your terminal window and navigate to the same directory as the respository. Once reached, type the following command:

```py

python VCSAT.py

```

This should launch the main python file that lets you navigate the options offered by this program.

## (1) Downloading NCBI Viral Protein Genomes


After choosing the virus to analyze from the left-hand side as well as selected filters, following the following steps in the ‘Download’ button on the top left-hand corner: 

Download -> Coding Regions -> All Records -> Default Naming Convention.
The output of this download is going to be a .fasta File. Open the file in .txt format and save that as a new file. This is going to be the starting file for the rest of the functions as part of the program.

**(2) Determining Protein Types Present in NCBI Fasta File**

Run VCSAT.py code on your terminal window. When prompted to, pick option (1). When prompted to input the specific file, copy and paste the entire pathway the .txt file from above. 

The output of this command is going to be a list of all the various naming conventions for protein types that are present in the FASTA download. Copy and paste this list in a different document and pick out the protein types that you wish to include for further analysis. A file named ‘All_protein_entry.txt’ will be made. DO NOT DELETE this file if you wish to proceed with the code.

**(3) Picking Protein Types to Analyze**

After choosing wish protein types to analyze further, run VCSAT.py code on your terminal window. When prompted to, pick option (2). When prompted to, paste all the various protein types that you would like to analyze further. Make sure to paste them in the same format as they are present in the list. If spaces and other miscellaneous characters are present in the list, make sure to type it out in the exact same manner. When prompted to, type out the file path for this code’s output file. In the output file, you will see only the protein types that you selected and their associated lengths. This is important for later in the program.


**(4) Determining Longest Protein Sequences**

Once finished with the previous step, it is time to pick the genotypes you want to analyze as well as finding the longest. Run VCSAT.py code on your terminal window. When prompted to, pick option (3) and go along with the code. A file for each genotype will be created in your desired folder.

**(5) Creating Files for Consensus Sequence Determination**

At this point of the program, you have one protein from each Accession number that represents the longest form of the proteins you want to analyze. To determine the consensus sequence, we will use the RIMD (Research Institute for Microbial Diseases at Osaka University) MAFFT Tool. In order to use this tool effectively, we need to place our sequences in a desired order and format. Run VCSAT.py code on your terminal window. When prompted to, pick option (4) and go along with the code. The consensus determination files will be created in your desired directory. Make sure to run this code one time for each genotype.

**(6) Using MAFFT Server**

In order to determine consensus sequence, use the following link: MAFFT Alignment Link. 

The bottom left pictures represent what should be selected with circles around the important portions in the menu options.

<img width="700" alt="image" src="https://github.com/user-attachments/assets/d0aaad7e-54d7-4e28-b6d8-db0d38c484d3" />

<img width="700" alt="image" src="https://github.com/user-attachments/assets/48c020d0-57ac-422a-be9f-1f45c362ec4a" />

After running the MAFFT command, results that resemble the pictures on the upper right should be made available to view. Select Export Alignment in FASTA → Convert FASTA to .txt file. To be easy to access, it is recommended to name this in the following convention: ‘(genotype)_consensus_file.txt’. 

**(7) Votes Calculation and Visualization**

In order to determine the variability of each position, we will use the post-MAFFT file in order to keep the lengths of all genomes the same. Run VCSAT.py code on your terminal window. When prompted to, pick option (6) and go along with the code. This will create a new file with percent variability and the different numbers of A, T, G and C for each position. In order to visualize the variability, run VCSAT.py code on your terminal window again. When prompted to, pick option (7) and go along with the code. This should create a file that represents the variability of the consensus sequence at each position of the nucleotide.

**(8) ORF Analysis**

After determining the variability of the entire consensus genome, it is important to find an ORF (Open Reading Frame) that includes both a start and stop codon. This truly represents a full protein sequence since it implies that the genetic code can be transcribed into a functional protein. In order to determine the most likely ORF, run VCSAT.py on your terminal window. When prompted to, pick option (8) and go along with the code. Create the file in the desired directory and run VCSAT.py. This time, pick option (9) and go along with the code. This allows you to analyze various options available for ORFs. When selecting likely candidate ORFs, make sure to account of the length of the desired protein and overall variability of the positions. More variability implies having a smaller sample size to perform analysis. 

**(9) Post-ORF Data Cleaning**

If at this stage of the code, there are likely Accession numbers that need to be deleted as well as genomes that could be modified due to the variability of the stop codon. To delete Accession numbers, run VCSAT.py on your terminal window. When prompted to, pick option (5) and go along with the code. The input for this code would be the max_lengths file that is the result at the end of Step 4. After deleting the said Accession numbers, run through steps 4-6. After getting the output of this step, run VCSAT.py on your terminal window. When prompted to, pick option (10) and go along with the code. Edit the sequences are required and continue with steps 7-8.


**(10) Creating Color Coded Variability Map**

To create a color-coded map of the variabilities present in the ORF chosen, run VCSAT.py on your terminal window. When prompted to, pick option (11) and go along with the code.
