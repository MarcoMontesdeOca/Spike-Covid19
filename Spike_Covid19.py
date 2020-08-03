## First Part: Data extraction
"""
Takes as input a GenPept file obtained from https://www.ncbi.nlm.nih.gov/protein/ after including the search terms:
"surface glycoprotein" and "Severe acute respiratory syndrome coronavirus 2". Such file can be downloaded from the "Send to:" link
at the results page. Once downloaded, the file extension should be manually changed from ".gp" to ".txt".

Outputs a .csv file including protein id, "country of origin", and "collection date" for further amino acid sequence 
extraction and formatting in the Second Part
"""

# Open GenPept full list
file = open("spike protein sequences.txt","r")
raw_string = file.read()
file.close()

# Extract country and collection date information
store = dict()
proteins = raw_string.split("LOCUS")
for protein in proteins:
    country_i = protein.find("country")
    if country_i != -1:     
        country = protein[country_i + 9 : country_i + 30]
        country = country.replace('"','')
        country = country.split("\n")
        country = country[0]
        col_i = protein.find("collection_date")
        col_date = protein[col_i + 17 : col_i + 27]        
        col_date = col_date.replace('"','')
        if col_date[:4] != "2020":
            if col_date[-3:] == "202":
                col_date = col_date[:7] + "2020"
            elif col_date[-3:] == "201":
                col_date = col_date[:7] + "2019"
        col_j = col_date.find("\n")
        if col_j != -1:
            col_date = col_date.split("\n")[0]
        id_i = protein.find("VERSION")
        pro_id = protein[id_i : id_i + 23]
        pro_id = pro_id.replace("VERSION"," ")
        pro_id = pro_id.replace(" ","")
        pro_id = pro_id.split("\n")
        pro_id = pro_id[0]
        source = [country, col_date]
        store[pro_id] = source

# Store results as .csv file
import csv
with open('protein_source.csv', mode='w') as csv_file:
    fieldnames = ['protein id', 'country', 'collection date']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    writer.writeheader()
    for key in store:
        protein = store[key]
        writer.writerow({'protein id': key, 'country': store[key][0], 'collection date': store[key][1]})

## Second Part: Data Formatting
"""
Takes as input .csv file including "protein id", "country of origin", and "collection date".
Extracts amino acid sequence and outputs Multi Fasta File.
"""

# Open accession id csv file
from Bio import Entrez, SeqIO
import csv
file_name = 'protein_source.csv' 
Ids = dict()
reader_object = open(file_name)
for row in csv.DictReader(reader_object): 
    date = row['collection date']
    origin = row['country']
    origin = origin.replace(" ", "_")
    id_ = row['protein id']
    sample = (date, origin)
    Ids[id_] = sample
reader_object.close()

# Extract Spike protein sequences from NCBI    
proteins = dict()  
Entrez.email = "emontesdeocas@gmail.com"
for key in Ids:
    protein_id = key   
    if protein_id != 'NA':
        handle = Entrez.efetch(db="protein", id=protein_id,
                               rettype="fasta")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        protein = str(record.seq)    
        new_key = str(Ids[key][0]) + '_' + Ids[key][1] + '_' + protein_id
        proteins[new_key] = protein

# Store all sequences as a Multi Fasta file
def writeseq(filename, Tdict):
    text_file = open(filename, 'a')
    for key in Tdict:        
        seq = str(Tdict[key])
        text_file.write(f'>{key}\n')
        text_file.write(seq + '\n')
    text_file.close()
writeseq("spike_samples.fasta", proteins)
