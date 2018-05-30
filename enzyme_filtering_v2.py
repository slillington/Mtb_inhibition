
"""
This script should accomplish the following:
1. Generate a list of InChi strings for all metabolites in a genome scale model
2. Translate those InChi strings into Mol objects in rdkit
3. Compute pairwise Tanimoto similarities between all metabolites
4. Compile BRENDA information on kinetic constants, inhibitors, and inhibitory constants from BRENDA
5. Identify lethal knockouts for Mtb under aerobic growth conditions
6. Compile a table of already observed inhibitors with their Ki's (if available) for the substrates of these lethal reactions

"""
import string
import hashlib
import pandas as pd
from SOAPpy import SOAPProxy ## for usage without WSDL file
from getInfoFromUniprot import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

def ecListbyOrg(organisms):
    
    #Compare these EC lists to those pulled from BRENDA
    #orgs = ["Mycobacterium tuberculosis","Staphylococcus aureus","Escherichia coli","Candida albicans","Homo sapiens"]
    EC_by_orgs = []
    wsdl = "http://www.brenda-enzymes.org/soap/brenda_server.php"
    pw = hashlib.sha256("synbiorox").hexdigest()
    client = SOAPProxy(wsdl)
    for o in organisms:
        parameters = "stephenlillington2017@u.northwestern.edu,"+pw+",#organism*"+o
        resultString = client.getOrganism(parameters)

        #Separate by the various entries
        entries = resultString.split('!')
    
        #Extract the EC numbers
        ecNums = []
        for e in entries:
            i = e.find("ecNumber*")
            h = e.find("#",i+15)
            ecNums.append(e[i+9:h])
        
        ecNums = list(set(ecNums))
        EC_by_orgs.append(ecNums) #Note some entries are byte strings and others are unicode strings

    #Write EC numbers to files
    file = open('C:\Github\pythonScripts\Mtb_inhibition\Mtb_inhibition\data\ec_numbers.txt','w')
    for elist in EC_by_orgs:
        file.write("%%%%%%%%%%%%%%%%NEW ORGANISM%%%%%%%%%%%%%%%%%%" + "\n\n\n")
        for a in elist:
            file.write(a + "\n")
    
    
    
    return EC_by_orgs


def main():
    #organisms = ["Mycobacterium tuberculosis","Staphylococcus aureus","Escherichia coli","Homo sapiens"]
    #ec_lists = ecListbyOrg(organisms)
    
    
    #Read in EC lists from excel file to pandas dataframe
    ec_df = pd.read_excel('C:\GitHub\pythonScripts\Mtb_inhibition\Mtb_inhibition\data\ec_numbers_list')
    print(ec_df)
    
    
    
    '''
    #Identify EC numbers that are common to both Mtb and humans
    mtb_list = ec_lists[0]
    hs_list = ec_lists[3]
    mtb_hs_common = [x for x in mtb_list if x in hs_list]
    filtered_targets = [y for y in mtb_list if not y in hs_list]
    filtered_info = [fgetInfoFromUniprot(z,"Mycobacterium+tuberculosis") for z in filtered_targets]
    del filtered_targets
    
    #Pull Uniprot information for those proteins. If the % identity (score/length of alignment)is below X, send Mtb dict to results
    #What sort of alignment algorithm to use? Start with +1 for match, -1 for mismatch, 0 for gap.
    matrix = matlist.blosum62
    for prot in mtb_hs_common:
        mtb_info = fgetInfoFromUniprot(prot,"Mycobacterium+tuberculosis")
        hs_info = fgetInfoFromUniprot(prot,"Homo+sapiens")
        
        if not hs_info['Sequence'] or hs_info['Sequence'] == []: #if Uniprot search for homo sapiens returns nothing
                filtered_targets.append(mtb_info)
                print("WARNING: " + hs_info['ECnumber'] + " has no sequence in Uniprot")
                
        elif not mtb_info['Sequence'] or mtb_info['Sequence'] == []: #if Uniprot search for mtb returns nothing
            print("WARNING: Mtb enzyme " + mtb_info['ECnumber'] + " has no sequence in Uniprot")
                
        else:       
            mtb_seq = mtb_info['Sequence']

            #print(mtb_seq)
            hs_seq = hs_info['Sequence']

            #print(hs_seq)
            
            #Perform sequence alignment
            tb_hs_alignment = pairwise2.align.globaldx(mtb_seq, hs_seq, matrix,score_only=True)
            tb_tb_alignment = pairwise2.align.globaldx(mtb_seq, mtb_seq, matrix,score_only=True)
            homology = tb_hs_alignment/tb_tb_alignment
            #print(homology)
            
            if homology < 0.42:
                filtered_info.append(mtb_info)
    #Need to change for Uniprot queries
    #organisms = ["Mycobacterium+tuberculosis","Staphylococcus+aureus","Escherichia+coli","Homo+sapiens"]
    
    #getInfoFromUniprot(elist,organisms)
    '''
    
if __name__ == "__main__":
    main()