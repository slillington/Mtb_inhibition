
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
#import hashlib
import pandas as pd
#from SOAPpy import SOAPProxy ## for usage without WSDL file
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

def filterByOrg(other_organism1,other_organism2=None,other_organism3=None):
#Read in EC lists from excel file to pandas dataframe
    ec_df = pd.read_excel(r"C:\GitHub\pythonScripts\Mtb_inhibition\Mtb_inhibition\data\ec_numbers_lists.xlsx")
    
    
    
    #Identify EC numbers that are common to both Mtb and humans
    mtb_list = list(ec_df['Mycobacterium tuberculosis'].dropna())
    other_list = list(ec_df[other_organism1].dropna())
    
    mtb_common = [x for x in mtb_list if x in other_list]
    filtered_targets = [y for y in mtb_list if not y in other_list]
    filtered_info = [fgetInfoFromUniprot(z,"Mycobacterium+tuberculosis") for z in filtered_targets]
    del filtered_targets
    
    #Pull Uniprot information for those proteins. If the % identity (score/length of alignment)is below X, send Mtb dict to results
    #What sort of alignment algorithm to use? Start with +1 for match, -1 for mismatch, 0 for gap.
    matrix = matlist.blosum62
    other_miscount = 0
    mtb_miscount = 0
    for prot in mtb_common:
        mtb_info = fgetInfoFromUniprot(prot,"Mycobacterium+tuberculosis")
        other_info = fgetInfoFromUniprot(prot,"Homo+sapiens")
        
        if not other_info['Sequence'] or other_info['Sequence'] == []: #if Uniprot search for homo sapiens returns nothing
            other_miscount = other_miscount+1
            #Check if also no sequence for Mtb
            if not mtb_info['Sequence'] or mtb_info['Sequence'] == []:
                filtered_info.append(mtb_info)
                
                
                
        elif not mtb_info['Sequence'] or mtb_info['Sequence'] == []: #if Uniprot search for mtb returns nothing
            #print("WARNING: Mtb enzyme " + mtb_info['ECnumber'] + " has no sequence in Uniprot")
            mtb_miscount = mtb_miscount+1
            #Write EC number to a txt file?
                
        else:       
            mtb_seq = mtb_info['Sequence']

            #print(mtb_seq)
            other_seq = other_info['Sequence']

            #print(hs_seq)
            
            #Perform sequence alignment
            try:
                tb_other_alignment = pairwise2.align.globaldx(mtb_seq, other_seq, matrix,score_only=True)
                tb_tb_alignment = pairwise2.align.globaldx(mtb_seq, mtb_seq, matrix,score_only=True)
                homology = tb_other_alignment/tb_tb_alignment
                #print(homology)
            
                if homology < 0.30:
                    filtered_info.append(mtb_info)
                    
            except:
                print(mtb_seq)
                print(other_seq)
                
    #Write the EC number, Protein name, gene name, and sequence to an excel file
    filtered_EC = [x['ECnumber'] for x in filtered_info]
    filtered_prot = [x['Protein Name'] for x in filtered_info]
    filtered_gene = [x['Gene'] for x in filtered_info]
    
    file = open(r"C:\Github\pythonScripts\Mtb_inhibition\Mtb_inhibition\data\filter_noHuman.tsv",'w')
    for i in range(0,len(filtered_EC)-1):
        file.write(filtered_EC[i] + '\t' + filtered_prot[i] + '\t' + filtered_gene[i] + '\n')
    file.close()                
    print("Number of EC numbers in Homo sapiens not matched with sequence: " + str(other_miscount))
    print("Number of EC numbers in Mycobacterium tuberculosis not matched with sequence: " + str(mtb_miscount))
    print(len(filtered_info))

    #print(filtered_info)
    
    #Write filtered_info to an excel file
    
    ##############################Now do Mtb, Human, Staphylococcus aureus############################################
    filtered_EC = [x['ECnumber'] for x in filtered_info]
    if not other_organism2 == None:
        other_list = list(ec_df[other_organism2].dropna())
        mtb_common = [x for x in filtered_EC if x in other_list]
        filtered_targets = [y for y in filtered_EC if not y in other_list]
        filtered_info = [fgetInfoFromUniprot(z,"Mycobacterium+tuberculosis") for z in filtered_targets]
        del filtered_targets
        
        #Pull Uniprot information for those proteins. If the % identity (score/length of alignment)is below X, send Mtb dict to results
        #What sort of alignment algorithm to use? Start with +1 for match, -1 for mismatch, 0 for gap.
        matrix = matlist.blosum62
        other_miscount = 0
        mtb_miscount = 0
        for prot in mtb_common:
            mtb_info = fgetInfoFromUniprot(prot,"Mycobacterium+tuberculosis")
            other_info = fgetInfoFromUniprot(prot,"Staphylococcus+aureus")
            
            if not other_info['Sequence'] or other_info['Sequence'] == []: #if Uniprot search for homo sapiens returns nothing
                other_miscount = other_miscount+1
                #Check if also no sequence for Mtb
                if not mtb_info['Sequence'] or mtb_info['Sequence'] == []:
                    filtered_info.append(mtb_info)
                    
                    
                    
            elif not mtb_info['Sequence'] or mtb_info['Sequence'] == []: #if Uniprot search for mtb returns nothing
                #print("WARNING: Mtb enzyme " + mtb_info['ECnumber'] + " has no sequence in Uniprot")
                mtb_miscount = mtb_miscount+1
                #Write EC number to a txt file?
                    
            else:       
                mtb_seq = mtb_info['Sequence']

                #print(mtb_seq)
                other_seq = other_info['Sequence']

                #print(hs_seq)
                
                #Perform sequence alignment
                try:
                    tb_other_alignment = pairwise2.align.globaldx(mtb_seq, other_seq, matrix,score_only=True)
                    tb_tb_alignment = pairwise2.align.globaldx(mtb_seq, mtb_seq, matrix,score_only=True)
                    homology = tb_other_alignment/tb_tb_alignment
                    #print(homology)
                
                    if homology < 0.30:
                        filtered_info.append(mtb_info)
                        
                except:
                    print(mtb_seq)
                    print(other_seq)

        #Write the EC number, Protein name, gene name, and sequence to an excel file
        filtered_EC = [x['ECnumber'] for x in filtered_info]
        filtered_prot = [x['Protein Name'] for x in filtered_info]
        filtered_gene = [x['Gene'] for x in filtered_info]
        
        print("Number of EC numbers in Staphylococcus aureus not matched with sequence: " + str(other_miscount))
        print(len(filtered_info))
        file = open(r"C:\Github\pythonScripts\Mtb_inhibition\Mtb_inhibition\data\filter_noHumanStaph.tsv",'w')
        for i in range(0,len(filtered_EC)-1):
            file.write(filtered_EC[i] + '\t' + filtered_prot[i] + '\t' + filtered_gene[i] + '\n')
        file.close()

    ##################Now Mtb but not in Humans, Staphylococcus aureus, or Escherichia coli########################
    filtered_EC = [x['ECnumber'] for x in filtered_info]
    if not other_organism3 == None:
        other_list = list(ec_df[other_organism3].dropna())
        mtb_common = [x for x in filtered_EC if x in other_list]
        filtered_targets = [y for y in filtered_EC if not y in other_list]
        filtered_info = [fgetInfoFromUniprot(z,"Mycobacterium+tuberculosis") for z in filtered_targets]
        del filtered_targets
        
        #Pull Uniprot information for those proteins. If the % identity (score/length of alignment)is below X, send Mtb dict to results
        #What sort of alignment algorithm to use? Start with +1 for match, -1 for mismatch, 0 for gap.
        matrix = matlist.blosum62
        other_miscount = 0
        mtb_miscount = 0
        for prot in mtb_common:
            mtb_info = fgetInfoFromUniprot(prot,"Mycobacterium+tuberculosis")
            other_info = fgetInfoFromUniprot(prot,"Escherichia+coli")
            
            if not other_info['Sequence'] or other_info['Sequence'] == []: #if Uniprot search for homo sapiens returns nothing
                other_miscount = other_miscount+1
                #Check if also no sequence for Mtb
                if not mtb_info['Sequence'] or mtb_info['Sequence'] == []:
                    filtered_info.append(mtb_info)
                    
                    
                    
            elif not mtb_info['Sequence'] or mtb_info['Sequence'] == []: #if Uniprot search for mtb returns nothing
                #print("WARNING: Mtb enzyme " + mtb_info['ECnumber'] + " has no sequence in Uniprot")
                mtb_miscount = mtb_miscount+1
                #Write EC number to a txt file?
                    
            else:       
                mtb_seq = mtb_info['Sequence']

                #print(mtb_seq)
                other_seq = other_info['Sequence']

                #print(hs_seq)
                
                #Perform sequence alignment
                try:
                    tb_other_alignment = pairwise2.align.globaldx(mtb_seq, other_seq, matrix,score_only=True)
                    tb_tb_alignment = pairwise2.align.globaldx(mtb_seq, mtb_seq, matrix,score_only=True)
                    homology = tb_other_alignment/tb_tb_alignment
                    #print(homology)
                
                    if homology < 0.30:
                        filtered_info.append(mtb_info)
                        
                except:
                    print(mtb_seq)
                    print(other_seq)
                    
        print("Number of EC numbers in Escherichia coli not matched with sequence: " + str(other_miscount))
        print(len(filtered_info))   
        
        #Write the EC number, Protein name, gene name, and sequence to an excel file
        filtered_EC = [x['ECnumber'] for x in filtered_info]
        filtered_prot = [x['Protein Name'] for x in filtered_info]
        filtered_gene = [x['Gene'] for x in filtered_info]
        
        file = open(r"C:\Github\pythonScripts\Mtb_inhibition\Mtb_inhibition\data\filtered_noHumanStaphEcoli.tsv",'w')
        
        for i in range(0,len(filtered_EC)-1):
            file.write(filtered_EC[i] + '\t' + filtered_prot[i] + '\t' + filtered_gene[i] + '\n')
        file.close()
    return
    
def main():
    #organisms = ["Mycobacterium tuberculosis","Staphylococcus aureus","Escherichia coli","Homo sapiens"]
    #ec_lists = ecListbyOrg(organisms)
    filterByOrg('Homo sapiens','Staphylococcus aureus','Escherichia coli')
    
    
    #Need to change for Uniprot queries
    #organisms = ["Mycobacterium+tuberculosis","Staphylococcus+aureus","Escherichia+coli","Homo+sapiens"]
    
    #getInfoFromUniprot(elist,organisms)

    
if __name__ == "__main__":
    main()