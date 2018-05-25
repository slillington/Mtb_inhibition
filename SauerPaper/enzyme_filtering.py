
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
from SOAPpy import SOAPProxy ## for usage without WSDL file


def enzyme_filtering():
    
    #Compare these EC lists to those pulled from BRENDA
    orgs = ["Mycobacterium tuberculosis","Staphylococcus aureus","Escherichia coli","Candida albicans","Homo sapiens"]
    EC_by_orgs = []
    wsdl = "http://www.brenda-enzymes.org/soap/brenda_server.php"
    pw = hashlib.sha256("synbiorox").hexdigest()
    client = SOAPProxy(wsdl)
    for o in orgs:
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
        #print(ecNums)
        #print(len(ecNums))
    #print(EC_by_orgs)
    
    #Filter out EC numbers from Mtb list that are present in other lists
    #Iterate through EC numbers in Mtb EC list and check if in any of the other lists
    Mtb_ecs = EC_by_orgs[0]
    SA_ecs = EC_by_orgs[1]
    EC_ecs = EC_by_orgs[2]
    CA_ecs = EC_by_orgs[3]
    HS_ecs = EC_by_orgs[4]
    master = SA_ecs + EC_ecs + CA_ecs + HS_ecs
    results = []
    for e in Mtb_ecs:
        if e in master:
            continue
        elif unicode(e) in master:
            continue
        else:
            results.append(e)
    print(results)
            
        
    #Get enzyme names for each EC numbers
    for r in results:
        parameters = "stephenlillington2017@u.northwestern.edu,"+pw+",ecNumber*" + r + "#organism*Mycobacterium tuberculosis"
        resultString = client.getEnzymeNames(parameters)
        print(resultString)
    
    
    return


def main():
    enzyme_filtering()
    
    
if __name__ == "__main__":
    main()