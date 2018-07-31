
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
from pandas import ExcelWriter

def getInfoFromBrenda():
    
    #Compare these EC lists to those pulled from BRENDA
    ec_df = pd.read_excel(r"C:\GitHub\pythonScripts\Mtb_inhibition\Mtb_inhibition\data\ec_numbers_lists.xlsx")
    ec_list = list(ec_df["Escherichia coli"].dropna())
    wsdl = "http://www.brenda-enzymes.org/soap/brenda_server.php"
    pw = hashlib.sha256("synbiorox").hexdigest()
    client = SOAPProxy(wsdl)
    
    #Information to extract
    ecNums = []
    inhibitors = []
    commentary = []
    #ki = []
    #Loop through EC list
    #ec_list = ec_list[0:100]
    for e in ec_list:
        parameters = "stephenlillington2017@u.northwestern.edu,"+pw+",ecNumber*"+e+"#organism*Escherichia coli"
        try:
            resultString = client.getInhibitors(parameters)

            #print(resultString)
            if resultString:
                #Separate by the various entries
                entries = resultString.split('!')
            

                for e in entries:
                    i = e.find("ecNumber*")
                    h = e.find("#",i+15)
                    ecNums.append(e[i+9:h])
                    
                    i = e.find("inhibitors*")
                    h = e.find("#",i+11)
                    inhibitors.append(e[i+11:h])
                    
                    i = e.find("commentary")
                    h = e.find("#",i+11)
                    commentary.append(e[i+11:h])
                    
        except:
            print("BRENDA access failed for EC number " + e)
            
            
                
    #Now write all of the lists to a dataframe
    results = pd.DataFrame(data={'EC Number':ecNums,'Inhibitors':inhibitors,'Commentary':commentary},columns=['EC Number','Inhibitors','Commentary'])
    #print(results)
    file = ExcelWriter("brenda_info_results.xlsx")
    results.to_excel(file)
    file.save()

    return


def main():
    getInfoFromBrenda()
    
    
if __name__ == "__main__":
    main()