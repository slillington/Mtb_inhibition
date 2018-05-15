#Performs everything I wrote in the Mtb_inhibition_calcs_nucs.py but uses the unprocessed BiGG genome scale model and interfaces with ChEBI and BRENDA to get better information.
"""
This script should accomplish the following:
1. Generate a list of InChi strings for all metabolites in a genome scale model
2. Translate those InChi strings into Mol objects in rdkit
3. Compute pairwise Tanimoto similarities between all metabolites
4. Compile BRENDA information on kinetic constants, inhibitors, and inhibitory constants from BRENDA
5. Identify lethal knockouts for Mtb under aerobic growth conditions
6. Compile a table of already observed inhibitors with their Ki's (if available) for the substrates of these lethal reactions

"""
import settings
from bigg import BiGG

bigg = BiGG()
mtb_model = BiGG()

BIGG_METABOLITESALL_FNAME = os.path.join(DATA_DIR, 'bigg_models_metabolites.txt')
#Generate a list of InChi strings for all metabolites in a BiGG model
def getInchis():

    #Get ChEBI IDs for the universal metabolite model in BiGG
    chebiIDs = bigg._get_metabolite_df(BIGG_METABOLITESALL_FNAME)
    chebi_inchis = settings.get_chebi_inchi_df() #a dataframe with ID and InChi string
    #print(chebi_inchis.get(['chebiID']))
    
    #Get ChEBI IDs for model of interest
    mdl = '
    
    #Extract InChis for ChEBI IDs that are in both dataframes
    model_inchis = []
    



def main():
    getInchis()
    
    
if __name__ == "__main__":
    main()