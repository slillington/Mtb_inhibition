#M tuberculosis FBA script

from __future__ import print_function
import cobra
import libsbml
from cobra.util.solver import linear_reaction_coefficients
import pandas as pd
from cobra.flux_analysis import single_reaction_deletion
'''Workflow for this script:
1. Pass file containing sorted compound IDs, InChiKeys, and name pairs with Tanimoto coefficients
2. For each metabolite in a pair (maybe exclude cofactors for now?) find the reactions in which the metabolite is a reactant
3. Can do this by model.metabolite.get_by_id(<cpdXXXXX_c0>).reactions
4. Look at the resulting reactions and pick the ones for which the coefficient for metabolite X is negative (aka is a reactant)
5. Set reaction flux to 0 to simulate knockout and store estimated growth rate (biomass production rate)
6. Identify targets that result in no or significantly reduced growth rate

###BUG - iNJ661 model when processed by KBASE has something wrong with it. It's most likely that a reaction is missing.

THIS IS RUNNING FOR E COLI RIGHT NOW!!!!!!!!!!


'''
def fba_on_sim_mets():
	'''
	#BiGG Model
	model = cobra.io.load_json_model('C:\Github\pythonScripts\Mtb_inhibition\Mtb_inhibition\iNJ661.json')
	model.reactions.get_by_id("EX_glc__D_e").upper_bound = -10
	model.reactions.get_by_id("EX_o2_e").upper_bound = -18
	model.optimize()
	print(model.summary())
	
	
	'''
	#KBASE Model
	model_KBASE = cobra.io.read_sbml_model('C:\GitHub\pythonScripts\Mtb_inhibition\model_objects\model_objects\iAF1260b.xml_model.sbml')
	model_BiGG = cobra.io.load_json_model('C:\GitHub\pythonScripts\Mtb_inhibition\Mtb_inhibition\iNJ661.json')
	#Add needed reactions
	'''
	#Set flux for glucose uptake to be 10
	#model.reactions.get_by_id("EX_cpd00027_e0").lower_bound = -1000
	model_KBASE.reactions.get_by_id("EX_cpd00027_e0").upper_bound = -5
	#model.reactions.get_by_id("GLCabc_c0").lower_bound = 10
	model_KBASE.reactions.get_by_id("EX_cpd00007_e0").upper_bound = -7
	#model.reactions.get_by_id("ATPS4r_c0").lower_bound = 0
	#model.objective = "BIOMASS_Mtb_9_60atp_c0"
	model_KBASE.objective = "bio1"
	#print(model.reactions.get_by_id("bio1").bounds)
	#print(model.reactions.get_by_id("BIOMASS_Mtb_9_60atp_c0").bounds)
	#model.reactions.get_by_id("BIOMASS_Mtb_9_60atp_c0").upper_bound = 0
	#print(model.reactions)
	#model.optimize()
	#print(model.summary())
	
	#Set flux for glucose uptake for BiGG model
	model_BiGG.reactions.get_by_id("EX_glc__D_e").upper_bound = -5
	model_BiGG.reactions.get_by_id("EX_o2_e").upper_bound = -7
	
	#Read in list of compound IDs
	xls = pd.ExcelFile(r'C:\GitHub\pythonScripts\Mtb_inhibition\Mtb_inhibition\Output_EcoliiAF1260b_AllMetabolites_5_3.xlsx')
	sheetX = xls.parse(0) #0 is the sheet number
	mets1 = sheetX['Met 1 cID']
	mets1 = list(set(mets1)) #Remove duplicates b/c we only need to test the knockouts once
	mets2 = sheetX['Met 2 cID']
	mets2 = list(set(mets2)) #Remove duplicates b/c we only need to test the knockouts once
	
	
	#Iterate through the list of pairs and see which ones give growth rates of 0
	#Remove any cofactors from the list (NAD, NADP, NADH, NADPH, AMP, ADP, ATP)
	rxns_list = []
	for m in mets1:
		if str(m) in ['nan','Nan']:
			continue
		else:
			a1 = model_KBASE.metabolites.get_by_id(m)
			for r in a1.reactions:
				temp = r.metabolites #Returns a dict with the metabolites in each reaction involving the specified metabolite and the coefficient.
				#print(temp)
				#print(temp[a1])
				if temp[a1] < 0:
					rxns_list.append(r)
	for m2 in mets2:
		if str(m2) in ['nan','Nan']:
			continue
		else:
			a2 = model_KBASE.metabolites.get_by_id(m2)
			for r2 in a2.reactions:
				temp = r2.metabolites #Returns a dict with the metabolites in each reaction involving the specified metabolite and the coefficient.
			
				if temp[a2] < 0:
					rxns_list.append(r2)
	'''			

	#Simulate single reaction knockouts for list of rxns_list
	#a = model.optimize()
	#print(a.fluxes)
	#print(model.summary())
	#table = single_reaction_deletion(model_KBASE,rxns_list)
	table2 = single_reaction_deletion(model_BiGG,model_BiGG.reactions)
	FBAoutput = pd.ExcelWriter('C:\Github\pythonScripts\Mtb_inhibition\Mtb_inhibition\FBAoutput_Mtb_all.xlsx')
	#table.to_excel(FBAoutput,'Sheet1',startcol=0)
	table2.to_excel(FBAoutput,'Sheet1',startcol=0)
	FBAoutput.save()
	
	return

def main():
	fba_on_sim_mets()
	
	
if __name__ == "__main__":
	main()
