#M tuberculosis FBA script

from __future__ import print_function
import cobra
import libsbml
from cobra.util.solver import linear_reaction_coefficients
import pandas as pd
'''Workflow for this script:
1. Pass file containing sorted compound IDs, InChiKeys, and name pairs with Tanimoto coefficients
2. For each metabolite in a pair (maybe exclude cofactors for now?) find the reactions in which the metabolite is a reactant
3. Can do this by model.metabolite.get_by_id(<cpdXXXXX_c0>).reactions
4. Look at the resulting reactions and pick the ones for which the coefficient for metabolite X is negative (aka is a reactant)
5. Set reaction flux to 0 to simulate knockout and store estimated growth rate (biomass production rate)
6. Identify targets that result in no or significantly reduced growth rate


'''

model = cobra.io.read_sbml_model('C:\GitHub\pythonScripts\Mtb_inhibition\model_objects\model_objects\iNJ661.xml_model.sbml')

#Read in list of compound IDs
xls = pd.ExcelFile(r'C:\GitHub\pythonScripts\Mtb_inhibition\Mtb_inhibition\toy_data_output_test.xlsx')
sheetX = xls.parse(0) #2 is the sheet number
mets1 = sheetX['Met 1 cID']
mets2 = sheetX['Met 2 cID']
	
#Iterate through the list of pairs and see which ones give growth rates of 0
idx = 0
while idx < len(mets1):
	a1 = model.metabolites.get_by_id(mets1[idx])
	a2 = model.metabolites.get_by_id(mets2[idx])
	
	#Look through Met 1 reactions
	for rxn in a1.reactions:
		temp = rxn.metabolites #Returns a dict with the metabolites in each reaction involving the specified metabolite and the coefficient.
		if temp[a1] < 0:
			#simulate FBA and get store growth rate
		
			with model:
				rxn.knock_out()
				print('Growth rate with ',rxn,' knocked out is ',model.optimize())
				
	#Look through Met 2 reactions
	for rxn in a2.reactions:
		temp = rxn.metabolites #Returns a dict with the metabolites in each reaction involving the specified metabolite and the coefficient.
		if temp[a2] < 0:
			#simulate FBA and get store growth rate
		
			with model:
				rxn.knock_out()
				print('Growth rate with ',rxn,' knocked out is ',model.optimize())
	
	idx = idx+1
