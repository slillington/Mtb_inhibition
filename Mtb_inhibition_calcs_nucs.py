''' Mass spec library curation using information derived from KBase

BUGS:
1. (stereo)Isomer SMILES keys are not differentiable (example maleylpyruvate and 3-fumarylpyruvate)
*KEGG IDs are not repeated when SMILES are repeated and InChiKeys are NOT repeated

Repeated SMILES (different InChiKeys):
maleylpyruvate, 3-fumarylpyruvate
2-keto-3-deoxygluconate_c0, 2-Dehydro-3-deoxy-D-galactonate_c0
R-acetoin, S-acetoin
D-mannitol-1-phosphate_c0, Galactitol 1-phosphate_c0
D-fructose-6-phosphate_c0, D-Tagatose 6-phosphate_c0, L-Tagatose 6-phosphate_c0
D-fructose-1,6-bisphosphate_c0, D-Tagatose 1,6-biphosphate_c0
L-Fuculose_c0, L-Rhamnulose_c0
L-Fuculose-1-phosphate_c0, L-Rhamnulose-1-phosphate_c0
Glucose-1-phosphate_c0, D-Mannose1-phosphate_c0, alpha-D-Hexose 1-phosphate_c0, beta-D-Glucose 1-phosphate_c0, D-Galactose 1-phosphate_c0
UDP-glucosamine_c0, UDP-galactosamine_c0
gluconate_c0, D-Altronate_c0, L-Idonate_c0, D-Mannonate_c0
LL-2,6-Diaminopimelate_c0, meso-2,6-Diaminopimelate_c0
D-Mucic acid_c0, D-Glucarate_c0
5-Dehydro-4-deoxy-D-glucarate_c0, 2-Dehydro-3-deoxy-D-glucarate_c0
Retinol_c0, 11-cis-Retinol_c0
Retinol palmitate_c0, 11-cis-Retinyl palmitate_c0
CDPparatose_c0, CDPtyvelose_c0
D-methylmalonyl-CoA_c0, L-methylmalonyl-CoA_c0
(R)-Allantoin_c0, Allantoin_c0
UDP-N-acetylglucosamine_c0, UDP-N-acetyl-D-galactosamine_c0, UDP-N-acetyl-D-mannosamine_c0, 
cellobiose, lactose
L-ribulose-5-phosphate_c0, D-Xylulose5-phosphate_c0, D-Ribulose5-phosphate_c0, L-Xylulose 5-phosphate_c0
D-Ribulose_c0, L-Ribulose_c0, D-Lyxulose_c0, L-Lyxulose_c0
(R)-1,2-Propanediol_c0, 1,2-Propanediol_c0
D-Lactaldehyde_c0, L-Lactaldehyde_c0
cis-Aconitate_c0, trans-Aconitate_c0
GDP-mannose, GDP-glucose
(R)-(+)-Citronellal_c0, (S)-(-)-Citronellal_c0


'''


import pandas as pd
import xlrd
import requests
import json
import pubchempy as pcp
from io import BytesIO
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from operator import itemgetter


		
def get_model_data(excel_file_path):

#get_model_data take the list of inchikeys and smiles from an excel file that contains the .tsv download from KBASE. This .tsv download has information
#on the metabolites in a genome-scale metabolic model.
#type input: string
#input: file path to excel file containing model metabolite information
#type output inchiKeys: list of InChiKeys (type string) in the model
#type output smiles: list of SMILES (type string) in the model   

    #xls = pd.ExcelFile('escherichia_coliMG1655_mets.xlsx')
	xls = pd.ExcelFile(excel_file_path)
	sheetX = xls.parse(1) #1 is the sheet number (c0 only)

	inchiKeys = sheetX['inchikey'] #Original length = 1275 InChiKeys for E coli model
	smiles = sheetX['smiles']

	inchiKeys = list(set(inchiKeys)) #After removing duplicates, 908 InChiKeys for E coli model. Assuming this removes <compound>_c0 and <compound>_e0 duplicates
	inchiKeys = [x for x in inchiKeys if not str(x) in ['NaN','nan']]
	smiles = [y for y in smiles if not str(y) in ['NaN','nan']]
	inchiKeys.sort()
	smiles = list(set(smiles))
	smiles.sort()
	
	return inchiKeys, smiles




	
def InChiKeyToInChi(keys):
#InChiKeyToInChi() take the set of InChiKeys passed as an argument and uses the ChemSpider InChiKey --> InChi functionality to pair a corresponding
#InChi to a given InChiKey. This take a while to run - on the order of 10 minutes for a list of ~900 InChiKeys.
#type input: list[string]
#input: list of InChiKey strings
#type output: list[dict]
#output: list of dictionaries with InChiKey and InChi value pairs for metabolites in the model
	

#For each InChiKey, search the ChemSpider database for a corresponding InChi string
	url = 'https://www.chemspider.com/InChI.asmx?op=InChIKeyToInChI'
	headers = {'content-type': 'text/xml'}
	dict_list = []
	
	#Report progress and get InChis from ChemSpider
	print('Getting InChis from ChemSpider.')
	for inChiKey in keys:
		body = """<?xml version="1.0" encoding="utf-8"?>
<soap:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
  <soap:Body>
    <InChIKeyToInChI xmlns="http://www.chemspider.com/">
      <inchi_key>""" +str(inChiKey)+"""</inchi_key>
    </InChIKeyToInChI>
  </soap:Body>
</soap:Envelope>"""
    #print(body)
		response = requests.post(url,data=body,headers=headers)
		txt = BytesIO(response.content)
		txt = txt.getvalue().decode('utf-8')
		idx = txt.find('InChI=')
		inchi = ''
		while True:
			if txt[idx] == '<':
				break
			else:
				inchi = inchi + txt[idx]
				idx = idx+1
		#if InChi was not found, value for "InChi" key is '<'
		temp = dict([("InChiKey",inChiKey),("InChi",inchi)])
		dict_list.append(temp)
		
	#Check the list for InChiKeys that were not paired with a corresponding InChi
	print('Checking ChemSpider results for InChiKeys with no InChi pair')
	trouble_keys = []
	inChis_inChiKeys_file = open(r'C:\Github\pythonScripts\Mtb_inhibition\inchis_inchikeys_iNJ661.txt','w+')
	#Tweak this one time to write to a file that I can use in the future
	for d in dict_list:
		if d["InChi"] == '>':
			url = 'http://cts.fiehnlab.ucdavis.edu/rest/convert/InChiKey/InChi%20Code/' + d["InChiKey"]
			r = requests.get(url)
			r = r.json()
			r = r[0]
			#r = r['result']
			if r['result']:
				d["InChi"] = r['result'][0]
				
				inChis_inChiKeys_file.write(str(d["InChiKey"]) + "	" + str(d["InChi"]) + "\n")
				print('Write success')
			else:
				try:
					print('No InChi from CTS. Trying PubChem')
					result = pcp.get_compounds(d["InChiKey"],'inchikey')
					d["InChi"] = (result[0].inchi)
					inChis_inChiKeys_file.write(str(d["InChiKey"]) + '	' + str(d["InChi"]) + '\n')
					print('Write success')
					
				except:
					trouble_keys.append(d["InChiKey"])
					print(d["InChiKey"])
					d["InChi"] = 'InChI=1S/Hg/q+2'
					inChis_inChiKeys_file.write(str(d["InChiKey"]) + '	' + str(d["InChi"]) + '\n')
					print('Write success')
		else:
			#print(type(d["InChiKey"]))
			#print(d["InChiKey"])
			#print(type(d["InChi"]))
			#print(d["InChi"])
			inChis_inChiKeys_file.write(str(d["InChiKey"]) + '	' + str(d["InChi"]) + '\n')
			print('Write success')
					
				
				
	#print (trouble_keys)
	#print(len(trouble_keys))
	trouble_file = open(r'C:\Github\pythonScripts\Mtb_inhibition\trouble_keys_file_toy.txt','w+')
	
	for k in trouble_keys:
		trouble_file.write(str(k)+'\n')
		
	return dict_list, trouble_keys

	
def inchiToMol(excel_file,sheet):
	mol_list = []
	data = pd.ExcelFile(excel_file)
	sheetX = data.parse(sheet) #is the sheet number (c0 only)
	#inchiKeys = sheetX['inchikey'] #Original length = 1275 InChiKeys for E coli model
	inchi = sheetX['inchi']
	#print(inchi)
	for s in inchi:	
		if not s == '>':
			mol = Chem.MolFromInchi(s,sanitize=True,removeHs=True,logLevel=None,treatWarningAsError=False)
			mol_list.append(mol)
		#else:
		#	mol_list.append(Chem.MolFromInchi('InChI=1S/Hg/q+2')) #InChiKeys with missing InChis are for now Mercury ion	
	return mol_list

	
def smilesToMol(smiles):
	mol_list = []
	for s in smiles:
		mol = Chem.MolFromSmiles(s)
		mol_list.append(mol)
	return mol_list

	
def fingerprint_and_tanimoto(mol_list):
	fps_rdk = [FingerprintMols.FingerprintMol(x,fingerprinter=Chem.RDKFingerprint) for x in mol_list]
	#fps_m = 
	print(len(fps_rdk))
	#Compute pairwise Tanimoto similarity for each pair of fingerprints
	#tanimotos is a list of lists
	
	#Loop through list in this way to not get repeat values or compare fingerprints to themselves
	tanimotos = []
	for i in range(len(fps_rdk)-1):
		for j in range(i+1,len(fps_rdk)):
			#temp_dict = dict([("Fingerprints",[fps[i],fps[j]]),("Tanimoto_coeff",DataStructs.FingerprintSimilarity(fps[i],fps[j]))])
			temp_dict = dict([("InChiKeys",[Chem.InchiToInchiKey(Chem.MolToInchi(mol_list[i])),Chem.InchiToInchiKey(Chem.MolToInchi(mol_list[j]))]),("Tanimoto_coeff",DataStructs.FingerprintSimilarity(fps_rdk[i],fps_rdk[j]))])
			tanimotos.append(temp_dict)
	
	#tanimotos is a list of dicts. Now want to search the tanimoto coeff key for highest values
	tanimotos.sort(key=lambda x:x['Tanimoto_coeff'])
	
	#Write the dictionary contents to a text file
	out = open(r'C:\Github\pythonScripts\Mtb_inhibition\Mtb_inhibition\tanimoto_output_Ecoli_all.txt','w')
	for t in tanimotos:
		key_pair = t["InChiKeys"]
		out.write(key_pair[0] + '	' + key_pair[1] + '	' + str(t["Tanimoto_coeff"]) + '\n')



		
def main():
	#inChiKeys,smiles = get_model_data(r'C:\GitHub\pythonScripts\Mtb_inhibition\Mtb_inhibition\mycobacterium_tuberculosis_mets.xlsx')


	source = input("Create Mol objects from InChi or SMILES?:   ")
	if source == 'InChi':
		#met_dicts, troubles = InChiKeyToInChi(inChiKeys)
		
		#Convert InChi to Mol
		mol_list = inchiToMol('C:\Github\pythonScripts\Mtb_inhibition\Mtb_inhibition\Mtb_iNJ661_inchis_and_inchikeys.xlsx',0)
	elif source == 'SMILES': #WARNING WILL NOT WORK RIGHT NOW
		mol_list = smilesToMol(smiles)
			
	else:
		raise IOError(source + ' is Not a valid input')

	#print(len(mol_list))
	#Convert mol objects to bit vectors and calculate tanimoto coefficients
	fingerprint_and_tanimoto(mol_list)	
	
		
		
if __name__ == "__main__":
	main()


