# Doing all tanimoto analysis

def fingerprint_and_tanimoto(mol_list):
	fps_rdk = [FingerprintMols.FingerprintMol(x,fingerprinter=Chem.RDKFingerprint) for x in mol_list]
	#fps_m = 
	#print(len(fps_rdk))
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



if __name__ == "__main__":
	main()