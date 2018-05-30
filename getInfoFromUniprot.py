#Testing Uniprot programmatic access and sequence alignment

import requests
#ec_nums = ["2.7.1.175","1.14.14.12","2.4.99.16"]


def fgetInfoFromUniprot(ec_number,organism):
#Returns a list of dictionaries with protein ID, Organism, Protein name, Gene name, and Sequence

    url = "https://www.uniprot.org/uniprot/?query=ec%3A" +ec_number+"+AND+organism%3A%22" +organism+"%22&sort=score&limit=1&columns=id,protein names,genes,sequence&format=fasta"
    r = requests.get(url)
    rtext = r.text
    
    #parse out text and put in dictionary format
    #keys: ID, Organism, Protein Name, Gene name, Sequence
    id_start = rtext.find('|')
    id_end = rtext.find('|',id_start+1)
    id = rtext[id_start+1:id_end]
    #print(id)
    
    name_end = rtext.find('OS=',id_end)
    name = rtext[id_end+1:name_end-1]
    #print(name)
    
    org_end = rtext.find('OX=',name_end)
    org = rtext[name_end+3:org_end-1]
    #print(org)
    
    gene_start = rtext.find('GN=',org_end)
    gene_end = rtext.find('PE=',gene_start)
    gene = rtext[gene_start+3:gene_end-1]
    #print(gene)
    
    sequence_start = rtext.find('\n',gene_end)
    sequence = rtext[sequence_start+1:-1]
    #print(sequence)
    sequence2 = sequence.replace(u"\n",u"")
    #print(sequence2)
    temp_dict = {'ECnumber':ec_number,'ID':id,'Organism':org,'Protein Name':name,'Gene':gene,'Sequence':sequence2}
        
    return temp_dict
        

