#Testing Uniprot programmatic access and sequence alignment

import requests
ec_nums = ["2.7.1.175"]

for e in ec_nums:
    url = "https://www.uniprot.org/uniprot/?query=ec%3A" +e+"+AND+organism%3A%22Mycobacterium+tuberculosis%22&sort=score&limit=1&columns=id,protein names,genes,sequence&format=fasta"
    #url = "https://www.uniprot.org/uniprot/?query=ec%3A2.7.1.175+AND+organism%3A%22Mycobacterium+tuberculosis%22&sort=score&columns=id,reviewed,protein names,genes,sequence&format=fasta"
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
    temp_dict = {'ID':id,'Organism':org,'Protein Name':name,'Gene':gene,'Sequence':sequence}
    print(temp_dict)


