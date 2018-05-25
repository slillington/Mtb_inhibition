#My settings file
#Contains necessary import statements, defined functions, etc.

import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import imp
import json
import inspect
import urllib 
from contextlib import closing

SCRIPT_DIR = os.getcwd()
BASE_DIR = os.path.join(*os.path.split(SCRIPT_DIR)[0:-1])

DATA_DIR = os.path.join(SCRIPT_DIR, 'data')
CACHE_DIR = os.path.join(SCRIPT_DIR, 'cache')
RESULT_DIR = os.path.join(SCRIPT_DIR, 'res')
CHEBI2INCHI_URL = os.path.join(DATA_DIR,'chebiId_inchi.tsv')
KEGG2CHEBI_FNAME = os.path.join(CACHE_DIR, 'kegg2chebi.csv')
BIGG_METABOLITE_FNAME = os.path.join(DATA_DIR, 'bigg_models_metabolites.txt')
BIGG_REACTION_FNAME = os.path.join(DATA_DIR, 'bigg_models_reactions.txt')
MTB_JSON_FNAME = os.path.join(DATA_DIR, 'iNJ661.json')
MTB_SBML_FNAME = os.path.join(DATA_DIR, 'iNJ661.xml')
ECOLI_JSON_FNAME = os.path.join(DATA_DIR,'iJO1366.json')
#ECOLI_JSON_FNAME = os.path.join(DATA_DIR,'iAF1260b.json')
STAPHA_JSON_FNAME = os.path.join(DATA_DIR,'iSB619.json')
HUMAN_JSON_FNAME = os.path.join(DATA_DIR,'iAT_PLT_636.json')
#HUMAN_JSON_FNAME = os.path.join(DATA_DIR,'RECON1.json')

def get_data_df(fname):
    return pd.DataFrame.from_csv(os.path.join(DATA_DIR, fname + '.csv'), header=0, index_col=None)

def write_cache(fname, df):
    df.to_csv(os.path.join(CACHE_DIR, fname + '.csv'))

def read_cache(fname):
    return pd.read_csv(os.path.join(CACHE_DIR, fname + '.csv'), index_col=0,
                       encoding='latin-1')

def get_reaction_table_from_xls():
    with open(ECOLI_XLS_FNAME) as fp:
        return pd.read_excel(fp, sheetname=2, header=0)

def plotdiag(lw=2, ax=None):
    if ax is None:
        ax = plt
    x1, x2, y1, y2 = ax.axis()
    minplot = np.min([x1, y1])
    maxplot = np.max([x2, y2])
    ax.plot([minplot,    maxplot],    [minplot,    maxplot],    'k-',  lw=lw)
    ax.plot([minplot,    maxplot/10], [minplot*10, maxplot],    'k--', lw=lw)
    ax.plot([minplot*10, maxplot],    [minplot,    maxplot/10], 'k--', lw=lw)

try:
    imp.find_module('cobra')
    cobrafound = True
except ImportError:
    cobrafound = False
    sys.stderr.write("WARNING: please install cobrapy to have full functionality")
if cobrafound:
    def get_ecoli_sbml():
        return cobra.io.read_sbml_model(ECOLI_SBML_FNAME)


def get_org_json(ORG_JSON_FNAME):
    with open(ORG_JSON_FNAME) as fp:
        model = json.load(fp)

    sparse = []
    for reaction in model['reactions']:
        for met, coeff in reaction['metabolites'].items():
            sparse.append([reaction['id'].lower(), met.lower(), coeff])

    sparse = pd.DataFrame(sparse, columns=['bigg.reaction',
                                           'bigg.metabolite', 'stoichiometry'])
    sparse.drop_duplicates(subset=['bigg.reaction','bigg.metabolite'],keep=False)
    #print(sparse)
    S = sparse.pivot(index='bigg.metabolite', columns='bigg.reaction', values='stoichiometry')
    S.fillna(0, inplace=True)
    return model, S


def get_chebi_inchi_df():
    #http = urllib3.PoolManager()
    with open(CHEBI2INCHI_URL,'r') as r:
        chebi_inchi_df = pd.read_csv(r, sep='\t')
    chebi_inchi_df['chebiID'] = chebi_inchi_df['CHEBI_ID'].apply(lambda c: 'CHEBI:%d' % c)
    chebi_inchi_df.rename(columns={'InChI': 'inchi'}, inplace=True)
    chebi_inchi_df.set_index('CHEBI_ID', inplace=True)
    return chebi_inchi_df


def savefig(fig, name, dpi=300):
    fig.savefig(os.path.join(RESULT_DIR, name + '.svg'))
    fig.savefig(os.path.join(RESULT_DIR, name + '.pdf'))
    fig.savefig(os.path.join(RESULT_DIR, name + '.png'), dpi=dpi)