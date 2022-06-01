# Downloads multiple molecule's substructures' CIDs and IsomericSMILES and creates files for each of them.
# Autor: Arthur A de Lacerda
# Autor: João Pedro Gonçalves Ribeiro
# Autor: Ramon Aragao
# Desenvolvimento: 20/05/2022
# Orientador: Edson Luiz Folador

import pandas as pd
from urllib.error import HTTPError
from urllib.request import urlopen
from urllib.parse import quote
from json import loads
from time import sleep
import os
from sys import argv


def get_result(url):
    """Accesses web api page"""
    try:
        connection = urlopen(url)
    except HTTPError:
        return None
    else:
        return connection.read().rstrip().decode('utf-8')


def get_listkey(smiles):
    """Gets list needed to get substructures"""
    result = get_result(
        f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/smiles/{smiles}/JSON")
    if not result:
        print("Could not get necessary listkey")
        return []

    listkey = loads(result)['Waiting']['ListKey']
    with open('listkey.txt', 'w') as f:
        f.write(listkey)


def listkey_to_substructures():
    """Gets substructures from listkey"""
    with open('listkey.txt', 'r') as list:
        listkey = list.read()
    result = get_result(
        f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{listkey}/property/IsomericSMILES/JSON")
    result = loads(result)['PropertyTable']['Properties']
    os.remove('listkey.txt')
    return result


def move_isosmiles():
    """Formats cid/smiles DataFrame"""
    global pubchem_df
    # Replaces cid as first column
    pubchem_df = pubchem_df[['IsomericSMILES'] + ['CID']]

    # Replaces main drug's cid for its name
    pubchem_df.loc[pubchem_df['IsomericSMILES'] == IsomericSMILES, 'CID'] = molecula

    # Replaces main drug as first row
    index = pubchem_df.index[pubchem_df['CID'] == molecula].to_list()
    idx = index + [i for i in range(len(pubchem_df)) if i != index[0]]
    pubchem_df = pubchem_df.iloc[idx]

    pubchem_df.reset_index(drop=True, inplace=True)


def create_files_pubchem():
    """Create pubchem files"""
    global pubchem_df
    # Creates file for main drug's smiles from pubchem
    with open(f'ligand/{molecula}-{molecula}.smi', 'w') as arqv:
        arqv.write(IsomericSMILES)

    # Creates files for each smiles from pubchem
    for i in range(1, len(pubchem_df)):
        with open(f'ligand/{molecula}-cid{pubchem_df["CID"][i]}.smi', 'w') as arqv:
            arqv.write(pubchem_df['IsomericSMILES'][i])


def create_files_zinc():
    """Create zinc files"""
    global zinc_df
    # Creates files for each smiles from zinc
    for i in range(len(zinc_df)):
        with open(f'ligand/{molecula}-{zinc_df["ZINC"][i]}.smi', 'w') as arqv:
            arqv.write(pubchem_df['IsomericSMILES'][i])


def _name(mol):
    """Gets initial information from the molecule's name"""
    global molecula, cid
    molecula = mol
    cid = get_result(
        f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{mol}/cids/TXT")
    if cid is not None:
        return True
    else:
        return False


def _cid(mol):
    """Gets initial information from a cid code"""
    global molecula, cid
    cid = mol
    molecula = get_result(
        f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{mol}/description/JSON").lower()
    if molecula is not None:
        molecula = loads(molecula)['InformationList']['Information'][0]['Title'].lower()
        return True
    else:
        return False


def _smiles(mol):
    """Gets initial information from a smiles code"""
    global molecula, cid
    smiles = get_result(
        f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/description/JSON?smiles={quote(mol)}")
    if smiles is not None:
        molecula = loads(smiles)['InformationList']['Information'][0]['Title'].lower()
        cid = loads(smiles)['InformationList']['Information'][0]['CID']
        return True
    else:
        return False


def _help():
    """Help function"""
    print(f"""
Usage: python3 {argv[0]} [-n] [-c] [-s] [-h] [molecule_descriptors]

Downloads multiple molecule's substructures' CIDs and IsomericSMILES and
creates files for each of them.

It takes each inputted molecule descriptor (one per molecule) and tries to
find information about said molecule on the Pubchem Compound Database and ZINC15.
If it finds them, it searches for similar ones by substructure and collects their
smiles and cid codes. Then it saves them  on the current path's ./ligand 
directory and saves each smiles retrieved a folder named after the input molecule,
each in a different ".smi" file.

Parameters:
  -n, --name                : Searches information on the base molecule by
                              name (Default).
  -c, --cid                 : Searches information on the base molecule by
                              cid.
  -s, --smiles              : Searches information on the base molecule by
                              smiles.
  -h, --help                : Prints this message.

At least one parameter must be passed.
Everything that is not one of these parameters is considered a molecule
descriptor.
If more than one of these parameters is passed, the latest one will take
effect. The exception is the help parameters, which will always print
this message and end the program regardless of position or other parameters.
""")
    exit(0)


# Checks if the required "ligand" directory exists in current directory.
if not os.path.isdir("./ligand"):
    print("Missing ./ligand directory")
    exit(1)

# Checks for conditions that would cause the program to print a help message and close
if ("-h" in argv) or ("--help" in argv) or (len(argv) < 2):
    _help()

key = "--name"  # Type of initial search to perform.
molecules = []  # List of molecules passed as parameters to search.

# Dictionary with flag parameters mapped to their respective search functions.
search_opts = {"-c": _cid, "--cid": _cid,
               "-n": _name, "--name": _name,
               "-s": _smiles, "--smiles": _smiles}

molecula = None  # Name of the initial molecule
cid = None  # CID of the initial molecule.

print()

# Iterates all parameters.
for param in range(1, len(argv)):
    # Checks if parameter is on Search Dictionary.
    if argv[param] in search_opts.keys():
        key = argv[param]
    else:
        molecules.append(argv[param])  # Appends molecules' descriptors to be processed to the list molecules

# Iterates all molecule descriptors in molecules playlist
for mol in molecules:
    print(f"{mol}:\n")

    # Executes search function and uses return to decide if it succeeded to retrieve the initial information.
    if search_opts[key](mol):
        IsomericSMILES = get_result(
            f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/TXT")
        print(f"\n{'-' * 30}\n")
        print(f'Drug:  {molecula}\n')
        print(f'CID:  {cid}\n')
        print(f'IsomericSMILES:  {IsomericSMILES}\n')
        print(f"{'-' * 30}\n\n")

        # Downloads substructures from PubChem
        get_listkey(IsomericSMILES)
        print('sleep 5sec')
        sleep(5)
        pubchem_df = pd.DataFrame(listkey_to_substructures())

        # Adjusts smiles on pubchem_df
        move_isosmiles()

        print(f'{len(pubchem_df.index)} substructures found on PubChem')

        # Creates file with all smiles and cid from pubchem
        pubchem_df.to_csv(f"./ligand/{molecula}-pubchem.txt", sep=' ', header=False, index=False)

        print("Creating pubchem .smi files...\n")

        create_files_pubchem()

        # Downloads substructures from ZINC15
        zinc = get_result(
            f'https://zinc15.docking.org/substances.smi?count=all&ecfp4_fp-tanimoto-30={quote(IsomericSMILES)}')

        # Creates file with all smiles and cid from ZINC15
        with open(f"./ligand/{molecula}-zinc.txt", 'w') as file:
            file.write(zinc)
        zinc_df = pd.read_csv(f'./ligand/{molecula}-zinc.txt', sep=' ', header=None)
        zinc_df.columns = ['IsomericSMILES', 'ZINC']
        print(f'{len(zinc_df.index)} substructures found on ZINC15')
        print("Creating zinc .smi files...\n")

        # Creates .smi files
        create_files_zinc()

        print('SCHLUSS!')
    else:
        print('Descriptor not in database')

    # Resets states of molecula and cid.
    molecula = None
    cid = None
    print(f"\n{'='*40}\n")
