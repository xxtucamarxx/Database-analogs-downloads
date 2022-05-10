import pandas as pd
from urllib.error import HTTPError
from urllib.request import urlopen
from json import loads
from time import sleep
import os
from sys import argv


def get_result(url):
    try:
        connection = urlopen(url)
    except HTTPError:
        return None
    else:
        return connection.read().rstrip().decode('utf-8')


def get_listkey(cid):
    result = get_result(f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/cid/{cid}/XML")
    if not result:
        return []
    listkey = ""
    for line in result.split("\n"):
        if line.find("ListKey") >= 0:
            listkey = line.split("ListKey>")[1][:-2]
    print(f'lista necessaria para as substruturas:\n{listkey}\n')
    with open('listkey.txt', 'w') as list:
        list.write(listkey)


def listkey_to_substructures():
    with open('listkey.txt', 'r') as list:
        listkey = list.read()
    result = get_result(
        f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{listkey}/property/IsomericSMILES/JSON")
    result = loads(result)['PropertyTable']['Properties']
    return result


def move_isosmiles(df):
    # Renomeia o codigo cid cujo isosmiles é igual ao do pubchem para o nome do farmaco
    df.loc[df['IsomericSMILES'] == IsomericSMILES, 'CID'] = molecula
    # Realoca o isosmiles do pubchem para primeiro da tabela
    index = df.index[df['CID'] == molecula].to_list()
    idx = index + [i for i in range(len(df)) if i != index[0]]
    df = df.iloc[idx]
    df.reset_index(drop=True, inplace=True)


def create_files():
    os.mkdir(f'ligand/{molecula}')
    os.remove('listkey.txt')
    with open(f'ligand/{molecula}/{molecula}-{molecula}.smi', 'w') as arqv:
        arqv.write(IsomericSMILES)
    for i in range(1, len(df)):
        with open(f'ligand/{molecula}/{molecula}-cid{df["CID"][i]}.smi', 'w') as arqv:
            arqv.write(df['IsomericSMILES'][i])


if len(argv) < 2:
    print("Você deve passar ao menos um parâmetro para este script.")
else:
    for mol in range(1, len(argv)):
        molecula = argv[mol]
        print(f"{molecula}:\n")

        cid = get_result(f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{molecula}/cids/TXT")

        if cid is not None:
            IsomericSMILES = get_result(
                f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/TXT")
            print(f'Farmaco: \n{cid}')
            print(f'Cid da farmaco:\n{cid}')
            print(f'Isosmiles da farmaco:\n{IsomericSMILES}')

            get_listkey(cid)
            sleep(5)
            df = pd.DataFrame(listkey_to_substructures())

            print("Execuntando script...\n")

            # Ajusta os isomiles
            move_isosmiles(df)

            # Cria arquivos .smi
            create_files()

            print('SCHLUSS!')
        else:
            print('Farmaco nao presente no banco de dados')

        print(f"\n{'-'*40}\n")
