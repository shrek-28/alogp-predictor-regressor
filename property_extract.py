import os 
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

sdf_file = "COCONUT_2022_01_2D_DS_CNP0000001-CNP0080000.sdf"

molecule_data = []

suppl = Chem.SDMolSupplier(sdf_file)

count = 0

for mol in suppl: 
    if mol is None:
        continue

    molecule_info = {
            "Coconut_ID": mol.GetProp("coconut_id") if mol.HasProp('coconut_id') else None,
            "InChI": Chem.MolToInchi(mol, options='/FixedH') if mol else None,
            "InChIKey": Chem.MolToInchiKey(mol) if mol else None,
            "SMILES": Chem.MolToSmiles(mol, isomericSmiles=False),
            "Sugar_Free_SMILES": mol.GetProp("sugar_free_smiles") if mol.HasProp("sugar_free_smiles") else None,
            "Molecular_Formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "Molecular_Weight": Descriptors.MolWt(mol),
            "Citation_DOI": mol.GetProp("citationDOI") if mol.HasProp("citationDOI") else None,
            "Taxa_Text": mol.GetProp("textTaxa") if mol.HasProp("textTaxa") else None,
            "Synonyms": mol.GetProp("synonyms") if mol.HasProp("synonyms") else None,
            "NPL_Score": float(mol.GetProp("NPL_score")) if mol.HasProp("NPL_score") else None,
            "Number_of_Carbons": int(mol.GetProp("number_of_carbons")) if mol.HasProp("number_of_carbons") else None,
            "Number_of_Nitrogens": int(mol.GetProp("number_of_nitrogens")) if mol.HasProp("number_of_nitrogens") else None,
            "Number_of_Oxygens": int(mol.GetProp("number_of_oxygens")) if mol.HasProp("number_of_oxygens") else None,
            "Total_Atom_Number": int(mol.GetProp("total_atom_number")) if mol.HasProp("total_atom_number") else None,
            "Bond_Count": int(mol.GetProp("bond_count")) if mol.HasProp("bond_count") else None,
            "Found_in_Databases": mol.GetProp("found_in_databases") if mol.HasProp("found_in_databases") else None,
            "Murko_Framework": mol.GetProp("murko_framework") if mol.HasProp("murko_framework") else None,
            "ALogP": Descriptors.MolLogP(mol),
            "APolarSurfaceArea": Descriptors.TPSA(mol),
            "Topo_PSA": Descriptors.TPSA(mol),
            "Clean_Energy": mol.GetProp("Clean Energy") if mol.HasProp("Clean Energy") else None
    }

    molecule_data.append(molecule_info)
    count+=1 

    if count % 1000 ==0: 
        print(f"{count} molecules processed")

df = pd.DataFrame(molecule_data)

output_csv = "molecule_details_fin.csv"
df.to_csv(output_csv, index=False)

print("molecule details saved.")