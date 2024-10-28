# Script for fetching ChEMBL data
from chembl_webresource_client.new_client import new_client
import csv

def fetch_chembl_data(target_name, output_file):
    target = new_client.target
    bioactivity = new_client.activity
    target_query = target.search(target_name)
    target_id = target_query[0]['target_chembl_id']
    bioactivities = bioactivity.filter(target_chembl_id=target_id, standard_type="IC50")

    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["SMILES", "IC50"])
        for entry in bioactivities:
            if entry['canonical_smiles']:
                writer.writerow([entry['canonical_smiles'], entry['standard_value']])

if __name__ == "__main__":
    fetch_chembl_data('BRAF', '../data/chembl_data.csv')
