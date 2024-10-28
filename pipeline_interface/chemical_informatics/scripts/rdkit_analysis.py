# Script for RDKit analysis and descriptor calculations
from rdkit import Chem
from rdkit.Chem import Descriptors
import h5py
import csv

def calculate_descriptors(input_file, output_hdf5):
    weights = []
    with open(input_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            mol = Chem.MolFromSmiles(row[0])
            if mol:
                weights.append(Descriptors.MolWt(mol))

    with h5py.File(output_hdf5, 'w') as f:
        f.create_dataset('mol_weights', data=weights)

if __name__ == "__main__":
    calculate_descriptors('../data/chembl_data.csv', '../data/chemical_data.h5')

