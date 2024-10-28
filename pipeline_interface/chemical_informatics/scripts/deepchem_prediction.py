# Script for DeepChem model training and prediction
import deepchem as dc
import numpy as np
import h5py

def load_data_from_hdf5(hdf5_file):
    with h5py.File(hdf5_file, 'r') as f:
        weights = np.array(f['mol_weights'])
    return weights

def train_predict_model(data):
    dataset = dc.data.NumpyDataset(X=data.reshape(-1, 1))
    model = dc.models.GraphConvModel(n_tasks=1, mode='regression')
    model.fit(dataset)
    predictions = model.predict(dataset)
    return predictions

if __name__ == "__main__":
    data = load_data_from_hdf5('../data/chemical_data.h5')
    predictions = train_predict_model(data)
    print("Predictions:", predictions)
