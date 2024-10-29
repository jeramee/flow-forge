# flow-forge

Container for Informatics Pipelines

## Informatics Project: Bioinformatics, Cheminformatics, and Machine Learning Integration

This project integrates bioinformatics, cheminformatics, and machine learning workflows, utilizing tools like `cbio_ml_pl`, HDF5, RDKit, DeepChem, Zarr, Biopython, scikit-allel, and R’s DESeq2 for a comprehensive analysis suite in cancer research and drug discovery.

## Step-by-Step Guide to Integrate R Packages in Conda Environment YAML

1. **Specify R and Bioconductor packages** in your environment YAML file, along with Python, rpy2, and any other dependencies.
2. **Install packages directly** via conda channels (usually `conda-forge` or `bioconda`). If direct installation isn’t possible, include an R script for installation.

### Modify the `environment.yml` File

Add R, rpy2, and Bioconductor packages to ensure installation within the Conda environment.

```yaml
name: bioinformatics_env  # Environment name
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.8
  - r-base              # Installs standalone R
  - r-essentials        # Core R packages
  - r-devtools          # Required for Bioconductor package management
  - r-rmarkdown         # For R Markdown (optional)
  - rpy2                # Python-R interface
  - bioconductor-deseq2 # Direct install of DESeq2
  - bioconductor-edger  # Direct install of edgeR
  - pandas
  - scikit-learn
  - h5py
  - numpy
```

## Project Structure

This directory structure organizes the workflows, data, and DESeq2 integration for bioinformatics.

```bash
informatics_project/
├── chemical_informatics/
│   ├── data/
│   │   ├── chemical_data.h5               # HDF5 for cheminformatics data
│   │   └── chembl_data.csv                # Queried data from ChEMBL
│   ├── scripts/
│   │   ├── fetch_chembl_data.py           # Fetch data from ChEMBL
│   │   ├── rdkit_analysis.py              # RDKit analysis and descriptor calculations
│   │   └── deepchem_prediction.py         # DeepChem model training and predictions
│   └── README.md                          # Documentation for chemical workflow
├── genomic_informatics/
│   ├── data/
│   │   ├── genomic_data.zarr              # Zarr storage for large genomic data
│   │   ├── cancer_sequences.fasta         # Cancer-related sequences in FASTA
│   │   └── cancer_variants.vcf            # VCF for variant data
│   ├── scripts/
│   │   ├── tcga_assembler_integration.py  # Download data from TCGA
│   │   ├── biopython_sequence_analysis.py # Sequence analysis using Biopython
│   │   └── scikit_allel_variant_analysis.py # Variant analysis with scikit-allel
│   └── README.md                          # Documentation for genomic workflow
├── r_integration/
│   ├── install_r_packages.R               # Script to install DESeq2 and edgeR in R
│   └── deseq2_analysis.py                 # Python code for DESeq2 integration with rpy2
├── config/
│   └── environment.yml                    # Conda environment configuration
└── README.md                              # Main project documentation
```

## Code Files

### Chemical Informatics Workflow

   1. `fetch_chembl_data.py` – Fetch data from ChEMBL.

   ```python
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
   ```

   2. `rdkit_analysis.py` – Uses RDKit for descriptor calculations and stores data in HDF5.

   ```python
   from rdkit import Chem
   from rdkit.Chem import Descriptors
   import h5py
   import csv

   def calculate_descriptors(input_file, output_hdf5):
       weights = []
       with open(input_file, 'r') as file:
           reader = csv.reader(file)
           next(reader)
           for row in reader:
               mol = Chem.MolFromSmiles(row[0])
               if mol:
                   weights.append(Descriptors.MolWt(mol))

       with h5py.File(output_hdf5, 'w') as f:
           f.create_dataset('mol_weights', data=weights)

   if __name__ == "__main__":
       calculate_descriptors('../data/chembl_data.csv', '../data/chemical_data.h5')
   ```

   3. `deepchem_prediction.py` – Trains a model with DeepChem.

   ```python
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
   ```

## Project Structure and Execution Guide

## Genomic Informatics Workflow

   1. `tcga_assembler_integration.py` – Downloads data from TCGA.

   ```bash
   # Clone TCGA-Assembler-2 and download LUAD data
   !git clone https://github.com/compgenomics/TCGA-Assembler-2.git
   !python TCGA-Assembler-2/Module1_DownloadData.py "LUAD"
   ```

   2. `biopython_sequence_analysis.py` – Analyzes sequences using Biopython.

   ```python
   from Bio import SeqIO

   def parse_sequences(fasta_file):
       for record in SeqIO.parse(fasta_file, "fasta"):
           print(f"ID: {record.id}")
           print(f"Sequence: {record.seq}")

   if __name__ == "__main__":
       parse_sequences('../data/cancer_sequences.fasta')
   ```

   3. `scikit_allel_variant_analysis.py` – Performs variant analysis using scikit-allel.

   ```python
   import allel

   def analyze_variants(vcf_file):
       callset = allel.read_vcf(vcf_file)
       allele_counts = callset['calldata/GT']
       freqs = allel.GenotypeArray(allele_counts).to_allele_counts().to_frequencies()
       print("Allele Frequencies:", freqs)
       
       high_freq_variants = [variant for variant, freq in zip(callset['variants/ID'], freqs) if freq[1] > 0.5]
       print("High-frequency variants:", high_freq_variants)

   if __name__ == "__main__":
       analyze_variants('../data/cancer_variants.vcf')
   ```

## R Integration for Differential Expression Analysis

   1. `install_r_packages.R` – Install DESeq2 and edgeR in R.

   ```r
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")

   BiocManager::install("DESeq2")
   BiocManager::install("edgeR")
   ```

   2. deseq2_analysis.py – Runs DESeq2 from Python using rpy2.

   ```python
   import rpy2.robjects as ro
   from rpy2.robjects import pandas2ri
   from rpy2.robjects.packages import importr

   # Activate pandas-to-R DataFrame conversion
   pandas2ri.activate()

   # Import DESeq2 in R
   base = importr('base')
   deseq2 = importr('DESeq2')

   # Assume `expression_data` and `sample_info` are pandas DataFrames
   expression_data_r = pandas2ri.py2rpy(expression_data)
   sample_info_r = pandas2ri.py2rpy(sample_info)

   # Run DESeq2 analysis
   ro.r('''
   dds <- DESeqDataSetFromMatrix(countData = expression_data,
                                 colData = sample_info,
                                 design = ~ condition)
   dds <- DESeq(dds)
   res <- results(dds)
   ''')

   # Convert results back to pandas DataFrame
   deseq2_results_r = ro.r('as.data.frame(res)')
   deseq2_results = pandas2ri.rpy2py(deseq2_results_r)
   print(deseq2_results.head())
   ```

## Execution Instructions

## 1. Run the Chemical Informatics Workflow

```bash
# Step 1: Fetch ChEMBL Data
python chemical_informatics/scripts/fetch_chembl_data.py

# Step 2: RDKit analysis and store in HDF5
python chemical_informatics/scripts/rdkit_analysis.py

# Step 3: Run DeepChem predictions
python chemical_informatics/scripts/deepchem_prediction.py
```

## 2. Run the Genomic Informatics Workflow

```bash
# Step 1: Download data from TCGA using TCGA-Assembler-2
bash genomic_informatics/scripts/tcga_assembler_integration.sh

# Step 2: Sequence Analysis
python genomic_informatics/scripts/biopython_sequence_analysis.py

# Step 3: Variant Analysis
python genomic_informatics/scripts/scikit_allel_variant_analysis.py
```

## 3. Run the Genomic Informatics Workflow

```bash
python r_integration/deseq2_analysis.py
```

This setup provides a complete, integrated bioinformatics and cheminformatics pipeline, allowing easy analysis of gene expression, chemical properties, and differential expression through a streamlined R-Python interface.

### License

This project is licensed under the MIT License.
