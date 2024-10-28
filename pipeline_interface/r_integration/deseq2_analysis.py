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
