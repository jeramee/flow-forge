# Script for scikit-allel variant analysis
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
