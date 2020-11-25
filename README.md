# HFEvsLFE_FisherExactTest
Analyzing Meta Data (HFEvsLFE) using FisherExactTest
Programs Mines All Significant Multi-Dimensional Adjusted Results
and performs a Fisher Exact Test on the Feed Efficiency
Status of the Sampels (HFE vs LFE) for all the Significant ASE samples
for a variant.

After implementing the Fisher Exact Test, p-values are adjusted using
the Benjamini-Hochberg method. Program then annotates all the results
with variant effects based on output from VEP (applied to testable
variants and non-testable). 


REQUIRED INPUT FILES

1. Testable_VEP_Results.txt\
- All the annotated data from VEP

2. sig_multi_dim_adj_results.txt
All the significant multi-dimensional adjusted results for a specific tissue

3. Meta Input File
Meta Data File about the chicken samples (tells which samples are HFE vs LFE)

REQUIRED INPUT SETTINGS

1. Project Name
Tells program how to parse the meta data file to identify the correct samples

2. Tissue
Tells program how to parse the meta data file to identify the correct samples

3. Sample Minimum
Minimum number of samples to considered when implementing the Fisher Exact Test 

OUTPUT FILES

1. fisher_exact_test_results.txt
Testable variants for Fisher Exact Test with corresponding p-values

2. non_testable_sig_ase_results.txt
Non-testable variants for the Fisher Exact Test
Fail due to minimal sample number required or resulting p-value is NaN

3. Final_Annot_Results_Fisher_Exact_Test.txt (FINAL FILE OF INTEREST)
Final testable variants with corresponding Fisher Exact Test p-values along with adjusted p-values (Benjamini Hochberg)
and corresponding variant effect information from VEP

4. ANNOT_non_testable_sig_ase_results.txt
Non-testable variants annotated with corresponding variant effect information 
