#!/usr/bin/env python3.6

import pandas as pd
import numpy as np
import scipy
from scipy import stats
import sys
import copy

"""

Programs Mines All Significant Multi-Dimensional Adjusted Results
and performs a Fisher Exact Test on the Feed Efficiency
Status of the Sampels (HFE vs LFE) for all the Significant ASE samples
for a variant.

After implementing the Fisher Exact Test, p-values are adjusted using
the Benjamini-Hochberg method. Program then annotates all the results
with variant effects based on output from VEP (applied to testable
variants and non-testable). 


REQUIRED INPUT FILES

1. Testable_VEP_Results.txt
- All the annotated data from VEP

2. sig_multi_dim_adj_results.txt
- All the significant multi-dimensional adjusted results for a specific tissue

3. Meta Input File
- Meta Data File about the chicken samples (tells which samples are HFE vs LFE)

REQUIRED INPUT SETTINGS

1. Project Name
- Tells program how to parse the meta data file to identify the correct samples

2. Tissue
- Tells program how to parse the meta data file to identify the correct samples

3. Sample Minimum
- Minimum number of samples to considered when implementing the Fisher Exact Test 

OUTPUT FILES

1. fisher_exact_test_results.txt
- Testable variants for Fisher Exact Test with corresponding p-values

2. non_testable_sig_ase_results.txt
- Non-testable variants for the Fisher Exact Test
- Fail due to minimal sample number required or resulting p-value is NaN

3. Final_Annot_Results_Fisher_Exact_Test.txt (FINAL FILE OF INTEREST)
- Final testable variants with corresponding Fisher Exact Test p-values along with adjusted p-values (Benjamini Hochberg)
and corresponding variant effect information from VEP

4. ANNOT_non_testable_sig_ase_results.txt
- Non-testable variants annotated with corresponding variant effect information 


"""


def create_chicken_meta_data_dict(meta_data_input_file):

    """

    Create the empty meta data dictionary for storing all
    the meta data results for analysis of the data

    : Param meta_data_input_file: Meta data file being examined

    : Return empty_meta_dict: Empty meta data dictionary to be filled


    """

    # Create lists to store variables
    project_list = []
    tissue_list = []

    # Open the input file
    input_file = open(meta_data_input_file, 'r')

    # Start looping over the file
    for line in input_file:

        # Skip the header line
        if line.startswith("Old"):
            continue

        # Deal with the rest of the data
        else:

            # Split the data
            split_data = line.split("\t")

            # Get the variables of interest
            project = split_data[2]
            tissue = split_data[3]

            # Add to the lists
            project_list.append(project)
            tissue_list.append(tissue)

    # Convert lists to sets and back to lists
    project_list = list(set(project_list))
    tissue_list = list(set(tissue_list))

    #### Create Inner Tissue Dictionary ###
    # Create empty inner dict
    inner_tissue_dict = {}

    # Loop over tissues
    for tissue in tissue_list:

        # Update the inner dictionary with empty dictionaries
        inner_tissue_dict.update({tissue: {}})
    ########################################

    ######### Meta Data Dictionary #########
    # Create a dictionary
    empty_meta_dict = {}
    
    # Loop over projects and in tissues for dictionary
    for project in project_list:

        # Make a deep copy of the
        deep_copy_inner_tissue_dict = copy.deepcopy(inner_tissue_dict)

        # Update the Meta Dictionary
        empty_meta_dict.update({project: deep_copy_inner_tissue_dict})
    #########################################

    # Close the file
    input_file.close()

    return(empty_meta_dict)


def  extract_chicken_meta_data(meta_data_input_file, empty_meta_data_dict):

    """

    Extract the meta data about the chickens with sample IDs as keys
    and feed efficiency status as values

    : Param meta_dta_input_file: Meta data file about the chickens
    : Param empty_meta_data_dict: Empty meta data dictionary to be filled

    : Return chicken_meta_dict: Filled in meta data dictionary

    """

    # Open the input file
    input_file = open(meta_data_input_file, 'r')

    # Rename Dictionary to prevent confusion
    chicken_meta_dict = empty_meta_data_dict

    # Start Looping over liens of the file
    for line in input_file:

        # Remove the new line for safety
        line = line.rstrip("\n")

        # Skip the header line
        if line.startswith("Old"):
            continue

        # Deal with actual data
        else:

            # Split the data
            split_data = line.split("\t")

            # Get the variables of interest
            correct_id = split_data[1]
            project = split_data[2]
            tissue = split_data[3]
            line = split_data[4]
            feed_status = split_data[6]

            # Update the chicken meta data dictionary
            chicken_meta_dict[project][tissue].update({correct_id: {'line': line,
                                                                     'feed_status': feed_status}})
            
    # Close the file
    input_file.close()

    return (chicken_meta_dict)


def extract_variant_information(vep_annotation_results_file_name):

    """

    Extract the variant information from the testable VEP results
    file created using VEP

    : Param vep_annotation_results_file_name: Name of file being parsed

    : Return variant_info_dict: Dictionary of variant information

    """

    # Create the dictionary to store the results
    variant_info_dict = {}

    # Open the file
    input_file = open(vep_annotation_results_file_name, 'r')

    # Loop over the file
    for line in input_file:

        # Remove the new line for safety
        line = line.rstrip("\n")

        # Skip the header line
        if line.startswith("#"):
            continue

        # Deal with Rest of Data
        else:

            # Split the lie
            split_line = line.split("\t")

            # Get variables of interest
            rs_id = split_line[0]
            location = split_line[1]
            consequence = split_line[3]
            impact = split_line[4]
            gene_symbol = split_line[5]

            # Add results to dictionary
            variant_info_dict.update({location: {'rs_id': rs_id,
                                                 'consequence': consequence,
                                                 'impact': impact,
                                                 'gene_symbol': gene_symbol}})

    # Close the file
    input_file.close()

    return(variant_info_dict)


def create_sig_counter_ase_dict(sample_names):

    """

    Create counter dictionary based on the total number
    of samples

    : Param sample_names: list of all the samples

    : Return counter_sig_ase_dict; Dictionary to count the number of Sig ASE samples

    """

    # Create a dictionary to store values
    counter_sig_ase_dict = {}

    # Get the total sample counter
    total_sample_count = len(sample_names)

    # Loop over the range
    for value in range(1, total_sample_count + 1):

        # Update the dictionary with values
        counter_sig_ase_dict.update({value: 0})

    return(counter_sig_ase_dict)


def create_sample_index_dict(sample_names):

    """

    Creates a sample index dictionary, scrubbing the names
    to remove extra information about tissue source, so names
    can be used to look up meta-data results

    : Param sample_names: List of sample names from file header

    : Return sample_index_dict: Dictionary of sample names (key-index, values-names)

    """

    # Create dictionary to store names
    sample_index_dict = {}

    # Start a index counter
    index_counter = 0

    # Loop over the samples
    for sample in sample_names:

        # Split the name on the underscore
        split_sample_name = sample.split("_")

        # Take the last entry which the sample ID
        sample_id = split_sample_name[-1]

        # Add to the sample index dictionary
        sample_index_dict.update({index_counter: sample_id})

        # Add to the index counter
        index_counter += 1
        
    return (sample_index_dict)


def analyze_ase_behavior(chicken_meta_dict, project, tissue,
                         multi_dim_adj_file_name, group_sample_min):

    """

    Analyze a multi-dimensional adjusted file with significant samples for ASE

    : Param chicken_meta_dict: Dictionary of all the meta data for the samples
    : Param project: Project name used to looking up values in the chicken meta data dictionary
    : Param tissue: Tissue being examined (used for dictionary lookup)
    : Param multi_dim_adj_file_name: Name of the file being examined
    : Param group_sample_min: Minimum number of samples that must be in each group to perform analysis (HFE vs LFE)

    : Return: output_file_name: Name of the results file from analysis, file needed
        for p-value adjustment and merging in variant meta data
    : Return: failure_output_file_name: Name of variants that could not be tested with hypothesis test

    """

    print ("Analyzing Multi-Dimensional Adjusted Samples")

    # Open the input file
    input_file = open(multi_dim_adj_file_name, 'r')

    # Create an output file for results
    output_file_name = 'fisher_exact_test_results.txt'

    # Create a failure output file
    failure_output_file_name = 'non_testable_sig_ase_results.txt'

    # Open the files for writing
    output_file = open(output_file_name, 'w')
    non_testable_output_file = open(failure_output_file_name, 'w')

    # Write Headers Testable Data
    output_file.write("Chrom\tPos\tID\tRef\tAlt\tTotal_Biallelic\tTotal_ASE\tASE_Alleles\tASE_Alleles_HFE\tASE_Alleles_LFE\t"
                      + "Sig_ASE_HFE\tNonSig_HFE_Biallelic\tSig_ASE_LFE\tNonSig_LFE_Biallelic\tFisherExactTest_Pvalue\n")

    # Write Header to Non_Testable Data
    non_testable_output_file.write("Chrom\tPos\tID\tRef\tAlt\tTotal_Biallelic\tTotal_ASE\tASE_Alleles\tASE_Alleles_HFE\tASE_Alleles_LFE\t"
                                   + "Sig_ASE_HFE\tNonSig_HFE_Biallelic\tSig_ASE_LFE\tNonSig_LFE_Biallelic\tFailure_Cause\n")
    
    # Start looping over the file
    for line in input_file:

        # Remove the new line from file
        line = line.rstrip("\n")
        # Hidden Tab After Last Entry (Bug in Coding)
        line = line.rstrip("\t")

        # Flag the header
        if line.startswith("#CHROM"):

            # Split the line
            line_split = line.split("\t")

            # Get sample names
            sample_names = line_split[9:]

            # Create sample index dictionary
            sample_index_dict = create_sample_index_dict(sample_names)

            # Create a Counter dictionary of ASE (use total samples count)
            counter_sig_ase_dict = create_sig_counter_ase_dict(sample_names)

        else:

            # Replaces all quotes in lines (Excels corrupts txts)
            line = line.replace('"', '')
            line = line.replace("'", "")

            # Split the line
            line_split = line.split("\t")

            # Get the variables of interest
            chromo = line_split[0]
            position = line_split[1]
            rs_id = line_split[2]
            reference = line_split[3]
            alternative = line_split[4]
            
            # Get the sample data
            samples_data = line_split[9:]

            # Create Counters
            sig_ASE_HFE = 0
            non_sig_Biallelic_HFE = 0
            total_Biallelic_HFE = 0

            sig_ASE_LFE = 0
            non_sig_Biallelic_LFE = 0
            total_Biallelic_LFE = 0

            total_Biallelic_Samples = 0
            total_sig_ASE_Samples = 0

            # Create lists to store directionality
            # (Convert lists to sets later to remove duplicates)
            sig_ASE_HFE_dir_list = []
            sig_ASE_LFE_dir_list = []

            # Keeps track of index for looking up sample name
            sample_index_counter = 0
            
            # Start looping over sample
            for sample in samples_data:

                # Get the sample ID name from Sample Index Dictionary (cleaned up ID)
                sample_ID = sample_index_dict[sample_index_counter]

                # Split the sample results
                sample_results = sample.split(":")

                # Get variables
                verdict = sample_results[0]
                counts = sample_results[2]
                significance_verdict = sample_results[4]

                # Filter out all non-biallelic samples
                if verdict != "Biallelic":

                    # Add to the sample index counter and continue
                    sample_index_counter += 1
                    continue

                # Biallelic Samples
                else:

                    # Get Ref and Alternative Counts
                    split_counts = counts.split(",")
                    ref_count = split_counts[0]
                    alt_count = split_counts[1]

                    # Get the Feed Efficiency status
                    feed_efficiency_status = chicken_meta_dict[project][tissue][sample_ID]['feed_status']

                    ################# High Feed Efficiency (HFE) Results ############################
                    # Add to counters based on results
                    if feed_efficiency_status == 'HFE' and significance_verdict == 'Fail':

                        # Add to specific counters
                        total_Biallelic_HFE += 1
                        total_Biallelic_Samples += 1
                        non_sig_Biallelic_HFE += 1
                        
                        # Add to the sample index counter and continue
                        sample_index_counter += 1
                        continue

                    elif feed_efficiency_status == 'HFE' and significance_verdict == 'Pass':

                        # Add to specific counters
                        total_Biallelic_HFE += 1
                        total_Biallelic_Samples += 1
                        # Sig Counters
                        sig_ASE_HFE += 1
                        total_sig_ASE_Samples += 1

                        # Get Directionality
                        # If Ref is the ASE Allele
                        if float(ref_count) > float(alt_count):
                            sig_ASE_HFE_dir_list.append("Ref")
                        # Alt Allele is the ASE Allele
                        else:
                            sig_ASE_HFE_dir_list.append("Alt")
                        
                        # Add to the sample index counter and continue
                        sample_index_counter += 1
                        continue

                    ################# Low Feed Efficiency (LFE) Results ############################
                    elif feed_efficiency_status == 'LFE' and significance_verdict == 'Fail':

                        # Add to specific counters
                        total_Biallelic_LFE += 1
                        total_Biallelic_Samples += 1
                        non_sig_Biallelic_LFE += 1
                        
                        # Add to the sample index counter and continue
                        sample_index_counter += 1
                        continue

                    elif feed_efficiency_status == 'LFE' and significance_verdict == 'Pass':

                        # Add to specific counters
                        total_Biallelic_LFE += 1
                        total_Biallelic_Samples += 1
                        # Sig Counters
                        sig_ASE_LFE += 1
                        total_sig_ASE_Samples += 1
                        

                        # Get Directionality
                        # If Ref is the ASE Allele
                        if float(ref_count) > float(alt_count):
                            sig_ASE_LFE_dir_list.append("Ref")
                        # Alt Allele is the ASE Allele
                        else:
                            sig_ASE_LFE_dir_list.append("Alt")
                        
                        # Add to the sample index counter and continue
                        sample_index_counter += 1
                        continue

                    # WARNING IN CASE OPTIONS ARE WRONG
                    else:

                        print ("Combo of feed efficiency status and significance not programmed correctly")
                        print ("Killing Program for Debugging")
                        sys.exit()

            # Writing Results to File Based on Findings
            # Validating Program is working correctly (comment out when done beta testing)

            ### Get the directionality of ASE alleles ###
            # Remove Duplicates and Convert To a List
            ASE_HFE_alleles_list = (list(set(sig_ASE_HFE_dir_list)))
            ASE_LFE_alleles_list = (list(set(sig_ASE_LFE_dir_list)))

            # Flag Empty Lists for HFE Alleles and Print NaN
            if len(ASE_HFE_alleles_list) == 0:
                ASE_HFE_alleles_str = 'NaN'
            # Non-Empty Lists Convert to a String
            else:
                ASE_HFE_alleles_str = ','.join(map(str, ASE_HFE_alleles_list))

            # Flag Empty Lists for LFE Alleles and Print NaN
            if len(ASE_LFE_alleles_list) == 0:
                ASE_LFE_alleles_str = 'NaN'
            # Non-Empty Lists Convert to a String
            else:
                ASE_LFE_alleles_str = ','.join(map(str, ASE_LFE_alleles_list))

            # Examine if the HFE vs LFE Alleles Match
            ase_alleles_list = list(set(sig_ASE_HFE_dir_list + sig_ASE_LFE_dir_list))
            ase_alleles_str = ','.join(map(str, ase_alleles_list))

            # Filter Out Data Based on the Minimum Sample Requirement for the Hypothesis Testing (Used Defined)
            if (total_Biallelic_HFE < group_sample_min) or (total_Biallelic_LFE < group_sample_min) :
                
                # Write Results to Non-Testable File
                non_testable_output_file.write(chromo + "\t" + position + "\t" + rs_id + "\t" + reference + "\t"
                                               + alternative + "\t" + str(total_Biallelic_Samples) + "\t" + str(total_sig_ASE_Samples) + "\t"
                                               + ase_alleles_str + "\t" + ASE_HFE_alleles_str + "\t" + ASE_LFE_alleles_str + "\t"
                                               + str(sig_ASE_HFE) + "\t" + str(non_sig_Biallelic_HFE) + "\t"
                                               + str(sig_ASE_LFE) + "\t" + str(non_sig_Biallelic_LFE) + "\t" + "Fail- Group Sample Min (HFE_Biallelic="
                                               + str(total_Biallelic_HFE) + " LFE_Biallelic=" + str(total_Biallelic_LFE)  + "\n")

                # Move Onto Next Line of Data
                continue

            # Passing Data Analyze and Produce a P-value
            else:

                # Perform the Fisher Exact Test
                # Using SciPy Built in Model: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
                fisher_exact_results = stats.fisher_exact([[sig_ASE_HFE, sig_ASE_LFE], [non_sig_Biallelic_HFE, non_sig_Biallelic_LFE]])
                p_value = fisher_exact_results[1]  
                
                # Convert p-value to a float for NumPy testing
                p_value = float(p_value)

                # Flag All P-values that are NaN
                if np.isnan(p_value) == True:

                    # Write Results to Non-Testable File
                    non_testable_output_file.write(chromo + "\t" + position + "\t" + rs_id + "\t" + reference + "\t"
                                                   + alternative + "\t" + str(total_Biallelic_Samples) + "\t" + str(total_sig_ASE_Samples) + "\t"
                                                   + ase_alleles_str + "\t" + ASE_HFE_alleles_str + "\t" + ASE_LFE_alleles_str + "\t"
                                                   + str(sig_ASE_HFE) + "\t" + str(non_sig_Biallelic_HFE) + "\t"
                                                   + str(sig_ASE_LFE) + "\t" + str(non_sig_Biallelic_LFE) + "\t" + "Fail- P-value is NaN" + "\n")
                    # Move Onto Next Line of Data
                    continue

                # Record all testable results (numeric p-values)
                else:
                    output_file.write(chromo + "\t" + position + "\t" + rs_id + "\t" + reference + "\t"
                                      + alternative + "\t" + str(total_Biallelic_Samples) + "\t" + str(total_sig_ASE_Samples) + "\t"
                                      + ase_alleles_str + "\t" + ASE_HFE_alleles_str + "\t" + ASE_LFE_alleles_str + "\t"
                                      + str(sig_ASE_HFE) + "\t" + str(non_sig_Biallelic_HFE) + "\t"
                                      + str(sig_ASE_LFE) + "\t" + str(non_sig_Biallelic_LFE) + "\t" + str(p_value) + "\n")
                    
    # Close the files
    input_file.close()
    output_file.close()
    non_testable_output_file.close()

    return (output_file_name, failure_output_file_name)
            

def extract_p_values_from_file(hypothesis_results_file_name):

    """

    Extracts the p-values from the hypothesis results file for
    adjustment

    : Param hypothesis_results_file_name: Name of the file being examined

    : Return list_pvalues: List of p-values for adjustment

    """

    # Create an empty list to store p-values
    list_pvalues = []

    # Open the file
    input_file = open(hypothesis_results_file_name, 'r')

    # Start looping over the lines of the file
    for line in input_file:

        # Remove the new line for safety
        line = line.rstrip("\n")

        # Skip the header line
        if line.startswith("Chrom"):
            continue

        # Deal with actual data
        else:

            # Split the line
            line_split = line.split("\t")

            # Get the p-value
            p_value = line_split[14]

            # Add the pvalue to the list
            list_pvalues.append(float(p_value))

    return (list_pvalues)
    

def create_final_results_file(testable_results_file_name, re_sorted_pvalues_list_dict, variant_info_dict):

    """

    Creates the final output file from the analysis, combing in all the variant data and adjusted
    p-value results

    : Param testable_results_file_name: Name of the results file with testable variants for Fisher Exact Test
    : Param re_sorted_pvalues_list_dict: Adjusted p-value dictionary
    : Param variant_info_dict: Dictionary of all the variant info to be merged into final data
    
    : Return None:

    """

    # Open the input file
    input_file = open(testable_results_file_name, 'r')

    # Output File Name
    output_file_name = "Final_Annot_Results_Fisher_Exact_Test.txt"

    # Open the output file
    output_file = open(output_file_name , 'w')

    # Write the header to the file
    # Write Headers Testable Data
    output_file.write("Chrom\tPos\tID\tGene_Symbol\tConsequence\tImpact\tRef\tAlt\tTotal_Biallelic\tTotal_ASE\tASE_Alleles\tASE_Alleles_HFE\tASE_Alleles_LFE\t"
                      + "Sig_ASE_HFE\tNonSig_HFE_Biallelic\tSig_ASE_LFE\tNonSig_LFE_Biallelic\tFisherExactTest_Pvalue\tAdjusted_Pvalue\n")

    # Create a variant Counter (Use for looking adjusted p-value)
    variant_index_counter = 0

    # Open the input file
    for line in input_file:

        # Remove the new line for safety reasons
        line = line.rstrip("\n")

        # Skip the header line
        if line.startswith("Chrom"):
            continue

        # Deal with rest of data
        else:

            # Split the line on tab
            split_line = line.split("\t")

            # Get Variables of Interest
            chromosome = split_line[0]
            position = split_line[1]
            rs_id = split_line[2]
            p_value = split_line[14]

            # Combine Chromosome with Position
            location = chromosome + ":" + position + "-" + position

            # Get all the rest of data as one list
            rest_of_data = split_line[3:]
            rest_of_data_str = "\t".join(map(str, rest_of_data))

            # Get Results from Variant Info Dictionary
            gene_symbol_var_dict = variant_info_dict[location]['gene_symbol']
            consequence_var_dict = variant_info_dict[location]['consequence']
            impact_var_dict = variant_info_dict[location]['impact']
            rs_id_var_dict = variant_info_dict[location]['rs_id']

            if rs_id != rs_id_var_dict:
                print ("Issues with RS IDs between files")
                print ("Killing Program")
                print ("")
                print ("Variant Line")
                print (line)
                sys.exit()

            # Get the p-value from the adjusted dictionary (list dictionary combo)
            p_value_results = re_sorted_pvalues_list_dict[variant_index_counter][1]

            # Break up the results
            original_pvalue_adj_dict = p_value_results['original_pvalue']
            adjusted_pvalue_adj_dict = p_value_results['adjusted_pvalue']

            # Double Check P-values Match Between Line of Data and Adjusted Dict
            if float(p_value) != float(original_pvalue_adj_dict):
                print ("Issues with Pvalues and P-value Adjustment Dictionary")
                print ("Killing Program")
                print ("")
                print ("Check variant counter is working correctly")
                print (line)
                print (p_value)
                print ("")
                print (p_value_results)
                print (original_pvalue_adj_dict)
                sys.exit()

            # Write Results to the file
            output_file.write(chromosome + "\t" + position + "\t" + rs_id + "\t" + gene_symbol_var_dict + "\t"
                              + consequence_var_dict + "\t" + impact_var_dict + "\t" + rest_of_data_str + "\t"
                              + str(adjusted_pvalue_adj_dict) + "\n")
            
            # Add to the variant counter
            variant_index_counter += 1

    # Close the files
    input_file.close()
    output_file.close()


def merge_annotation_failure_file(failure_results_file_name, variant_info_dict):

    """

    Merge the annotation information into the failure variants (not testable for hypothesis test)

    : Param failure_results_file_name: Name of the failure file
    : Param variant_info_dict: Variant information

    : Return None:

    """

    # Create new output file name
    output_file_name = 'ANNOT_non_testable_sig_ase_results.txt'

    # Open the output file
    output_file = open(output_file_name, 'w')

    # Open the input file
    input_file = open(failure_results_file_name, 'r')

    # Write Headers Testable Data
    output_file.write("Chrom\tPos\tID\tGene_Symbol\tConsequence\tImpact\tRef\tAlt\tTotal_Biallelic\tTotal_ASE\t"
                      + "ASE_Alleles\tASE_Alleles_HFE\tASE_Alleles_LFE\t"
                      + "Sig_ASE_HFE\tNonSig_HFE_Biallelic\tSig_ASE_LFE\tNonSig_LFE_Biallelic\tFailure_Cause\n")

    # Open the input file
    for line in input_file:

        # Remove the new line for safety reasons
        line = line.rstrip("\n")

        # Skip the header line
        if line.startswith("Chrom"):
            continue

        # Deal with rest of data
        else:

            # Split the line on tab
            split_line = line.split("\t")

            # Get Variables of Interest
            chromosome = split_line[0]
            position = split_line[1]
            rs_id = split_line[2]

            # Combine Chromosome with Position
            location = chromosome + ":" + position + "-" + position

            # Get all the rest of data as one list
            rest_of_data = split_line[3:]
            rest_of_data_str = "\t".join(map(str, rest_of_data))

            # Get Results from Variant Info Dictionary
            gene_symbol_var_dict = variant_info_dict[location]['gene_symbol']
            consequence_var_dict = variant_info_dict[location]['consequence']
            impact_var_dict = variant_info_dict[location]['impact']
            rs_id_var_dict = variant_info_dict[location]['rs_id']

            if rs_id != rs_id_var_dict:
                print ("Issues with RS IDs between files")
                print ("Killing Program")
                print ("")
                print ("Variant Line")
                print (line)
                sys.exit()

            # Write Results to the file
            output_file.write(chromosome + "\t" + position + "\t" + rs_id + "\t" + gene_symbol_var_dict + "\t"
                              + consequence_var_dict + "\t" + impact_var_dict + "\t" + rest_of_data_str + "\n")
            
    # Close the files
    input_file.close()
    output_file.close()


    return()


##################################################################################################################################
############################################## FDR CODE BUILT INTO PROGRAM #######################################################     
##################################################################################################################################


def create_pvalues_fdr_results_dict(pvalues_list):


    '''
    Creates a dictionary inside a list to allow corrected
    FDR pvalues to be stored and also allow for values to changed
    and updated in the inside dictionary
    : parameters pvalues_list: list of input pvalues
    : return pvalues_fdr_results_list: list dictionary where the list
        value is the original p-value, followed index value then dictionary
    '''

    # Create list dictionary combo to store results
    pvalues_list_dict = []

    # Start index counter
    x = 0

    # Start looping over pvalues
    for pvalue in pvalues_list:

        # For each entry add a value as an outside list with a dictionary inside
        # Allow sorting of the outside value, very important
        pvalues_list_dict.append([pvalue, x, {'index_in_list': x,
                                               'original_pvalue': pvalue,
                                               'adjusted_pvalue': 'nan'}])
        # Add to the index counter
        x += 1

    return (pvalues_list_dict)


def reorder_sorted_pvalues_list_dict(sorted_pvalues_list_dict):

    """
    Takes in the final list-dictionary with corrected p-values
    and replaces pvalue "list" value with the original list index value
    so the list-dictionary can be re-sorted to match the original list
    order.
    : parameters sorted_pvalues_list_dict: list-dictionary with all entries
        where the list value is the original p-value
    : return re_sorted_pvalues_list_dict: list-dictionary where the
        original the list value is the original index value and has been
        re-sorted by that value
    """

    # Loop over all the entries in the list dictionary
    for entry in sorted_pvalues_list_dict:

        # Delete the first entry in the list
        # To get to the original index value
        del entry[0]

    # Sort the list (outside value using original index value)
    re_sorted_pvalues_list_dict = sorted(sorted_pvalues_list_dict)

    return (re_sorted_pvalues_list_dict)


def fdr_correction(pvalues_list):

    """
    Function performs FDR correction of the p-values
    using Benjamini-Hochberg (1995) which sorts the list of pvalues
    and then determines the p-value correction based on the rank and following
    equation (p-value x NumbTest / p-value_Rank)
    :param values: list of pvalues to correct
    :return p_value_dict: dictionary of corrected pvalues
    
    """

    # Create a list of pvalues with a dictionary inside each entry
    pvalues_list_dict = create_pvalues_fdr_results_dict(pvalues_list)

    # Sort the list in Reverse (outside value uses raw pvalue)
    # if a duplicate is found, goes to next value (aka index value),
    # index value required to prevent sort from breaking
    reverse_sorted_pvalues_list_dict = sorted(pvalues_list_dict, reverse=True)

    #Position Movement Counter
    i = 0

    # Get the total number of pvalues analyzed
    total_pvalues = len(pvalues_list_dict)

    # Start looping over list dictionary 
    for entry in reverse_sorted_pvalues_list_dict:

        # Get the sorted position of p-value (most significant to least significant),
        # opposite the current list
        p_value_sorted_position = total_pvalues - i

        # Last value in list (no need for following calculation)
        if i == 0:
            
            # Get corrected pvalue for position i
            # Important: No correction occurs for position because last value in sorted list
            fdr_adj_pvalue = round((entry[0] * total_pvalues / p_value_sorted_position), 8)

            # Adjust p-values may be greater than 1
            if fdr_adj_pvalue > 1:
                fdr_adj_pvalue = 1

            # Update the dictionary value with corrected pvalue
            entry[2]['adjusted_pvalue'] = fdr_adj_pvalue

            # Increment i value
            i+=1

        else:
            
            # Get corrected pvalue for current entry
            fdr_adj_pvalue = round((entry[0] * total_pvalues / p_value_sorted_position), 8)

            # Get pvalue for the prior entry, which has a less significant original pvalue, which after
            # adjustment could be more significant and needs to replace the adjusted pvalue
            # of the current entry if this occurs
            prior_entry_adj_pvalue = reverse_sorted_pvalues_list_dict[i -1][2]['adjusted_pvalue']

            # Correct adjusted pvalue if lower ranked pvalue is more signficant
            if prior_entry_adj_pvalue < fdr_adj_pvalue:
                fdr_adj_pvalue = prior_entry_adj_pvalue

            # Update the dictionary value with corrected pvalue
            entry[2]['adjusted_pvalue'] = fdr_adj_pvalue

            # Increment i value
            i+=1

    # Replace first value in list with index to re-order values to match orginal data
    re_sorted_pvalues_list_dict = reorder_sorted_pvalues_list_dict(reverse_sorted_pvalues_list_dict)

    return(re_sorted_pvalues_list_dict)


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################


def main():

    ##################### Input Files/Settings ##################
    meta_data_input_file = r'C:\Users\mjtom\Desktop\Final_VADT_Results_For_Paper\Masked_Genome_Results\Proportional_Hypothesis_Testing\Meta_Data_Chickens_No_Duplicates.txt'

    vep_annotation_results_file_name = 'Testable_VEP_Results_Breast_Muscle.txt'

    multi_dim_adj_file_name = 'BM_sig_multi_dim_adj_results.txt'

    project = 'Feed Efficiency'

    tissue = 'Breast Muscle'

    # Minimum Number of Samples for Each Group Tested
    group_sample_min = 5

    
    ##################################################################################################################################
    ################################################## DO NOT CHANGE BELOW ###########################################################
    ##################################################################################################################################
    
    print ("Starting Fisher Exact Test Program")

    # Create Meta Data Dict for Storing Data
    empty_meta_data_dict = create_chicken_meta_data_dict(meta_data_input_file)

    # Extract the meta data about the chickens with IDs as dictionary keys
    chicken_meta_dict = extract_chicken_meta_data(meta_data_input_file, empty_meta_data_dict)

    # Extract Variant Information
    variant_info_dict = extract_variant_information(vep_annotation_results_file_name)

    # Loop over analyze multi-dimensional adjusted results
    analysis_results = analyze_ase_behavior(chicken_meta_dict, project, tissue, multi_dim_adj_file_name, group_sample_min)
    testable_results_file_name = analysis_results[0]
    failure_results_file_name = analysis_results[1]

    # Extract P-values for Adjusting
    list_pvalues = extract_p_values_from_file(testable_results_file_name)

    # Perform FDR Correction on list
    re_sorted_pvalues_list_dict = fdr_correction(list_pvalues)

    # Merge all Variant Info and Adjusted P-values into the Results File
    create_final_results_file(testable_results_file_name, re_sorted_pvalues_list_dict, variant_info_dict)

    # Merge Variant Results into Non-Testable Results for Record Keeping
    merge_annotation_failure_file(failure_results_file_name, variant_info_dict)
    
    print ("Done Running Program")

                
main()























