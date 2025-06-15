# ESL-PSC functions

import os, subprocess, math, re, time, datetime, shutil, argparse, sys
from collections import defaultdict, Counter
from Bio import SeqIO
import numpy as np
from . import sps_density
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations, chain
from statistics import median


def parse_args_with_config(parser, raw_args=None):
    '''will parse args in the esl-psc_config.txt file if it exists and also
    get args from sys.argv[1:] and return the args namespace'''
    cmdline = raw_args if raw_args is not None else sys.argv[1:]
    if os.path.exists("esl_psc_config.txt"): # check for the config file
        print("getting args from esl_psc_config.txt...")
        with open("esl_psc_config.txt") as file:
            args, remaining = parser.parse_known_args(file.read().split()
                                                      + cmdline)
            if len(remaining) > 0:
                print("unrecognized args: ", remaining) 
    else: 
        print("did not find an esl_psc_config.txt in this directory")
        args = parser.parse_args(cmdline)
    # Point esl_main_dir at the *project root* (one level above this package)
    this_dir = os.path.dirname(os.path.abspath(__file__))
    args.esl_main_dir = os.path.abspath(os.path.join(this_dir, os.pardir))
    # sanity-check: make sure the ESL binaries are actually there
    if not os.path.exists(os.path.join(args.esl_main_dir, "bin", "preprocess")):
        raise FileNotFoundError(
            f"Expected ESL binaries in {args.esl_main_dir!s}/bin but none were found"
        )
    # Ensure esl_inputs_outputs_dir is set if not already provided
    if not hasattr(args, 'esl_inputs_outputs_dir') or not args.esl_inputs_outputs_dir:
        args.esl_inputs_outputs_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "preprocessed_data_and_outputs/")

    # Conditional error for missing species_pheno_path when plot generation is requested
    plot_requested = getattr(args, 'make_sps_plot', False) or getattr(args, 'make_sps_kde_plot', False)
    species_pheno_path_missing = not getattr(args, 'species_pheno_path', None)

    if plot_requested and species_pheno_path_missing:
        raise ValueError("Error: A species phenotype file is required to make species plots or KDE plots.")

    # Check for relative paths and adjust them to be absolute
    path_args = [
        'esl_inputs_outputs_dir', 'species_pheno_path', 'prediction_alignments_dir', 
        'output_dir', 'canceled_alignments_dir', 'response_file', 'species_groups_file', 
        'limited_genes_list', 'alignments_dir', 'response_dir', 'esl_main_dir'
    ]

    # Convert path arguments to absolute paths
    for path_arg in path_args:
        path = getattr(args, path_arg)
        if path:
            setattr(args, path_arg, os.path.abspath(path))

    return args

def is_fasta(file_name):
    return file_name.endswith(('.fa','.fas','.fasta'))

def file_lines_to_list(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    line_list = [line.strip() for line in lines]
    return line_list

def get_species_to_check(response_matrix_path, check_order = False):
    '''takes a full path to a response matrix and returns a list of species
    in which the even (0-based) indices should be +1 and odds -1.
    If check_order is True, assertions will ensure the 1, -1, 1, -1 order.
    '''
    species_list = []
    response_lines = file_lines_to_list(response_matrix_path)
    for line_index, line in enumerate(response_lines):
        species_name, response_number = line.split('\t')

        # check to make sure ordering is correct
        if check_order:
            if line_index % 2 == 0:
                # even indices are 0, 2, 4, i.e. 1st, 3rd, 5th
                assert int(response_number) == 1 or int(response_number) == 0
            else:
                assert int(response_number) == -1 or int(response_number) == 0
        #this shouldn't happen unless it's a pair
        if int(response_number) != 0: 
            species_list.append(species_name)
    return species_list

def make_taxa_list(alignments_dir_path):
    '''takes a directory path and checks files in the directory for
    lines starting with '>' and assumes those are species IDs and returns
    a list of all of the unique ones found in all files, thus all taxa in
    appearing in the alignments will be represented'''
    previous_dir = os.getcwd()
    os.chdir(alignments_dir_path)
    subprocess.run('grep -h ">" * | sort | uniq | sed "s/>//" > temp___.txt',
                   shell=True, check=True)
    taxa_list = file_lines_to_list('temp___.txt')
    os.remove('temp___.txt')
    os.chdir(previous_dir)
    return taxa_list

def get_pos_and_neg_list(position_list):
    ''' takes an ordered list of positions at a site in which the order is
    1, -1, 1, -1, etc  and returns two lists, a list of the positive class
    AAs and one for the negative AAs in the same order as input
    '''
    # split into trait-positive and trait negative lists
    pos_position_list = position_list[::2]
    neg_position_list = position_list[1::2]
    return pos_position_list, neg_position_list

def make_input_pheno_dict(list_of_species_names):
    ''' takes an ordered list of species assuming 1, -1, 1, -1, etc
    returnes a dictionary of speices namesas keys and their response values
    '''
    pheno_dict = {}
    for index, species_name in enumerate(list_of_species_names):
        if index % 2 == 0: # even indices are 0, 2, 4, i.e. 1st, 3rd, 5th...
            pheno_dict[species_name] = 1
        else:
            pheno_dict[species_name] = -1
    return pheno_dict

def penalty_function_maker(penalty, penalty_type = "linear"):
    '''define how to calculate the variable penalty could be square, cube,
    etc.
    '''
    # will loop through different powers of the num variable sites
    if penalty_type == "linear":
        def temp_function_version(number_variable_sites):
            return number_variable_sites + penalty
    elif penalty_type == "sqrt":
        def temp_function_version(number_variable_sites):
            return pow(number_variable_sites, 0.5) # does sqrt penalty 
    return temp_function_version # function to generate penalty

def get_gene_names(dir_or_file_path):
    '''takes a path to a directory or a text file (an ESL path file). If a
    path file is given, returns the names of the genes from that, and if a dir
    path is given, it is assumed to be an alignment directory and the gene
    names are taken from fasta file names. returns a list of gene names.'''
    if os.path.isdir(dir_or_file_path): # if a directory path is given
        return [file[:-4] for file in os.listdir(dir_or_file_path)
                if file.endswith('.fas')]
    else: # must be a path file (text file)
        gene_name_list = []
        alignment_path_lines = file_lines_to_list(dir_or_file_path)
        # loop through lines in the path file and extract the file name
        for line in alignment_path_lines:
            gene_name = os.path.split(line)[1].strip()
            gene_name_list.append(gene_name[:-4]) #remove the .fas
        return gene_name_list
        
def get_seq_records_in_order(fasta_file, species_list):
    '''takes a relative fasta file path and ordered species list and
    returns a list of seq_records in order
    '''
    records = [] # list of seqrecord for the species in this group
    # make a dict of records in order to index by species id
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # loop to get seqrecords for this gene alignment and species combo
    for species in species_list:
        records.append(record_dict[species]) # keep species in same order
    return records
    
def is_variable_site(residues):
    '''takes a tuple or list of residues; returns 1 if variable or 0 if not
    can also tak a Counter object from the collections module
    '''
    if type(residues) == Counter:
        counted = residues # if a Counter object has been given use that
    else:
        counted = Counter(residues) 
    if '-' in counted:
        del counted['-'] #1st remove counts for gaps from the allele counter
    num_aas_left_after_gaps_removed = sum(counted.values())
    if len(counted) == 0 or len(counted) == 1:
        return 0 # go on if all the same or were all gaps and thus removed
    if len(counted) == num_aas_left_after_gaps_removed:
        return 0 # skip if all singletons
    if max(counted.values()) == num_aas_left_after_gaps_removed - 1:
        return 0 # skip if one singleton
    return 1 # this means site is variable

def count_var_sites(records):
    '''takes a list of sequence record objects from an alignment and returns a
    count of the number of variable sites in the alignment using the
    is_variable_site function above
    '''
    zipped_seqs = zip(*records) # this zips together the sequences
    num_var_sites = 0
    for residues_at_site in zipped_seqs: # check positions
        num_var_sites += is_variable_site(residues_at_site)
    return num_var_sites

def get_median_var_sites(alignment_dir):
    '''takes a full path to a directory of alignments and returns the median
    number fo variable sites, ignoring alignments with zero variable sites.
    The median is rounded down to the nearest integer.
    '''
    previous_dir = os.getcwd()
    os.chdir(alignment_dir)
    alignment_list = [file for file in os.listdir() if file.endswith('.fas')]
    numbers_of_var_sites = []
    for alignment in alignment_list:
        num_var = count_var_sites(list(SeqIO.parse(alignment, "fasta")))
        if num_var == 0:
            continue # skip if no variable sites, it won't count toward median
        numbers_of_var_sites.append(num_var)
    os.chdir(previous_dir)
    return int(median(numbers_of_var_sites))

def get_pheno_dict(species_pheno_csv_path, str_phenos = False):
    '''takes a full path to a csv file that has lines like: "species, 1"
    and returns a dictionary with keys: species, values: phenotype values. If
    str_phenos = True, the pheno type values will be returned as strings.
    '''
    pheno_dict = {}
    pheno_lines = file_lines_to_list(species_pheno_csv_path)
    for line in pheno_lines:
        # split each line and assign values
        species, pheno_value = [item.strip() for item in line.split(',')]
        pheno_dict[species] = (int(pheno_value) if not str_phenos
                               else pheno_value) # give int pheno by default
    return pheno_dict

def report_elapsed_time(start_time):
    '''takes a start time (from time.time() and prints elapsed time since'''
    secs_since_start = time.time() - start_time
    print('\n', 'time elapsed: ' +
          str(datetime.timedelta(seconds = secs_since_start)) + '\n')
    return

def parse_ESL_weight_line(line: str):
    """
    Given *one* line from the new-style out_feature_weights.txt, return

        label : str   – path_without_index **or** the literal string 'y_intercept'
        weight: float – the model weight (or the intercept value)

    The input line is always exactly two tab-separated fields.
    """
    label, weight = line.rstrip("\n").split("\t")
    return label, float(weight)

def parse_ESL_weight_label(label):
    '''takes a label (str) parsed from a line of the text versions of the
    feature weights output from ESL and returns a 3-tuple of the gene name
    (str), the position (int), the amoni acid that is weighted at that
    position (str).
    '''
    # label will end in a '_' followed by the site then '_' and the AA 
    gene_name, position, aa_to_check_for = label.rsplit("_", 2)
    position = int(position)
    return gene_name, position, aa_to_check_for

def make_path_file(alignments_dir):
    '''makes a pathfile for a folder of alignemnts and leaves the file in the
    same folder. returns the path to the path file
    '''
    alignment_file_list = [file for file in os.listdir(alignments_dir)
                           if file.endswith('.fas')]
    path_file_path = os.path.join(alignments_dir, 'paths.txt')
    with open(path_file_path, 'w') as file:
        file.write('\n'.join(alignment_file_list))
    return path_file_path

def clear_existing_folder(directory_path):
    '''checks if a directory exists and then completely deletes it if it does'''
    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        shutil.rmtree(directory_path)
    elif os.path.exists(directory_path) and not os.path.isdir(directory_path):
        raise Exception("trying to check if this directory exists "
                        + directory_path + " but this is not a directory")

def make_response_files(response_dir, list_of_species_combos):
    '''takes a directory path, either existing or not, and a list of species
    combos in 1, -1, 1, -1 order and generates an ESL response file for each
    combo. returns a list of response file paths
    '''
    response_file_list = []
    clear_existing_folder(response_dir)
    os.mkdir(response_dir)
    # combinations are given sequential codes: 0, 1, 2...
    #   the same codes will be used in gap-canceled alignments & preprocess
    for index1, combo in enumerate(list_of_species_combos): 
        response_lines = []
        for index2, species in enumerate(combo):
            # assumes species greoups are listed in "+1, -1, +1, -1" order
            response = 1 if index2 % 2 == 0 else -1
            response_lines.append(species + '\t' + str(response))
        # response file names are combo_1.txt, combo_2.txt etc.
        file_path = (os.path.join(response_dir,
                                  'combo_' + str(index1) + '.txt'))
        with open(file_path, 'w') as file:
            file.write('\n'.join(response_lines))
        response_file_list.append(file_path) # add response file to list
    print('response matrix files created in: ' + response_dir)
    return response_file_list

def make_null_combos(list_of_species_combos):
    '''takes a list of species combinations (must be even number of species)
    and generates all possible null versions by flipping the order of half of
    the pairs and this reversing the responsevalues they will recieve in esl-psc
    Returns a list of all null combos for all input combos'''
    new_list_of_species_combos = []
    for combo_num, in_combo in enumerate(list_of_species_combos):
        if len(in_combo) % 4 != 0:
            raise ValueError('Must be an even number of pairs for null')
        # get a list of pairs with each one as 2-tuples 
        pair_list = [in_combo[n:n+2] for n in range(0, len(in_combo), 2)]
        # loop through all combinations of 1/2 of the indices up to
        #  the length of in_combo
        num_pairs = len(pair_list)
        # get all combinations of half of the indices up to num_pairs
        all_index_combos = combinations(range(num_pairs), num_pairs//2)
        # remove mirror image combos from the index_combos
        index_combos = []
        # for different species combos, if the same species are always the ones
        #   that get reversed it can cause problems, so flip everyother combo
        if combo_num % 2 != 0:
            all_index_combos = reversed(list(all_index_combos))
        for index_combo in all_index_combos:
            mirror_combo = tuple(set(range(num_pairs)) - set(index_combo))
            if mirror_combo not in index_combos:
                index_combos.append(index_combo)
        for index_combo in index_combos:
            # index_combo will be like (0,2,6) or (0,1,2) etc.
            # reverse the order of the indexed tuples in the combination
            null_combo = [tuple(reversed(pair)) if i in index_combo else pair
                          for i, pair in enumerate(pair_list)]
            # chain species together in a single tuple
            null_combo = list(chain(*null_combo))
            new_list_of_species_combos.append(null_combo)
    return new_list_of_species_combos
            
def run_preprocess(esl_dir_path, response_matrix_file_path, path_file_path,
                   preprocessed_input_folder, esl_inputs_folder_name,
                   use_is = True):
    '''run the esl preprocessor and move resulting folder, whose name is
    preprocessed_input_folder, to its destination,
    which is esl_inputs_folder_name
    '''
    print("Running ESL preprocess...")
    input_folder_path = os.path.split(path_file_path)[0] #path file folder
    preexisting_cwd = os.getcwd() # record current directory
    print(input_folder_path)
    os.chdir(input_folder_path) # change to correct directory for esl
    if os.path.exists(preprocessed_input_folder):
        shutil.rmtree(preprocessed_input_folder) # this shouldn't exist here
    # next check to make sure the new folder doesn't already exist in the
    #   destination esl_inputs_folder_name where it needs to be moved to
    if os.path.exists(os.path.join(esl_inputs_folder_name,
                                   preprocessed_input_folder)):
        raise Exception("A folder with the same name as the intended name for "
                        "the new preprocess folder already exists at "
                        + esl_inputs_folder_name + " where the preprocess is "
                        "supposed to be moved to. So that won't work")
    #construct command to run ESL preprocess
    print(os.getcwd())
    preprocess_command_list = [os.path.join(esl_dir_path,
                               'bin/preprocess'),
                               response_matrix_file_path,
                               path_file_path,
                               preprocessed_input_folder]
    if use_is:
        preprocess_command_list.append("is") # add this to ignore singletons
    # make sure the input file names are right including ".txt" or get seg fault
    print(' '.join(preprocess_command_list))
    subprocess.run(' '.join(preprocess_command_list), shell=True, check=True)

    # move the input folder from preprocess to its folder
    clear_existing_folder(os.path.join(esl_inputs_folder_name,
                                       preprocessed_input_folder))
    shutil.move(preprocessed_input_folder, esl_inputs_folder_name)
    os.chdir(preexisting_cwd)
    return

def rmse_range_pred_plots(pred_csv_path, title, pheno_names = None,
                          min_genes = 0, plot_type = 'violin'):
    '''calls sps_density.create_sps_plot for each of a series of RMSE ranks.
    pheno_names is a tuple with the +1 pheno name first and the -1 pheno 2nd,
    or if it isn't given the defaults will be used.  plot type can be 'kde' or
    violin (violin will assign anything > 1 or < -1 to 1 and -1 respectively
    by default to make it easier to see the region of overlap.
    '''
    # code borrowed largely from Louise, with some tweaks
    start_time = time.time()
    print("making sps density plot figure...")
    rmse_cutoffs = [.05, .1]
    if plot_type == 'violin':
        fig, axes = plt.subplots(ncols = len(rmse_cutoffs), figsize=(8, 7),
                             constrained_layout=True)
    elif plot_type == 'kde':
        fig, axes = plt.subplots(ncols = 1, nrows = len(rmse_cutoffs),
                                 figsize=(8, 7))
    
    df = pd.read_csv(pred_csv_path, dtype = {'species':str,
                                             'SPS':float,
                                             'input_RMSE':float,
                                             'true_phenotype':str})
    df = df[['species', 'SPS', 'num_genes', 'input_RMSE', 'true_phenotype']]
    # save fig in same folder with the predictions CSV file
    fig_path = os.path.join(os.path.split(pred_csv_path)[0],
                            title + '_pred_sps_plot.svg')
    if not pheno_names: # set phenotype names
        pheno_names = (1, -1)
    pos_pheno_name = pheno_names[0]
    neg_pheno_name = pheno_names[1]
        
    for index, rmse_cutoff in enumerate(rmse_cutoffs):
        # create each plot
        print("making plot with MFS cutoff: " + str(rmse_cutoff))
        if plot_type == 'kde':    
            sps_density.create_sps_plot(df = df,
                                    RMSE_rank = rmse_cutoff,
                                    title = title,
                                    neg_pheno_name = neg_pheno_name,
                                    pos_pheno_name = pos_pheno_name,
                                    axes = axes[index],
                                    min_genes = min_genes)
        elif plot_type == 'violin':
            sps_density.create_sps_plot_violin(df = df,
                                    RMSE_rank = rmse_cutoff,
                                    title = title,
                                    neg_pheno_name = neg_pheno_name,
                                    pos_pheno_name = pos_pheno_name,
                                    axes = axes[index],
                                    min_genes = min_genes)
            axes[index].set_ylim([-1.1, 1.1])
            axes[index].axhline(y=0, linestyle='--', color='lightgray')
    fig.set_tight_layout(True)
    plt.savefig(fig_path)
    report_elapsed_time(start_time)
    plt.show()
    print("\npredictions density plot figure saved at: " + fig_path)

    

class ESLRunFamily():
    '''a class to contain a family of runs with the same inputs.
    the preprocessed input, input pheno dict (i.e. response matrix), and
    the variable part of the group penalty calculation (i.e. the penalty
    term will all be the same for each run within the family. however
    the sparsity parameters (lambdas) will be allowed to vary.  In order to
    use this class to run ESL, the esl_main_dir (containing the ESL executable
    must be set as a class attribute.
    '''
    
    def __init__(self, args, input_pheno_dict,
                 preprocessed_input_folder, penalty_term, input_alignments_dir,
                 label = '', gene_objects_dict = None):
        self.args = args
        self.input_alignments_dir = input_alignments_dir
        self.gene_objects_dict = gene_objects_dict # must have for preds & genes
        self.input_pheno_dict = input_pheno_dict # keys: species, values: phenos
        self.preprocessed_input_folder = preprocessed_input_folder # just name
        self.penalty_term = penalty_term # just used to label the family's term
        self.runs_list = []
        self.label = label
        self.species_combo_tag = self.get_species_tag()

    def __str__(self):
        return ('ESLRunFamily: ' + 'penalty_term: ' + str(self.penalty_term)
                + ' species: ' + self.get_species_tag() + ' '
                + (self.label if self.label != self.get_species_tag() else ''))

    def __repr__(self):
        return self.__str__()

    def get_species_tag(self):
        input_species_list = list(self.input_pheno_dict.keys())
        if len(input_species_list) > 12:
            return str(self.label) 
        elif '_' in input_species_list[0]:
            return '.'.join([name.split('_')[1][:3] for
                             name in input_species_list])
        else:
            return '.'.join([name[:3] for name in input_species_list])

    def gridsearch(self):
        if self.args.use_logspace: # use logspace if this option is given
            self.do_logspace_grid_search(self.args.initial_lambda1,
                                         self.args.initial_lambda2,
                                         self.args.final_lambda1,
                                         self.args.final_lambda2,
                                         self.args.num_log_points,
                                         only_lambda1 = self.args.lambda1_only)
        else: # use normal linear grid with lambda step increment
            self.do_esl_grid_search(self.args.initial_lambda1,
                                    self.args.initial_lambda2,
                                    self.args.final_lambda1,
                                    self.args.final_lambda2,
                                    self.args.lambda_step)

    def do_grid_search(self, initial_lambda1, initial_lambda2,
                           final_lambda1, final_lambda2, lambda_step):
        '''generates esl runs over ranges of both sparsity parameters'''
        lambda1 = initial_lambda1
        lambda2 = initial_lambda2
        while lambda1 <= final_lambda1:
            while lambda2 <= final_lambda2:
                # create an esl_run object and append it to the family run list
                self.runs_list.append(ESLRun(lambda1, lambda2, self))
                #increment lambda2
                lambda2 = round(lambda2 + lambda_step, 3)
            #increment lambda1
            lambda1 = round(lambda1 + lambda_step, 3)
            #set lambda2 back to initial
            lambda2 = initial_lambda2
        #now that all run objects are created, run esl for each
        self.run_all_runs()

    def do_logspace_grid_search(self, initial_lambda1, initial_lambda2,
                                final_lambda1, final_lambda2, num_values,
                                only_lambda1 = False):
        '''generates esl runs over logspace of both sparsity parameters.
        initial_value is the base 10 log of the 1st point. For readability,
        lambda values are rounded to 8 decimal places'''
        digits_to_round = abs(int(np.log10(initial_lambda1))) + 5 
        lambda1_values = np.logspace(np.log10(initial_lambda1),
                                     np.log10(final_lambda1),
                                     num_values)
        if only_lambda1:
            for lambda1 in lambda1_values:
                # create an esl_run object and append it to the family run list
                self.runs_list.append(ESLRun(round(lambda1, digits_to_round),
                                             0, # lambda 2 is always zero
                                             self))
            #now that all run objects are created, run esl for each
            self.run_all_runs()
            return
        lambda2_values = np.logspace(np.log10(initial_lambda2),
                                     np.log10(final_lambda2),
                                     num_values)
        for lambda1 in lambda1_values:
            for lambda2 in lambda2_values:
                # create an esl_run object and append it to the family run list
                self.runs_list.append(ESLRun(round(lambda1, digits_to_round),
                                             round(lambda2, digits_to_round),
                                             self))
        #now that all run objects are created, run esl for each
        self.run_all_runs()

    def run_all_runs(self):
        for run in self.runs_list:
            print('Lambda1: ' + str(run.lambda1) +
                  '  Lambda2: ' + str(run.lambda2))
            run.run_lasso()

    def do_all_calculations(self):
        num_runs = str(len(self.runs_list))
        for run_num, run in enumerate(self.runs_list, 1):
            print('Calculating predictions and/or weights for: ' + str(run)
                  + '\nrun ' + str(run_num) + ' of ' + num_runs
                  + ' in current grid;  time: '
                  + time.strftime("%H:%M:%S", time.localtime()))
            run.calc_preds_and_weights(only_pos_gss = self.args.only_pos_gss,
                                       skip_pred = self.args.no_pred_output,
                                       sites = self.args.show_selected_sites)


class ESLRun():
    '''a class to track info and results from each run'''
    def __init__(self, lambda1, lambda2, run_family):
        self.run_family = run_family
        self.lambda1 = lambda1
        self.lambda2 = lambda2
        self.output = (self.run_family.preprocessed_input_folder +
                       self.get_lambda_tag() + '_out_feature_weights.txt')
        self.num_included_genes = 0
        self.input_rmse = 0
        self.y_intercept = 0
        self.species_scores = None

    def __str__(self):
        return ('run_object: l1: ' +
                str(self.lambda1) + ' l2: ' + str(self.lambda2) +
                ' group_pen_term: ' + str(self.run_family.penalty_term)
                + ' ' + str(self.run_family.label))

    def __repr__(self):
        return self.__str__()

    def get_lambda_tag(self):
        '''make string to add to the output name to identify run parameters'''
        return ('_l1_' + str(round(self.lambda1,5)) +
                '_l2_' + str(round(self.lambda2,5)))

    def run_lasso(self):
        '''runs logistic lasso. assumes a relative path to the preprocess so be
        in the inputs+outputs folder. this will generate an xml file and
        a txt file with the feature weights. esl_main_dir is the directory
        that contains the bin/sg_lasso executable
        '''
        preprocessed_dir_name = self.run_family.preprocessed_input_folder
        output_name = (self.run_family.preprocessed_input_folder +
                       self.get_lambda_tag())
        # generate the command for calling ESL logistic lasso
        esl_command_list = [os.path.join(self.run_family.args.esl_main_dir,
                            'bin/sg_lasso'),
                            '-f', preprocessed_dir_name + '/feature_' +
                            preprocessed_dir_name + '.txt',
                            '-z', str(self.lambda1),
                            '-y', str(self.lambda2),
                            '-n', preprocessed_dir_name + '/group_indices_' +
                            preprocessed_dir_name + '.txt',
                            '-r', preprocessed_dir_name + '/response_' +
                            preprocessed_dir_name + '.txt',
                            '-w', output_name + '_out_feature_weights']
        # run esl
        subprocess.run(' '.join(esl_command_list), shell=True, check=True)
        
        # command to pull out feature weights and create text files
        grep_command = (r'grep -P "<item>.*</item>" ' +
                        preprocessed_dir_name + self.get_lambda_tag() +
                        r'_out_feature_weights.xml | sed' +
                        r' -re "s/.*<item>(.*)<\/item>.*/\1/" > ' +
                        r'temp_out_feature_weights.txt') 
        # creates 'temp_out_feature_weights.txt' (this is fast)          
        subprocess.run(grep_command, shell=True, check=True)
        
        # generate text output from temp file by removing zero lines
        with open('temp_out_feature_weights.txt','r') as temp_file:
            temp_lines = temp_file.readlines()
        with open(preprocessed_dir_name + '/feature_mapping_' +
                  preprocessed_dir_name + '.txt', 'r') as map_file:
            map_lines = map_file.readlines()

        # paste lines together but filter out zeros
        output_line_list = list()
        for mapping, weight in zip(map_lines[1:], temp_lines):      # skip header
            if weight.strip() == "0.00000000000000000e+00":
                continue
            feature_path = mapping.strip().split('\t')[-1]          # <- index gone
            output_line_list.append(f"{feature_path}\t{weight}")
        
        # ---- intercept ----------------------------------------------------
        xml_path = (
            f"{preprocessed_dir_name}{self.get_lambda_tag()}_out_feature_weights.xml"
        )
        with open(xml_path) as xf:
            xml_txt = xf.read()
        intercept_val = float(
            re.search(r"<intercept_value>(.*?)</intercept_value>", xml_txt).group(1)
        )
        #add intercept as last line
        output_line_list.append(f"Intercept\t{intercept_val:.15e}\n")

        # ---- write file & clean up ----------------------------------------
        with open(self.output, 'w') as oh:
            oh.writelines(output_line_list)

        os.remove('temp_out_feature_weights.txt') # temp weights
        os.remove(xml_path) # XML now redundant
        return

    def calc_preds_and_weights(self, only_pos_gss = False,
                                     skip_pred = False, sites = False):
        '''tally gene GSSs and species prediction scores'''
        # create dict for tracking GSSs; keys: gene objects, values: GSSs 
        included_gene_weights = defaultdict(lambda : 0) 
        # species prediction scores
        self.species_scores = defaultdict(lambda: 0.0)

        ###### tally all selected sites ######
        # get lines of esl model feature weights from the lasso output txt file
        weights_file = open(self.output,'r') # model weights for sites
        for position_line in weights_file:
            label, weight = parse_ESL_weight_line(position_line)
            # --- intercept line -------------------------------------------
            if label == "Intercept":
                self.y_intercept = weight # remember for later
                continue # skip normal processing
            if weight == 0.0:
                continue # skip zero weights which will be most of them
            gene_name, position, aa_to_check_for = parse_ESL_weight_label(label)

            ### first adjust the gene's weight sum (tabulating the GSS) ###
            # get gene obj
            current_gene_obj = self.run_family.gene_objects_dict[gene_name] 
            if only_pos_gss:
                if weight > 0:
                    included_gene_weights[current_gene_obj] += weight
            else: # normally this.  only_pos_gss is False by default
                included_gene_weights[current_gene_obj] += abs(weight)

            ### add site to gene's sites ###
            if (current_gene_obj.selected_sites[position] < weight
                and weight> 0 and sites):
                current_gene_obj.selected_sites[position] = weight
            
            ### Predictions  ###
            # adjust species scores tally by checking sequences
            if not skip_pred:
                alignment_file = open(os.path.join(
                    self.run_family.args.prediction_alignments_dir,
                    gene_name + '.fas'))
                for line in alignment_file: # loop through alignment lines
                    if line[0] == '>': 
                        species = line.strip('>\n') # get species name
                    else: # all other lines are sequence lines
                        if species in self.run_family.input_pheno_dict:
                            continue # skip input species here
                        if line[position] == aa_to_check_for: # 0-indexed
                            self.species_scores[species] += weight #add the site
                alignment_file.close()

                # now get scores from input alignments to calculate input RMSE
                input_alignment_file = open(os.path.join(
                    self.run_family.input_alignments_dir,
                    gene_name + '.fas'))
                for line in input_alignment_file: # loop through alignment lines
                    if line[0] == '>': 
                        species = line.strip('>\n') # get species name
                    else: # all other lines are sequence lines
                        if species not in self.run_family.input_pheno_dict:
                            continue # only look at input species here
                        if line[position] == aa_to_check_for: # 0-indexed
                            self.species_scores[species] += weight #add the site
                input_alignment_file.close()
        weights_file.close()
        for sp in self.species_scores:
            self.species_scores[sp] += self.y_intercept #add intercept to all scores
        
        ###### calculate input species RMSE (MFS) ######
        if not skip_pred:
            sum_of_squared_diffs = 0 #sum of (predicted - observed)^2
            n_species = len(self.run_family.input_pheno_dict.keys())
            for species in self.run_family.input_pheno_dict.keys():
                predicted = self.run_family.input_pheno_dict[species]
                observed = self.species_scores[species]
                sum_of_squared_diffs += (predicted - observed) ** 2
            # calc RMSE for run
            self.input_rmse = math.sqrt(sum_of_squared_diffs/ n_species)

        ###### update scores for all included genes ######
        # get GSS ranks for run
        ranked_genes = sorted(included_gene_weights.items(),
                              key = lambda gene_gss : gene_gss[1],
                              reverse = True)
        self.num_included_genes = len(ranked_genes) # record num genes included
        for rank, gene_gss in enumerate(ranked_genes, 1):
            gene_obj, gss = gene_gss # gene_gss is tuple of gene object and GSS
            if not gene_obj.best_rank: # best_rank will initially be None
                gene_obj.best_rank = rank
            elif rank < gene_obj.best_rank: # update gene's best rank
                gene_obj.best_rank = rank
            if gss > gene_obj.highest_gss: # update gene's highest gss
                gene_obj.highest_gss = gss
        return  


class GeneObject():
    '''a class to track data about genes'''
    def __init__(self, name):
        self.name = name
        self.num_var_sites = 0
        self.length = 0
        
    def __str__(self):
        return self.name + '\t'

    def __repr__(self):
        return 'gene_object: ' + self.name

class ESLGeneObject(GeneObject):
    '''a class to track weights and ranks in ESL runs for each gene'''
    def __init__(self, name):
        GeneObject.__init__(self, name)
        self.selected_sites = defaultdict(lambda:0) #key=position, val=top score
        self.highest_gss = 0
        self.best_rank = None
        self.highest_ever_gss = 0 # for multi matrix runs
        self.best_ever_rank = None # for multi matrix runs
        self.num_combos_ranked = 0 # for multi matrix runs
        self.num_combos_ranked_top = 0 # times ranked above top_rank_threshold

class SiteCounterGene(GeneObject):
    '''a class to track data about each gene from results'''
    def __init__(self, name, sequence_records):
        GeneObject.__init__(self, name)
        self.sequence_records = sequence_records
        self.variable_sites = [] # num non singleton sites
        self.length = len(self.sequence_records[0])

class SiteObject():
    '''a class to track sites for site counting'''
    def __init__(self, position, num_alleles = None, converge_degree = None,
                 num_non_gap_pairs = None):
        self.position = position
        self.num_alleles = num_alleles
        self.converge_degree = converge_degree
        self.num_non_gap_pairs = num_non_gap_pairs
        self.conv_score = converge_degree / (num_non_gap_pairs * 2)

    def __str__(self):
        return self.name + '\t'

    def __repr__(self):
        return 'SiteObject: ' + str(self.position)

class ESLGeneDict(dict):
    '''a dictionary of ESLGeneObjects for ESL integration runs'''
    def __init__(self, list_of_gene_names):
        if list_of_gene_names:
            for gene_name in list_of_gene_names:
                # add a new gene_object to self dictionary for each gene in list
                self[gene_name] = ESLGeneObject(gene_name)

    def get_sorted_list(self, multimatrix = False):
        '''returns a list of the gene objects sorted by best rank and then
        highest GSS. if multimatrix is true then it will sort by
        best_ever_rank and highest_ever_gss.  It sorts such that all the
        genes where rank is None are at the bottom. if multimatrix is True
        the genes will be sorted by num combos ranked and top ranked after this.
        '''
        # sort the list of gene objects by best rank and highest GSS
        #   the key tuple is sorted item by item, the first item is inf if
        #   the bestrank is None, and inf > any int so those genes go to the
        #   end, then for ties in best rank, it sorts by reverse (negative)
        #   of highest GSS
        rank_var = 'best_rank' if not multimatrix else 'best_ever_rank'
        gss_var = 'highest_gss' if not multimatrix else 'highest_ever_gss'
        sorted_genes = sorted(list(self.values()),key = lambda gene:
                              (float('inf') if vars(gene)[rank_var] is None
                    else vars(gene)[rank_var],
                    - vars(gene)[gss_var])) # sort by - of highest_gss 2nd
        if multimatrix:
            sorted_genes = sorted(sorted_genes, key = lambda gene :
                          (gene.num_combos_ranked, gene.num_combos_ranked_top),
                          reverse = True)
            
        return sorted_genes
