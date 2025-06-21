# a unified script for running integrated ESL analyses

import argparse, os, time
from . import esl_psc_functions as ecf
from collections import defaultdict


def get_esl_args(parser = None):
    '''this function takes an optional premade argparser object and then adds
    these standard ESL-PSC arguments and returns a the parser
    object.'''
    # Create the parser if an existing one hasn't been passed in
    if not parser:
        parser = argparse.ArgumentParser(description = 'ESL-PSC integrator')
    # ************Add arguments************

    ######### Directories and Paths #########
    group = parser.add_argument_group('Directories and Paths')

    help_txt = '''The full path to the folder where the preprocessed input
    folder will be. this same folder will be where the output feature weights
    files are put. 
    '''
    group.add_argument('--esl_inputs_outputs_dir', help = help_txt,
                           type = str)
    help_txt = '''The full path to the species phenotypes file which has the 
    species name then a comma and then a 1 or -1 for the phenotype class.
    Any species that is not in the phenotype file will not be included in the
    predictions output even if it was in the prediction alignments. If the
    phenotype file is not included, all species in the alignments will get
    SPSs but no true phenotypes will be listed.
    '''
    group.add_argument('--species_pheno_path', help = help_txt,
                           type = str, default = None)
    help_txt = '''The full path to the directory that contains the full
    alignments that include the species not in the response matrix for which
    predictions will be generated. Required for prediction output.
    '''
    group.add_argument('--prediction_alignments_dir', help = help_txt,
                           type = str, required = False)
    help_txt = '''* The name of this overall run to use in the output files'''
    group.add_argument('--output_file_base_name', help = help_txt,
                           type = str, required = True)
    help_txt = '''where the output will be saved. default = ESL-PSC directory'''
    group.add_argument('--output_dir', help = help_txt, type = str)

    ######### Hyperparameters #########
    group = parser.add_argument_group('Hyperparameters')
    help_txt = '''initial lambda 1 value (lambda1 is the position sparsity
    parameter) default = .01
    '''
    group.add_argument('--initial_lambda1', help = help_txt,
                        type = float, default = .01)
    help_txt = '''final lambda 1 value (lambda1 is the position sparsity
    parameter) default = .99
    '''
    group.add_argument('--final_lambda1', help = help_txt,
                        type = float, default = .99)
    help_txt = '''initial lambda 2 value (lambda2 is the group sparsity
    parameter) default = .01
    '''
    group.add_argument('--initial_lambda2', help = help_txt,
                        type = float, default = .01)
    help_txt = '''final lambda 2 value (lambda2 is the group sparsity
    parameter) default = .99
    '''
    group.add_argument('--final_lambda2', help = help_txt,
                        type = float, default = .99)
    help_txt = '''the increment to increase the lambda values with each step'''
    group.add_argument('--lambda_step', help = help_txt,
                        type = float, default = .05)
    help_txt = '''group penalty calculation type ("sqrt", "median", and "linear"
    are available) note that if "sqrt" is chosen there will be no group penalty
    range. We recommend "median" where a single linear penalty with the median
    number of variable sites across all alignments with at least one variable
    site as the constant term will be used. In this case, any other settings
    for the group penalty terms (final, initial, step) will be ignored. The
    type can also be "std" which results in the lasso default sqrt penalty.
    Beware when reusing existing preprocess in which the group penalties were
    rewritten by a previous run with a different penalty calculation type. In
    this case, the preprocess must be repeated to get std penalties back.
    '''
    group.add_argument('--group_penalty_type', help = help_txt,
                        type = str, default = "median")
    help_txt = '''group penalty constant term initial value'''
    group.add_argument('--initial_gp_value', help = help_txt,
                        type = int, default = 1)
    help_txt = '''group penalty constant term final value'''
    group.add_argument('--final_gp_value', help = help_txt,
                        type = int, default = 1)
    help_txt = '''group penalty constant term increment'''
    group.add_argument('--gp_step', help = help_txt, type=int, default = 6)
    help_txt = '''the number of points per lambda parameter in a logspace
    of values to test. This number squared will be the total grid points.'''
    group.add_argument('--num_log_points', help = help_txt,
                        type = int, default = 20)
    help_txt = '''the names of the two phenotypes seperated by a space with
    the convergent phenotype coming first.'''
    group.add_argument('--pheno_names',
                        nargs = 2,
                        metavar = ('convergent_phenotype',
                                   'nonconvergent phenotype'),
                        help = help_txt,
                        type = str)
    help_txt = '''minimum number of genes a model must have in order for that
    model to be included in the prediction scores plots. default = 0'''
    group.add_argument('--min_genes', help = help_txt,
                        type = int, default = 0)
   

    ######### Options #########
    group = parser.add_argument_group('Options')
    help_txt = '''Instead of using initial and final lambda values with a lambda
    step to generate a linear grid of points, this option will use a log space
    of points. The number of points to use can be given as a hyperparameter.
    the same defaults are used for the initial and final lambda values.
    '''
    group.add_argument('--use_logspace', help = help_txt,
                        action = 'store_true', default = False)
    help_txt = '''don't run preprocess, assume the folder is there already
    and use it
    '''
    group.add_argument('--use_existing_preprocess', help = help_txt,
                        action = 'store_true', default = False)
    help_txt = '''don't replace group penalties. This is automatically set to
    True if the group_penalty_type is "std"'''
    group.add_argument('--use_default_gp', help = help_txt,
                        action = 'store_true', default = False)
    help_txt = '''don't delete the raw model output files for each run'''
    group.add_argument('--keep_raw_output', help = help_txt,
                        action = 'store_true', default = False)
    help_txt = '''print a dictionary of all selected sites with their highest
    model score for every gene in the gene_ranks output file.'''
    group.add_argument('--show_selected_sites', help = help_txt,
                        action = 'store_true', default = False)
    help_txt = ''' don't output a gene ranks file'''
    group.add_argument('--no_genes_output', help = help_txt,
                        action = 'store_true', default = False)
    group2 = group.add_mutually_exclusive_group() #can't make plot with no pred
    help_txt = '''don't output a species predictions file'''
    group2.add_argument('--no_pred_output', help = help_txt,
                        action = 'store_true', default = False)
    help_txt = '''make a violin plot showing sps density for each true
    phenotype. By default this will display SPS of > 1 or < -1 as 1 amd -1'''
    group2.add_argument('--make_sps_plot', help = help_txt,
                        action = 'store_true', default = False)
    help_txt = '''make a kde plot showing sps density for each true
    phenotype.'''
    group2.add_argument('--make_sps_kde_plot', help = help_txt,
                        action = 'store_true', default = False)
    
    
    # return the parser, the args still need to be parsed
    return parser

def get_single_run_esl_args(parser, input_alignments_req = True):
    '''this function takes a premade argparser group object and then adds
    these esl-psc arguments for one esl run. In, e.g., a multi matrix run,
    these won't be necesary.
    '''
    group = parser.add_argument_group('ESL Input Files')
    help_txt = '''*the full path to the response matrix file'''
    group.add_argument('--response_matrix_path', help = help_txt,
                       type = str, required = True)
    help_txt = '''the full path to the alignment path file. if not given
    one will be generated.
    '''
    group.add_argument('--path_file_path', help = help_txt,
                       type = str, default = None)
    help_txt = '''the name (not full path) of the folder of preprocessed
    input for lasso. if not given one will be generated
    '''
    group.add_argument('--preprocessed_dir_name', help = help_txt,
                       type = str, default = None)
    help_txt = '''the full path to the directory of alignments with of only
    input species'''
    if input_alignments_req: 
        help_txt = '*' + help_txt #add a * to indicate required
    group.add_argument('--input_alignments_dir', help = help_txt,
                       type = str, required = input_alignments_req)


def replace_group_penalties(esl_inputs_outputs_dir, gene_objects_dict,
                            penalty_function, input_species_list,
                            preprocessed_dir_name, input_alignments_dir):
    ''' takes esl args and a penalty function and modifies the group indices
    file in the preprocessed input to change the group penalties.
    the penalty function operates on the number of variable sites.
    '''
    # go to alignment folder for this species combo
    os.chdir(input_alignments_dir)
    gene_list = list(gene_objects_dict.keys())
    # loop through alignment files in order
    new_penalties = [] # list of numbers
    alignment_file_list = [name + '.fas' for name in gene_list]
    for file_path in alignment_file_list:
        records = ecf.get_seq_records_in_order(file_path,
                                               input_species_list)  
        # calculate new penalty based on num variable sites and add to list
        new_penalties.append(penalty_function(ecf.count_var_sites(records)))
    # now modify penalties in group index file
    group_indices_file = (preprocessed_dir_name + '/group_indices_'
                               + preprocessed_dir_name + '.txt')
    group_indices_file = os.path.join(esl_inputs_outputs_dir,
                                      group_indices_file)
    with open(group_indices_file, 'r') as file:
        group_indices_lines = file.readlines() # get existing lines
    # now just replace the old penalties
    group_indices_lines[2] = ','.join([str(penalty)
                                        for penalty in new_penalties]) 
    with open(group_indices_file, 'w') as file:
        file.write(''.join(group_indices_lines)) #overwrite new penalties
    return



def esl_integration(args,
                    input_species_list,
                    preprocessed_dir_name,
                    input_alignments_dir,
                    gene_objects_dict,
                    label = ''):
    """runs ESL for given inputs across lambda parameters and group penalties
    and returns dictionaries of objects for runs and genes. gene_objects_dict is
    a dictionary of ESLGeneObject as values and gene names as keys. label is
    used as a name for the integration, e.g. a species combo code for a
    multimatrix integration.
    """
    
    # check for errors
    assert(not os.path.isabs(preprocessed_dir_name)) # not full path!
    assert(os.path.exists(args.esl_inputs_outputs_dir)) # check that this exists
    ###

    # prepare variables for this integrated set of runs
    
    # dictionary of input species phenotypes
    input_pheno_dict = ecf.make_input_pheno_dict(input_species_list)
    print('species being checked:\n', ", ".join(input_species_list))

    # create a list of esl_run objects for the whole integration
    esl_run_list = []

    # set some args
    args.only_pos_gss = False
    args.lambda1_only = False

    # if using the median penalty, modify args for penalty term accordingly
    if args.group_penalty_type == 'median':
        median_gp = ecf.get_median_var_sites(input_alignments_dir)
        initial_penalty = median_gp
        final_penalty = median_gp
        penalty_type = 'linear'
    else:
        initial_penalty = args.initial_gp_value
        final_penalty = args.final_gp_value
        penalty_type = args.group_penalty_type

    ############## Loop through penalty terms ##############
    penalty_term = initial_penalty
    while penalty_term <= final_penalty:
        print('Penalty term: ' + str(penalty_term))

        
        ########## Replace Group Penalties ##########
        # count variable sites in each alignment file and replace group
        # penalties in the group indices file
        if (args.use_default_gp or (args.group_penalty_type == "std")):
            pass # skip replacing penalties and use std sqrt penalties 
        else:
            # make penalty function which will take num var sites as its arg
            penalty_function = ecf.penalty_function_maker(penalty_term,
                                                          penalty_type)
            replace_group_penalties(args.esl_inputs_outputs_dir,
                                    gene_objects_dict,
                                    penalty_function,
                                    input_species_list,
                                    preprocessed_dir_name,
                                    input_alignments_dir)

        ########## Run ESL grid search ##########
        os.chdir(args.esl_inputs_outputs_dir)
        # make a run family object for this grid search
        grid_run_family = ecf.ESLRunFamily(args,
                                           input_pheno_dict,
                                           preprocessed_dir_name,
                                           penalty_term,
                                           input_alignments_dir,
                                           label,
                                           gene_objects_dict =
                                           gene_objects_dict)
        # run all esl runs with all parameters generating all out put files
        grid_run_family.gridsearch() # gets needed parameters from args
        
        # do predictions and update
        grid_run_family.do_all_calculations()
        
        # add runs to full run_list
        esl_run_list.extend(grid_run_family.runs_list) 
                                         
        # remove out_feature_weights.txt files to not contaminate next iteration
        if not args.keep_raw_output: # skip if this option is given
            for file_name in os.listdir():
                if file_name[-3:] == 'txt' or file_name[-3:] == 'xml':
                    os.remove(file_name)

        # End of while loop to loop through penalty terms; increment term
        penalty_term += args.gp_step   
    
    # end of esl_integration() fuinction 
    return gene_objects_dict, esl_run_list

def generate_gene_ranks_output(gene_objects_dict, output_dir, output_file_name,
                               show_sites = False, multimatrix = False):
    """generate a database of genes with integrated ESL ranks, max GSS, etc.
    with an option to show the selected sites for every gene in a dict format
    with positions as keys and highest weight for that position (any residue)
    """
    print("\nGenerating gene ranks output file...")
    # set the data to display depending on whether or not it's for multimatrix 
    if not multimatrix: # for single matrix integration runs output
        column_names = ['gene_name', 'highest_gss', 'best_rank']
        def get_line(gene): 
            return [gene.name, str(gene.highest_gss), str(gene.best_rank)]
    else: # for multimatrix output
        column_names = ['gene_name',
                    'num_combos_ranked',
                    'num_combos_ranked_top',
                    'highest_ever_gss',
                    'best_ever_rank']
        def get_line(gene):
            return [gene.name,
                    str(gene.num_combos_ranked),
                    str(gene.num_combos_ranked_top),
                    str(gene.highest_ever_gss),
                    str(gene.best_ever_rank)]
    if show_sites:
        column_names.extend(['num_selected_sites',
                             'sites_and_weights'])
    # sort the list of gene objects, if multimatrix it will sort accordingly
    sorted_genes = gene_objects_dict.get_sorted_list(multimatrix = multimatrix)
        
    # loop through and add data
    headers = ','.join(column_names)
    output_data_lines = [headers] # keep data lines in a list. start w headers
    for gene in sorted_genes:
        line = get_line(gene) # this is defined above for single or multimatrix
        if show_sites: # if true, add the number of sites and the sites dict
            line.extend([str(len(gene.selected_sites)),
                         str(dict(gene.selected_sites))]) 
        output_data_lines.append(','.join(line))
    # make output file path
    output_path = os.path.join(output_dir,
                               output_file_name + '_gene_ranks.csv')
    # write output
    with open(output_path, 'w') as file:
        file.write('\n'.join(output_data_lines))
    print("Gene ranks file written as: " + output_path)
    return

def generate_predictions_output(esl_run_list, output_path, phenofile = None):
    """generate a database output file with all species predictions from run
    list data. output csv file will be saved in output_dir
    """
    print("Generating predictions output file...")
    column_names = ['species_combo', 'lambda1','lambda2', 'penalty_term',
                    'num_genes', 'input_RMSE', 'species', 'SPS']
    # if an all_species_phenotype file has been included then read that
    if phenofile: 
        column_names.append('true_phenotype')
        # read the file
        all_species_pheno_dict = ecf.get_pheno_dict(phenofile)
    else:
        # this will just let it add nothing to each line if no pheno file 
        all_species_pheno_dict = defaultdict(lambda:'') 

    #loop through runs and add lines of data
    headers = ','.join(column_names)
    output_data_lines = [headers] # keep data lines in a list. start w headers
    for run in esl_run_list:
        for species, score in run.species_scores.items():
            if (species in run.run_family.input_pheno_dict or 
                (phenofile and species not in all_species_pheno_dict)):
                # don't include predictions of species in the input set
                # also, if there is a phenotype file (phenofile != None)
                # then if the species is not in the pheno file exclude it
                continue 
            true_pheno = str(all_species_pheno_dict[species]) # get true pheno
            # construct output line
            line = [run.run_family.species_combo_tag, # species used
                    str(run.lambda1), # lambda 1
                    str(run.lambda2), # lambda 2
                    str(run.run_family.penalty_term), # group penalty term
                    str(run.num_included_genes), #num genes select
                    str(run.input_rmse),
                    species, # species name
                    str(score)] #prediction score (SPS)
            if phenofile:
                line.append(true_pheno)
            output_data_lines.append(','.join(line))
    
    # write output
    with open(output_path, 'w') as file:
        file.write('\n'.join(output_data_lines))
    print("Predictions file written as: " + output_path)
    return output_path


if __name__ == '__main__':
    start_time = time.time()

    # get parser
    parser = get_esl_args()

    # get paths to response file, path file, and input alignment directory
    get_single_run_esl_args(parser)  

    args = ecf.parse_args_with_config(parser)

    # set output_dir
    if not args.output_dir:
        args.output_dir = args.esl_main_dir

    # set group penalty std if necessary
    if args.group_penalty_type == "std":
        args.use_default_gp = True

    # generate preprocess name and path file if they were not given in args
    if not args.path_file_path:
        args.path_file_path = ecf.make_path_file(args.input_alignments_dir)

    if not args.preprocessed_dir_name:
        # construct name from response matrix and alignments name
        args.preprocessed_dir_name = os.path.join(args.esl_inputs_outputs_dir,
                     os.path.basename(args.input_alignments_dir)
                     + "_" + os.path.basename(args.response_matrix_path))
        
    
    # list of input species
    input_species_list = ecf.get_species_to_check(args.response_matrix_path)
    
    # get a list of all gene names from alignment file names in the path file
    gene_list = ecf.get_gene_names(args.path_file_path)
    # make initial blank gene_objects_dict
    gene_objects_dict = ecf.ESLGeneDict(gene_list)

    # clear preexisting output files in the inputs folder
    os.chdir(args.esl_inputs_outputs_dir)
    for file in os.listdir():
        if file[-4:] == ".txt" or file[-4:] == ".xml":
            os.remove(file)

    # Run preprocess if necessary
    #   get full absolute path to preproces directory
    preprocess_dir_full_path = os.path.join(args.esl_inputs_outputs_dir,
                                            args.preprocessed_dir_name)
    if not args.use_existing_preprocess:
        # if the folder of preprocess files already exists, delete it
        ecf.clear_existing_folder(preprocess_dir_full_path)

        # run ESL preprocess once
        ecf.run_preprocess(args.esl_main_dir,
                           args.response_matrix_path,
                           args.path_file_path,
                           args.preprocessed_dir_name,
                           args.esl_inputs_outputs_dir,
                           use_is = True)
    else: # if we are skipping the preprocess step, the folder should exist
        assert os.path.exists(preprocess_dir_full_path)    

    # Call the main function which should return dictionaries of runs and genes
    # assign its output to these two variables
    gene_objects_dict, esl_run_list = esl_integration(args,
                                                    input_species_list,
                                                    args.preprocessed_dir_name,
                                                    args.input_alignments_dir,
                                                    gene_objects_dict)
    
    
    # Print a clear, multi-line summary so the GUI output isn't squished
    print("\nESL-PSC integration finished!")
    print(f"A total of {len(esl_run_list)} ESL models were built")
    print("The arguments for this integration run were:")
    for key, value in vars(args).items():
        print(f"{key} = {value}")
    
    # call output functions which should generate output text files
    if not args.no_genes_output: # skip this output if flag is true
        generate_gene_ranks_output(gene_objects_dict, args.esl_main_dir,
                                   args.output_file_base_name,
                                   show_sites = args.show_selected_sites)
    print('\n')
    if not args.no_pred_output: # skip this output if flag is true
        # make full file path of output predictions file
        preds_output_path = os.path.join(args.esl_main_dir,
                                         args.output_file_base_name
                                         + '_species_predictions.csv')
        generate_predictions_output(esl_run_list,
                                    preds_output_path,
                                    args.species_pheno_path)
        
    ecf.report_elapsed_time(start_time) # print time before making plot

    if args.make_sps_plot or args.make_sps_kde_plot:
        plot_type = 'violin' if args.make_sps_plot else 'kde'
        # generate and show density plots of predictions
        # call sps_density.create_sps_plots for various rmse cutoffs
        ecf.rmse_range_pred_plots(preds_output_path,
                                  args.output_file_base_name,
                                  args.pheno_names,
                                  args.min_genes,
                                  plot_type)

