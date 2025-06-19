# A script to automate ESL-PSC integration experiments for multiple input
#  matrices of species

import argparse, os, time, shutil, random, itertools
from . import esl_integrator as esl_int
from . import esl_psc_functions as ecf
from . import deletion_canceler as dc
from Bio import SeqIO
from Bio.Seq import Seq

def validate_specific_paths(args):
    '''
    Check that certain known arguments (if provided) are valid paths.
    '''
    paths_to_check = [
        'esl_main_dir',
        'alignments_dir',
        'prediction_alignments_dir',
        'response_dir',
        'species_groups_file',
        'species_pheno_path',
        'canceled_alignments_dir'
        # Add or remove as needed for your script
    ]

    problem_paths = []

    for attr_name in paths_to_check:
        path_value = getattr(args, attr_name, None)
        # If user didn't supply it, or it's None, skip it
        if path_value is None:
            continue

        # Check existence
        if not os.path.exists(path_value):
            problem_paths.append(f"{attr_name} = '{path_value}'")

    if problem_paths:
        error_msg = (
            "The following paths were provided but do not exist on disk:\n"
            + "\n".join(problem_paths)
            + "\nPlease check these paths for errors."
        )
        raise ValueError(error_msg)

def check_species_in_alignments(list_of_species_combos, alignments_dir):
    '''
    Checks that every species in each combo from the species_groups_file
    is actually present in at least one alignment file in alignments_dir.
    Raises a ValueError if any species is missing.
    '''
    # Collect every species from the combos into a set
    all_species = set()
    for combo in list_of_species_combos:
        all_species.update(combo)
    
    # Collect all species actually found in the alignment files
    found_species = set()
    if not os.path.isdir(alignments_dir):
        raise ValueError(f"Alignment directory '{alignments_dir}'"
                         "does not exist or is not a directory.")
    
    for file_name in os.listdir(alignments_dir):
        if ecf.is_fasta(file_name):
            file_path = os.path.join(alignments_dir, file_name)
            for record in SeqIO.parse(file_path, 'fasta'):
                # record.id is your species ID/header
                found_species.add(record.id)
    
    # Determine if any species are missing
    missing_species = all_species - found_species
    if missing_species:
        raise ValueError(
            f"The following species from the species groups file were not "
            f"found in any alignment in '{alignments_dir}':\n"
            f"{', '.join(sorted(missing_species))} Double check the spelling"
        )

def randomize_alignments(original_alignments_directory, species_list):
    '''generates randomly pair-flipped alignments from deletion-canceled
    input alignment files
    '''
    # make new dir path, which will always be called the same thing
    parent_directory, _ = os.path.split(original_alignments_directory)
    scrambled_alignments_dir = os.path.join(parent_directory,
                                            "scrambled_alignments/")
    # clear existing directory and create new one
    ecf.clear_existing_folder(scrambled_alignments_dir)
    os.mkdir(scrambled_alignments_dir)
    
    # loop through alignments in original_alignments_directory and randomize
    for file in os.listdir(original_alignments_directory):
        if not ecf.is_fasta(file):
            continue
##        if file_count % 1000 == 0:
##            print('scanning file number ' + str(file_count))
        os.chdir(original_alignments_directory) #be in original files directory
        
        # get seq records in order of the species combo i.e. 1, -1, 1, -1 etc.
        records = ecf.get_seq_records_in_order(file, species_list)

        # ***Now do the scrambling
        # make a list of the sequences which are converted from str to lists
        # so now we have a list of lists of AAs
        seq_list = [list(record.seq) for record in records]

        # loop through all positions, scramble AAs
        for index in range(len(seq_list[0])):
            # make list of AAs at this position
            position_list = [seq[index] for seq in seq_list]

            # scramble the list
            # split into list of lists with one pair per sublist
            pair_list = [position_list[n:n+2] for n in
                         range(0,len(position_list),2)]
            
            # loop through pairs and shuffle them
            for pair_index in range(len(pair_list)):
                random.shuffle(pair_list[pair_index])
                
            # flatten the list back out
            position_list = list(itertools.chain(*pair_list))

            # change each sequence to the AA from the scrambled position list
            for seq_num in range(len(seq_list)):
                seq_list[seq_num][index] = position_list[seq_num]

        # now write new fasta file with modified sequences
        # first modify seq_records
        for seq_num, record in enumerate(records):
            record.seq = Seq(''.join(seq_list[seq_num]))
        # change to new alignment files directory to write new file
        os.chdir(scrambled_alignments_dir)
        
        # write new file
        with open(file, "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta-2line")
    
    # return path to new directory
    return scrambled_alignments_dir
    

def run_multi_matrix_integration(args, list_of_species_combos,
                                 response_file_list):
    '''Run esl_integrator for each of a list of species combos with
    existing response matrix files in response_file_list.
    species_combo_list is a list of tuples that each has one species combo
    response_file_list is a list of full paths to response matrix files in the
        same order at the corresponding species combos in species_combo_list
    '''
    # Make a gene_objects_dict
    # look in one of the canceled alignemnt directories
    if len(list_of_species_combos) > 1 and not args.use_uncanceled_alignments:
        alignment_sub_dir = os.listdir(args.canceled_alignments_dir)[0]
        alignment_sub_dir = os.path.join(args.canceled_alignments_dir,
                                         alignment_sub_dir)
    else: # if its just one matrix being run with multimatrix (or uncanceled)
        alignment_sub_dir = args.canceled_alignments_dir
    gene_name_list = ecf.get_gene_names(alignment_sub_dir)
    # make gene_objects_dict from name list
    gene_objects_dict = ecf.ESLGeneDict(gene_name_list)

    # set top_rank_threshold (rank abovewhich to count genes as top genes)
    top_rank_threshold = max(1, len(gene_objects_dict) * args.top_rank_frac)

    master_run_list = [] # a list of all runs from all matrices

    # loop through combos
    for combo_num, response_file in enumerate(response_file_list):
        combo = list_of_species_combos[combo_num]
        combo_name = 'combo_' + str(combo_num)
        total_combos = len(list_of_species_combos)
        
        # MODIFICATION: Announce combo processing in a parsable way
        print(f"\n--- Processing combo {combo_num + 1} of {total_combos} ({combo_name}) ---")
        
        response_path = response_file_list[combo_num]

        # get name of one alignment directory in the canceled_alignments_dir
        if (len(list_of_species_combos) > 1
            and not args.use_uncanceled_alignments): 
            # use combo name
            gap_canceled_alignments_path = (
                os.path.join(args.canceled_alignments_dir,
                             combo_name + '-alignments'))
        else: # if its just one matrix being run or uncanceled alignments
            #in this case the canceled_alignments_dir has the files in itself
            gap_canceled_alignments_path = args.canceled_alignments_dir
        # generate a path file
        path_file_path = ecf.make_path_file(gap_canceled_alignments_path)

        # name of preprocess directory (not full path)
        preprocess_dir_name = args.output_file_base_name + '_' + combo_name
        
        if not args.make_pair_randomized_null_models:
            # ***Do a Normal Multimatrix Integration***
            # run preprocess if needed
            if not args.use_existing_preprocess:
                # if the folder is there already remove it first
                ecf.clear_existing_folder(
                    os.path.join(args.esl_inputs_outputs_dir,
                                 preprocess_dir_name))
                ecf.run_preprocess(args.esl_main_dir, response_path,
                                   path_file_path, preprocess_dir_name,
                                   args.esl_inputs_outputs_dir, use_is = True)
                
            # run esl integration
            _, run_list = esl_int.esl_integration(args,
                                                  combo,
                                                  preprocess_dir_name,
                                                  gap_canceled_alignments_path,
                                                  gene_objects_dict,
                                                  combo_name)
            # gene_objects_dict is an object and only a reference is passed to
            #   each run so same object persists and accumulates all the data
            master_run_list.extend(run_list)
        else:
            # ***Do a Randomized Alignment Null Multimatrix Integration***
            for run_num in range(args.num_randomized_alignments):
                # first randomize the alignment and get path to new align dir
                rand_aligns = randomize_alignments(gap_canceled_alignments_path,
                                                   combo) # species list 
                path_file_path = ecf.make_path_file(rand_aligns)
                # then re-run preprocess
                ecf.clear_existing_folder(
                    os.path.join(args.esl_inputs_outputs_dir,
                                 preprocess_dir_name))
                ecf.run_preprocess(args.esl_main_dir, response_path,
                                   path_file_path, preprocess_dir_name,
                                   args.esl_inputs_outputs_dir, use_is = True)
                # then repeat integration
                _, run_list = esl_int.esl_integration(args,
                                                  combo,
                                                  preprocess_dir_name,
                                                  rand_aligns,
                                                  gene_objects_dict,
                                                  combo_name
                                                      + '_' + str(run_num))
                master_run_list.extend(run_list)
        

        # update gene variables to track best scores and num combos ranked etc.
        for gene in gene_objects_dict.values():
            if gene.best_rank: # count if ranked at all
                gene.num_combos_ranked += 1
            else:
                continue # if never ranked the gene won't have any GSS 
            if gene.best_rank <= top_rank_threshold: # count if top gene
                gene.num_combos_ranked_top += 1              
            # best gene rank
            if not gene.best_ever_rank: 
                gene.best_ever_rank = gene.best_rank # set if still None
            elif gene.best_rank < gene.best_ever_rank: # check for best rank
                gene.best_ever_rank = gene.best_rank
            # highest gss
            if gene.highest_gss > gene.highest_ever_gss: # check for highest gss
                gene.highest_ever_gss = gene.highest_gss
            # reset variables for each individual combo
            gene.highest_gss = 0
            gene.best_rank = None
                        
        # delete path file and preprocess unless --preserve_preprocess
        os.remove(path_file_path)
        if args.delete_preprocess: # delete preprocess to keep folder clean
            shutil.rmtree(os.path.join(args.esl_inputs_outputs_dir,
                                       preprocess_dir_name))

    return gene_objects_dict, master_run_list


def main(raw_args=None):
    start_time = time.time()

    desc_text = '''This will run ESL integrations for many species combinations.
                All necessary args for esl_integrator.py must be included to
                specify how each integration run will be performed. Alignments
                must be in 2-line fasta format and file names end in ".fas".
                If no species groups file is given, existing response matrices
                must be given.  An * indicates required arguments. args can be
                given in a config file called esl_ct_config.txt with 1 per line.
                '''
    parser = argparse.ArgumentParser(description = desc_text)

    esl_int.get_esl_args(parser)

    dc.get_deletion_canceler_args(parser)

    #### Get args specific to multimatrix esl ####
    group = parser.add_argument_group('Multimatrix-specific Optional Arguments')
    group.add_argument('--alignments_dir',
                    help = ('Full path to the original alignments directory. '
                            'if this is not provided, the '
                            'prediction_alignments_dir will be used instead'),
                    type = str)
    group.add_argument('--top_rank_frac',
                        help = 'Fraction of genes to count as "top genes"',
                        type = float, default = .01)
    
    # Options 
    group.add_argument('--response_dir',
                        help = ('Folder with response matrices. Any txt file in'
                                'this folder is assumed to be a response '
                                'matrix file'),
                        type = str)
    group.add_argument('--use_uncanceled_alignments',
                        help = ('Use the alignments_dir alignments '
                                'for all matrices without doing gap canceling'),
                        action = 'store_true', default = False)
    group.add_argument('--use_existing_alignments',
                        help = 'Use existing files in canceled_alignments_dir ',
                        action = 'store_true', default = False)
    group.add_argument('--delete_preprocess',
                        help = 'Clear preprocess folders after each matrix run',
                        action = 'store_true', default = False)
    group.add_argument('--make_null_models',
                        help = ('Make null response-flipped ESL-PSC models. '
                                'must have an even number of pairs. All '
                                'balanced flippings of the responce values will'
                                ' be generated for each combo and all will be '
                                'run and aggregated to maximally decouple true '
                                'convergences'),
                        action = 'store_true', default = False)
    group.add_argument('--make_pair_randomized_null_models',
                        help = ('Make null pair randomized ESL-PSC models. '
                                'a copy of input deletion-canceled alignment '
                                'will, for each variable site, be randomized '
                                'such that the residues of each contrast pair '
                                'will be either flipped or not and the '
                                'esl integration will be repeated for each '
                                'one. the results are then aggregated for all'),
                        action = 'store_true', default = False)
    group.add_argument('--num_randomized_alignments',
                        help = 'number of pair-randomized alignments to make',
                        type = int, default = 10)
    
    # Ensure we have sensible, project-root defaults for these two paths
    args = ecf.parse_args_with_config(parser, raw_args) # checks for config file
    
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

    # 1) Where all intermediate ESL artefacts live
    if not getattr(args, "esl_inputs_outputs_dir", None):
        args.esl_inputs_outputs_dir = os.path.join(project_root,
                                                   "preprocessed_data_and_outputs")

    # 2) Root of the ESL install (used when spawning helper binaries);
    #    only override if the user did *not* supply a value.
    if not getattr(args, "esl_main_dir", None):
        args.esl_main_dir = project_root

    # Guarantee this folder exist
    os.makedirs(args.esl_inputs_outputs_dir, exist_ok=True)

    if not args.species_groups_file and not args.response_dir:
        # if 1 response file is given instead of species groups it won't work
        raise argparse.ArgumentTypeError('must include a species groups file '
                                         'or a response_dir with matrices')

    if args.alignments_dir and args.prediction_alignments_dir:
        # Both arguments have been provided, no action is required
        pass
    elif args.alignments_dir:
        # Only alignments_dir has been provided, assign
        # prediction_alignments_dir to the same value
        args.prediction_alignments_dir = args.alignments_dir
    elif args.prediction_alignments_dir:
        # Only prediction_alignments_dir has been provided, assign
        # alignments_dir to the same value
        args.alignments_dir = args.prediction_alignments_dir
    elif args.use_existing_alignments and args.no_pred_output:
        pass
    else:
        # Neither argument has been provided, raise an error
        raise ValueError("At least one of --alignments-dir or "
                         "--prediction-alignments-dir must be provided.")

    validate_specific_paths(args) #verify that dir and file paths are real

    # set output_dir
    if not args.output_dir:
        args.output_dir = args.esl_main_dir
    # ensure the output directory actually exists
    os.makedirs(args.output_dir, exist_ok=True)
    
    # if use_uncanceled_alignments option is given, set canceled alignments dir
    if args.use_uncanceled_alignments:
        args.canceled_alignments_dir = args.alignments_dir
        args.use_existing_alignments = True # set this for convenience
        
    
    # 1) Generate a List of Response Matrix Files. this can either be from a
    #   species group file or from a provided folder of response files

    # master list of species combinations
    #   a list of response matrix files serves as the main list of combos.
    #   they can be generated from the species groups file.
    #   The index of the file name in the list can be the combo code to
    #   link to the correct preprocess folder and gap-canceled alignments folder

    if args.response_dir: # this means there is a directory of response matrices
        # sort them to get deterministic ordering
        response_file_list = sorted(os.listdir(args.response_dir))
        response_file_list = [os.path.join(args.response_dir, file) for
                              file in response_file_list] # full paths
        response_dir = args.response_dir # for use later
        # generate list_of_species_combos
        list_of_species_combos = []
        for response_file in response_file_list:
            list_of_species_combos.append(
                ecf.get_species_to_check(response_file))
    elif args.species_groups_file: # a species group file was given
        response_file_list = []
        list_of_species_combos = dc.parse_species_groups(
            args.species_groups_file)
        # check that species names are all in alignments (correct spellings)
        check_species_in_alignments(list_of_species_combos, args.alignments_dir)
        # now make a directory of response files under dir with group file in it
#        group_file_parent_dir = os.path.split(args.species_groups_file)[0]
        response_dir = os.path.join(args.output_dir,
                                    os.path.split(args.species_groups_file)[1]
                                    + '_response_matrices')
        response_dir = response_dir.replace('.txt', '') #delete group file ext
        if not args.make_null_models: # this will be done later if doing nulls
            response_file_list = ecf.make_response_files(response_dir,
                                                         list_of_species_combos)
    else: 
        raise ValueError("must give either response_dir or species_groups_file")
    # we now have a response_file_list and a response_dir with the files in it

    if args.make_null_models:
        # if this option is chosen, we will take all of the species combos and
        #   and generate a new species combo list in which the pairs response
        #   values are flipped or not flipped in all possible ways
        list_of_species_combos = ecf.make_null_combos(list_of_species_combos)
        # now fix the response directory
        response_file_list = ecf.make_response_files(response_dir,
                                                     list_of_species_combos)
        

    # 2) Generate Gap-canceled Alignments
    if not args.use_existing_alignments: # skip if using existing alignments
        if not args.canceled_alignments_dir: # new alignments path not given
            # if no canceled_alignments_dir was given, generate a name for one
            aligns_parent_dir = os.path.split(args.alignments_dir)[0]
            if args.species_groups_file: # use species groups file name
                dir_name = args.species_groups_file.replace('.txt','')
            else: # use output base name if no groups file
                dir_name = args.output_file_base_name.replace('.txt','')
            dir_name += '_gap-canceled_alignments'
            args.canceled_alignments_dir = os.path.join(aligns_parent_dir,
                                                    dir_name)
        # if the folder already exists, remove it, but check with user first
        elif os.path.exists(args.canceled_alignments_dir):
            print('Canceled alignment directory: '
                  + args.canceled_alignments_dir)
            Q = input("the named gap-canceled alignments folder exists already."
                      " Are you sure you want to delete it and generate new "
                      "alignments? enter y to continue or anything else to "
                      "quit")
            if Q != 'y': # give the user a chance to quit
                raise Exception("Canceled by user. use "
                                "--use_existing_alignments option to use the "
                                "existing alignments but they must have "
                                "correct subfolder names")             
        ecf.clear_existing_folder(args.canceled_alignments_dir)
        # generate new alignments
        dc.generate_gap_canceled_alignments(args, list_of_species_combos,
                                            enumerate_combos = True,
                                            limited_genes_list =
                                            args.limited_genes_list)
    else: # if using existing alignments check to make sure a folder is there
        if not args.canceled_alignments_dir: 
            raise Exception("When the use_existing_alignments option is in "
                            "use, then a canceled_alignments_dir must be given")

    # 3) Loop Through the Response Files and Run Preprocess and Integration
    # clear preexisting output files in the inputs outputs folder folder
    previous_working_dir = os.getcwd() # record this
    os.chdir(args.esl_inputs_outputs_dir) # go to folder with raw output
    for file in os.listdir():
        if file[-4:] == ".txt" or file[-4:] == ".xml":
            os.remove(file) 
    os.chdir(previous_working_dir) # go back to original working dir

    # run multimatrix integration
    gene_objects_dict, master_run_list = run_multi_matrix_integration(
                                                    args,
                                                    list_of_species_combos,
                                                    response_file_list)
    # 4) generate output
    print("\nmultimatrix integration finished! ",
          "A total of " + str(len(master_run_list))
          + " ESL models were built\n",
          "The arguments for this integration run were:\n")
    for key, value in vars(args).items(): # repeat the input of args at the end
          print(str(key) + ' = ' + str(value))

    # print these paths so they don't get lost
    print('\nResponse matrices directory: ' + response_dir,
          '\nGap-canceled alignments directory: '
          + args.canceled_alignments_dir)
    
    # call output functions which should generate output files
    if not args.no_genes_output: # skip genes output if flag is true
        esl_int.generate_gene_ranks_output(gene_objects_dict, args.output_dir,
                                      args.output_file_base_name,
                                      show_sites = args.show_selected_sites,
                                           multimatrix = True)
    print('\n')
    if not args.no_pred_output: # skip this output if flag is true
        # make full file path of output predictions file
        preds_output_path = os.path.join(args.output_dir,
                                         args.output_file_base_name
                                         + '_species_predictions.csv')
        esl_int.generate_predictions_output(master_run_list,
                                            preds_output_path,
                                            args.species_pheno_path)

    ecf.report_elapsed_time(start_time) # print time taken for execution

    if args.make_sps_plot or args.make_sps_kde_plot:
        plot_type = 'violin' if args.make_sps_plot else 'kde'
        # generate and show density plots of predictions
        # call sps_density.create_sps_plots for various rmse cutoffs
        ecf.rmse_range_pred_plots(preds_output_path,
                                  args.output_file_base_name,
                                  args.pheno_names,
                                  args.min_genes,
                                  plot_type)

if __name__ == '__main__':
    main()