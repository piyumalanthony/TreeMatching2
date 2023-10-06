import argparse
import os

from constants import MIN_TAXA, MAX_TAXA, MIN_SEQ_LEN, MAX_SEQ_LEN, TAXA_GAP, NUM_REPLICATES, IQTREE_PREV_RELEASE_PATH, \
    IQTREE_PATH
from data_generation import simulate_tree_topology, run_alisim, run_iqtree, run_mcmctree, run_alisim_aa
from error_analysis import process_data_for_plots, process_rmse_data, rmse_box_plot_generation, calculate_correlation
from iqtree_release_test import generate_tree_replicates_seq_len_fixed, generate_alignments_fixed_seq_len, \
    run_iqtree_fixed_seq_len, run_iqtree_emperical

parser = argparse.ArgumentParser()
parser.add_argument('--data_type', help='Sequence data type', choices=['DNA', 'AA'])
parser.add_argument('--model', help='Substitution model for simulations',
                    choices=['JC', 'HKY', 'WAG', 'JC_Gamma', 'HKY_Gamma', 'WAG_Gamma', 'Poisson'])
# parser.add_argument('--multiple_models', help='use all substitution model for simulations',
#                     choices=[True, False])
parser.add_argument('--file_path', help='file path to store experimental data')
parser.add_argument('--iqtree_path', help='IQ-TREE path for experiments')
parser.add_argument('--generate_data', help='Flag to generate data for the simulations')
parser.add_argument('--min_taxa', help='Minimum number of taxa')
parser.add_argument('--max_taxa', help='Maximum number of taxa')
parser.add_argument('--min_seq_len', help='Minimum sequence length')
parser.add_argument('--max_seq_len', help='Maximum sequence length')
parser.add_argument('--gap', help='Gap of number of taxa')
parser.add_argument('--err', help='Error analysis flag for comparison')
parser.add_argument('--model_list', help='List of models for error analysis')

# Params for performance testing
parser.add_argument('--iqtree_test', help='flag for performance testing compared to previous version')
parser.add_argument('--num_rep', help='Number of replicates for performance test')
parser.add_argument('--seq_len', help='Sequence length for fixed seq_len test')

# Empirical data test for IQ-TREE
parser.add_argument('--empirical_test', help='flag for empirical testing compared to previous version')
parser.add_argument('--num_emp_rep', help='Number of replicates for empirical test')
parser.add_argument('--dataset', help='Name of the dataset')
parser.add_argument('--dataset_path', help='Path to the dataset')

args = parser.parse_args()

min_taxa = MIN_TAXA
max_taxa = MAX_TAXA
min_seq_len = MIN_SEQ_LEN
max_seq_len = MAX_SEQ_LEN
gap = TAXA_GAP
num_replicates = NUM_REPLICATES

if args.min_taxa is not None:
    min_taxa = int(args.min_taxa)
if args.max_taxa is not None:
    max_taxa = int(args.max_taxa)
if args.min_seq_len is not None:
    min_seq_len = int(args.min_seq_len)
if args.max_seq_len is not None:
    max_seq_len = int(args.max_seq_len)
if args.gap is not None:
    gap = int(args.gap)
if args.num_rep is not None:
    num_replicates = int(args.num_rep)

# generate data for DNA data
if args.generate_data and args.data_type == 'DNA':
    dir_model = f'{args.file_path}/{args.model}'
    dir_iqtree_output = f'{args.file_path}/{args.model}/iqtree_output'
    dir_mcmctree_output = f'{args.file_path}/{args.model}/mcmctree_output'
    dir_pickle = f'{args.file_path}/{args.model}/pickle'
    dir_plots = f'{args.file_path}/{args.model}/plots'
    dir_simulated_data = f'{args.file_path}/{args.model}/simulated_data'

    if not (os.path.exists(dir_model)):
        os.mkdir(dir_model)

    if not (os.path.exists(dir_iqtree_output)):
        os.mkdir(dir_iqtree_output)

    if not (os.path.exists(dir_mcmctree_output)):
        os.mkdir(dir_mcmctree_output)

    if not (os.path.exists(dir_pickle)):
        os.mkdir(dir_pickle)

    if not (os.path.exists(dir_plots)):
        os.mkdir(dir_plots)

    if not (os.path.exists(dir_simulated_data)):
        os.mkdir(dir_simulated_data)

    simulate_tree_topology(min_taxa=min_taxa, max_taxa=max_taxa, gap=gap,
                           file_str=dir_model)

    model_str = args.model
    if model_str == 'JC_Gamma':
        s_model = 'JC+G5{0.5}'
    elif model_str == 'HKY_Gamma':
        s_model = 'HKY+G5{0.5}'
    else:
        s_model = model_str
    run_alisim(path=dir_model, min_num_taxa=min_taxa, max_num_taxa=max_taxa, min_seq_len=min_seq_len,
               max_seq_len=max_seq_len, gap=gap, model=s_model)
    run_iqtree(file_path=dir_model, model=s_model, min_taxa=min_taxa, max_taxa=max_taxa,
               min_seq_len=min_seq_len, max_seq_len=max_seq_len, gap=gap)

    run_mcmctree(file_path=dir_model, model=s_model, min_taxa=min_taxa, max_taxa=max_taxa, min_seq_len=min_seq_len,
                 max_seq_len=max_seq_len, gap=gap, data_type=args.data_type)

if args.generate_data and args.data_type == 'AA':
    dir = f'{args.file_path}/AA'
    dir_model = f'{args.file_path}/AA/{args.model}'
    dir_iqtree_output = f'{args.file_path}/AA/{args.model}/iqtree_output'
    dir_mcmctree_output = f'{args.file_path}/AA/{args.model}/mcmctree_output'
    dir_pickle = f'{args.file_path}/AA/{args.model}/pickle'
    dir_plots = f'{args.file_path}/AA/{args.model}/plots'
    dir_simulated_data = f'{args.file_path}/AA/{args.model}/simulated_data'

    if not (os.path.exists(dir)):
        os.mkdir(dir)

    if not (os.path.exists(dir_model)):
        os.mkdir(dir_model)

    if not (os.path.exists(dir_iqtree_output)):
        os.mkdir(dir_iqtree_output)

    if not (os.path.exists(dir_mcmctree_output)):
        os.mkdir(dir_mcmctree_output)

    if not (os.path.exists(dir_pickle)):
        os.mkdir(dir_pickle)

    if not (os.path.exists(dir_plots)):
        os.mkdir(dir_plots)

    if not (os.path.exists(dir_simulated_data)):
        os.mkdir(dir_simulated_data)

    simulate_tree_topology(min_taxa=min_taxa, max_taxa=max_taxa, gap=gap,
                           file_str=dir_model)

    model_str = args.model
    if model_str == 'JC_Gamma':
        s_model = 'JC+G5{0.5}'
    elif model_str == 'HKY_Gamma':
        s_model = 'HKY+G5{0.5}'
    elif model_str == 'WAG_Gamma':
        s_model = 'WAG+G5{0.5}'
    else:
        s_model = model_str
    run_alisim_aa(path=dir_model, min_num_taxa=min_taxa, max_num_taxa=max_taxa, min_seq_len=min_seq_len,
               max_seq_len=max_seq_len, gap=gap, model=s_model)
    run_iqtree(file_path=dir_model, model=s_model, min_taxa=min_taxa, max_taxa=max_taxa,
               min_seq_len=min_seq_len, max_seq_len=max_seq_len, gap=gap)

    run_mcmctree(file_path=dir_model, model=s_model, min_taxa=min_taxa, max_taxa=max_taxa, min_seq_len=min_seq_len,
                 max_seq_len=max_seq_len, gap=gap, data_type=args.data_type)

if args.err:
    model_str = args.model_list
    model_list = model_str.strip().split(",")
    print(f'*********** Set of models to be analyzed:{model_list} ***********')

    file_str_list = [f'{args.file_path}/{model}' for model in model_list]
    for model, file_path in zip(model_list, file_str_list):
        print(f'*********** Processing data for model:{model} *********** \n\n')
        process_data_for_plots(file_path, min_taxa, max_taxa, min_seq_len, max_seq_len, gap)
        process_rmse_data(file_path, min_taxa, max_taxa, min_seq_len, max_seq_len, gap)
        rmse_box_plot_generation(file_path, model, min_taxa, max_taxa, gap)
        calculate_correlation(file_path, min_taxa, max_taxa, min_seq_len, max_seq_len, gap)

if args.iqtree_test:
    dir_experiment = f'{args.file_path}'
    seq_len = int(args.seq_len)
    generate_tree_replicates_seq_len_fixed(num_replicates, min_taxa, max_taxa, gap, dir_experiment)
    generate_alignments_fixed_seq_len(dir_experiment, num_replicates, seq_len, min_taxa, max_taxa, gap)
    print(f'|| *********** Running IQ-TREE Release version  *********** ||')
    run_iqtree_fixed_seq_len(dir_experiment, IQTREE_PREV_RELEASE_PATH, min_taxa, max_taxa, num_replicates, gap, "v1")
    print(f'|| *********** Running IQ-TREE New version  *********** ||')
    run_iqtree_fixed_seq_len(dir_experiment, IQTREE_PATH, min_taxa, max_taxa, num_replicates, gap, "v2")

if args.empirical_test:
    dir_experiment = f'{args.file_path}'
    # dataset_name = str(args.dataset)
    data_path = f'{args.dataset_path}'
    dataset1_path = f'{data_path}/74.aa'
    print(f'|| *********** Running IQ-TREE Release version  *********** ||')
    run_iqtree_emperical(dir_experiment, IQTREE_PREV_RELEASE_PATH, 'aa_74_MF', dataset1_path, 'v1')
    print(f'|| *********** Running IQ-TREE New version  *********** ||')
    run_iqtree_emperical(dir_experiment, IQTREE_PATH, 'aa_74_MF', dataset1_path, 'v2')

    dataset2_path = f'{data_path}/218.dna'
    print(f'|| *********** Running IQ-TREE Release version  *********** ||')
    run_iqtree_emperical(dir_experiment, IQTREE_PREV_RELEASE_PATH, 'dna_218_MF', dataset2_path, 'v1')
    print(f'|| *********** Running IQ-TREE New version  *********** ||')
    run_iqtree_emperical(dir_experiment, IQTREE_PATH, 'dna_218_MF', dataset2_path, 'v2')
    print(f'|| *********** Completed the empirical experiment  *********** ||')



