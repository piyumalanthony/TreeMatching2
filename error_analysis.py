import os
import pickle
import matplotlib.pyplot as plt
import scipy.stats as stat
import pandas as pd
import numpy as np

from constants import MIN_TAXA, MAX_TAXA, MIN_SEQ_LEN, MAX_SEQ_LEN
from tree_comp import generate_br_plots, generate_hessian_plots, generate_mle_plots, load_rmse_data, \
    generate_br_plots_compatible_tree, generate_hessian_plots_compatible_tree, generate_mle_plots_compatible_tree


def process_data_for_plots(file_path, min_taxa, max_taxa, min_seq, max_seq, gap, model, compatible_tree):
    mcmctree_output_file = 'out.BV'
    if model == 'WAG_Gamma':
        mcmctree_output_file = 'rst2'
    if not (os.path.exists(f'{file_path}/plots/br_plots')):
        os.mkdir(f'{file_path}/plots/br_plots')

    if not (os.path.exists(f'{file_path}/plots/hessian_plots')):
        os.mkdir(f'{file_path}/plots/hessian_plots')

    if not (os.path.exists(f'{file_path}/plots/likelihood_plots')):
        os.mkdir(f'{file_path}/plots/likelihood_plots')

    for num_taxa_iter in range(min_taxa, max_taxa + gap, gap):
        for seq_len_iter in range(min_seq, max_seq + min_seq, min_seq):
            hessian_list = []
            tree_string = ""
            branch_length_vector = ""
            gradient_vector = ""
            with open(f'{file_path}/mcmctree_output/{num_taxa_iter}/{seq_len_iter}/{mcmctree_output_file}') as f:
                for k, line in enumerate(f.readlines()):
                    if k == 3:
                        tree_string = line
                        # print(i, tree_string)
                    if k == 5:
                        branch_length_vector = line
                    if k == 7:
                        gradient_vector = line
                    if k > 11:
                        # print(i, line)
                        hessian_list.append(line)
                with open(f'{file_path}/mcmctree_output/{num_taxa_iter}/{seq_len_iter}/tree.nwk', "w") as f1:
                    f1.write(tree_string)
                with open(f'{file_path}/mcmctree_output/{num_taxa_iter}/{seq_len_iter}/brLengths.txt', "w") as f2:
                    f2.write(branch_length_vector)
                with open(f'{file_path}/mcmctree_output/{num_taxa_iter}/{seq_len_iter}/hessian.txt', "w") as f3:
                    f3.write("".join(hessian_list))
                with open(f'{file_path}/mcmctree_output/{num_taxa_iter}/{seq_len_iter}/gradient.txt', "w") as f4:
                    f4.write(gradient_vector)
            print(f'-------------- {num_taxa_iter} NUM TAXA, {seq_len_iter} SEQ LENGTH ----------------------------')
            if compatible_tree:
                generate_br_plots_compatible_tree(file_path, num_taxa_iter, seq_len_iter)
                print("                Completed Branch plot generation")
                generate_hessian_plots_compatible_tree(file_path, num_taxa_iter, seq_len_iter)
                print("                Completed Hessian plot generation")
                generate_mle_plots_compatible_tree(file_path, num_taxa_iter, seq_len_iter)
                print("                Completed MLE plot generation\n\n")
            else:
                generate_br_plots(file_path, num_taxa_iter, seq_len_iter)
                print("                Completed Branch plot generation")
                generate_hessian_plots(file_path, num_taxa_iter, seq_len_iter)
                print("                Completed Hessian plot generation")
                generate_mle_plots(file_path, num_taxa_iter, seq_len_iter)
                print("                Completed MLE plot generation\n\n")


def process_rmse_data(file_path, min_taxa, max_taxa, min_seq, max_seq, gap):
    if not (os.path.exists(f'{file_path}/err')):
        os.mkdir(f'{file_path}/err')

    rmse_br = []
    rmse_hessian_diag = []
    rmse_hessian_off_diag = []
    rmse_likelihood = []

    for num_taxa in range(min_taxa, max_taxa + gap, gap):
        rmse_br_local = []
        rmse_hessian_diag_local = []
        rmse_hessian_off_diag_local = []
        rmse_likelihood_local = []
        for seq_len in range(min_seq, max_seq + min_seq, min_seq):
            br_err, hessian_diag_err, hessian_off_diag_err, likelihood_err = load_rmse_data(file_path, num_taxa,
                                                                                            seq_len)
            rmse_br_local.append(br_err)
            rmse_hessian_diag_local.append(hessian_diag_err)
            rmse_hessian_off_diag_local.append(hessian_off_diag_err)
            rmse_likelihood_local.append(likelihood_err)
        else:
            rmse_br.append(rmse_br_local)
            rmse_hessian_diag.append(rmse_hessian_diag_local)
            rmse_hessian_off_diag.append(rmse_hessian_off_diag_local)
            rmse_likelihood.append(rmse_likelihood_local)

    with open(f'{file_path}/err/error_analysis.pkl', 'wb') as f:
        pickle.dump(
            {'rmse_br': rmse_br, 'rmse_hessian_diag': rmse_hessian_diag, 'rmse_hessian_off_diag': rmse_hessian_off_diag,
             'rmse_likelihood': rmse_likelihood}, f)


def rmse_box_plot_generation(file_path, model, min_num_taxa, max_num_taxa, gap):
    with open(f'{file_path}/err/error_analysis.pkl', 'rb') as f:
        pickle_data = pickle.load(f)

    rmse_br = pickle_data['rmse_br']
    rmse_hessian_diag = pickle_data['rmse_hessian_diag']
    rmse_hessian_off_diag = pickle_data['rmse_hessian_off_diag']
    rmse_likelihood = pickle_data['rmse_likelihood']
    lables = [i for i in range(min_num_taxa, max_num_taxa + gap, gap)]

    plt.boxplot(rmse_br, labels=lables)
    plt.title(f'RMSE Branch Lengths for {model} Substitution model')
    plt.xlabel('Number of Taxa')
    plt.ylabel('RMSE')
    plt.savefig(f'{file_path}/err/RMSE_br_{model}_box_plot.png')
    plt.clf()

    plt.boxplot(rmse_hessian_diag, labels=lables)
    plt.title(f'RMSE Hessian diagonal for {model} Substitution model')
    plt.xlabel('Number of Taxa')
    plt.ylabel('RMSE')
    plt.savefig(f'{file_path}/err/RMSE_hessian_diag_{model}_box_plot.png')
    plt.clf()

    plt.boxplot(rmse_hessian_off_diag, labels=lables)
    plt.title(f'RMSE Hessian off-diagonal for {model} Substitution model')
    plt.xlabel('Number of Taxa')
    plt.ylabel('RMSE')
    plt.savefig(f'{file_path}/err/RMSE_hessian_off_diag_{model}_box_plot.png')
    plt.clf()

    plt.boxplot(rmse_likelihood, labels=lables)
    plt.title(f'RMSE likelihood for {model} Substitution model')
    plt.xlabel('Number of Taxa')
    plt.ylabel('RMSE')
    plt.savefig(f'{file_path}/err/RMSE_likelihood_{model}_box_plot.png')
    plt.clf()


def calculate_correlation(file_path, min_taxa, max_taxa, min_seq, max_seq, gap):
    iqtree_global_br = []
    iqtree_global_h_diagonal = []
    iqtree_global_h_off_diagonal = []
    iqtree_global_likelihood = []
    mcmctree_global_br = []
    mcmctree_global_h_diagonal = []
    mcmctree_global_h_off_diagonal = []
    mcmctree_global_likelihood = []
    corr_df = pd.DataFrame(
        columns=['num_taxa', 'seq_len', 'br_corr', 'h_diagonal_corr', 'h_off_diagonal_corr', 'likelihood_corr'])
    index = 0
    for num_taxa in range(min_taxa, max_taxa + gap, gap):
        for seq_len in range(min_seq, max_seq + min_seq, min_seq):
            with open(f'{file_path}/pickle/br_idx_mapping_{num_taxa}_{seq_len}.pkl', 'rb') as f:
                pickle_data = pickle.load(f)
            iqtree_br = np.array(pickle_data["iqtree_blengths"])
            baseml_br = np.array(pickle_data["mcmctree_blengths"])

            with open(f'{file_path}/pickle/hessian_data_{num_taxa}_{seq_len}.pkl', 'rb') as f:
                pickle_data_h = pickle.load(f)
            iqtree_diagonal_h = np.array(pickle_data_h["iqtree_hessian_diagonal"])
            baseml_diagonal_h = np.array(pickle_data_h["baseml_hessian_diagonal"])
            iqtree_off_diagonal_h = np.array(pickle_data_h["iqtree_hessian_off_diagonal"])
            baseml_off_diagonal_h = np.array(pickle_data_h["baseml_hessian_off_diagonal"])

            with open(f'{file_path}/pickle/likelihood_data_{num_taxa}_{seq_len}.pkl', 'rb') as f:
                pickle_data_l = pickle.load(f)
            iqtree_likelihood = np.array(pickle_data_l["iqtree_mle"])
            baseml_likelihood = np.array(pickle_data_l["baseml_mle"])

            iqtree_global_br = iqtree_global_br + list(pickle_data["iqtree_blengths"])
            mcmctree_global_br = mcmctree_global_br + list(pickle_data["mcmctree_blengths"])

            iqtree_global_h_diagonal = iqtree_global_h_diagonal + list(pickle_data_h["iqtree_hessian_diagonal"])
            iqtree_global_h_off_diagonal = iqtree_global_h_off_diagonal + list(
                pickle_data_h["iqtree_hessian_off_diagonal"])
            mcmctree_global_h_diagonal = mcmctree_global_h_diagonal + list(pickle_data_h["baseml_hessian_diagonal"])
            mcmctree_global_h_off_diagonal = mcmctree_global_h_off_diagonal + list(
                pickle_data_h["baseml_hessian_off_diagonal"])

            iqtree_global_likelihood = iqtree_global_likelihood + list(pickle_data_l["iqtree_mle"])
            mcmctree_global_likelihood = mcmctree_global_likelihood + list(pickle_data_l["baseml_mle"])

            local_br_corr = stat.pearsonr(iqtree_br, baseml_br)[0]
            local_h_diagonal_corr = stat.pearsonr(iqtree_diagonal_h, baseml_diagonal_h)[0]
            local_h_off_diagonal_corr = stat.pearsonr(iqtree_off_diagonal_h, baseml_off_diagonal_h)[0]
            local_likelihood_corr = stat.pearsonr(iqtree_likelihood, baseml_likelihood)[0]

            corr_df.loc[index] = [num_taxa, seq_len, local_br_corr, local_h_diagonal_corr, local_h_off_diagonal_corr,
                                  local_likelihood_corr]
            index = index + 1

    else:
        global_br_corr = stat.pearsonr(iqtree_global_br, mcmctree_global_br)[0]
        global_h_diagonal_corr = stat.pearsonr(iqtree_global_h_diagonal, mcmctree_global_h_diagonal)[0]
        global_h_off_diagonal_corr = stat.pearsonr(iqtree_global_h_off_diagonal, mcmctree_global_h_off_diagonal)[0]
        global_likelihood_corr = stat.pearsonr(iqtree_global_likelihood, mcmctree_global_likelihood)[0]

        corr_df.loc[index] = ["ALL", "ALL", global_br_corr, global_h_diagonal_corr, global_h_off_diagonal_corr,
                              global_likelihood_corr]

        model = file_path.strip().split("/")[-1]
        corr_df.to_csv(f'{file_path}/err/{model}_corr.csv')
