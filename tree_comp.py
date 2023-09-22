import os

import matplotlib.pyplot as plt
from ete3 import Tree
import numpy as np
import random
import copy
import pickle

from likelihood import calculate_approx_likelihood

random.seed(1)


def rmse(arr1, arr2):
    arr_1 = np.array(arr1)
    arr_2 = np.array(arr2)
    err = arr_1 - arr_2
    err_sqr = np.power(err, 2)
    err_sqr_mean = np.sum(err_sqr) / np.shape(err_sqr)[0]
    return np.sqrt(err_sqr_mean)


def generate_br_plots(file_path, num_taxa, seq_len):
    iqtree_br = f'{file_path}/iqtree_output/{num_taxa}/{seq_len}/output_blengths.gh'
    mcmctree_br = f'{file_path}/mcmctree_output/{num_taxa}/{seq_len}/brLengths.txt'
    with open(iqtree_br) as f:
        data = f.readline()
        iqtree_blengths = data.strip().split()
        iqtree_blengths = [float(i) for i in iqtree_blengths]

    with open(mcmctree_br) as f:
        data = f.readline()
        mcmctree_blengths = data.strip().split()
        mcmctree_blengths = [float(i) for i in mcmctree_blengths]

    iqtree_file = f'{file_path}/iqtree_output/{num_taxa}/{seq_len}/output.treefile'
    mcmctree_file = f'{file_path}/mcmctree_output/{num_taxa}/{seq_len}/tree.nwk'

    tree_iqtree = Tree(iqtree_file)

    iqtree_node_traversal = []
    for node in tree_iqtree.traverse("preorder"):
        edge = node.get_edges()[0]
        clade_1 = list(edge[0])
        clade_2 = list(edge[1])
        node_name = clade_1 + clade_2
        node_name = [node.name for node in node_name]
        node_name.sort()
        node_name = ':'.join(node_name)
        node.name = node_name
        iqtree_node_traversal.append(node.name)

    tree_mcmctree = Tree(mcmctree_file)

    mcmctree_node_traversal = []
    for node in tree_mcmctree.traverse("preorder"):
        edge = node.get_edges()[0]
        clade_1 = list(edge[0])
        clade_2 = list(edge[1])
        node_name = clade_1 + clade_2
        node_name = [node.name for node in node_name]
        node_name.sort()
        node_name = ':'.join(node_name)
        node.name = node_name
        mcmctree_node_traversal.append(node.name)

    iqtree_diff = []
    mcmctree_diff = []
    for u in mcmctree_node_traversal:
        if u not in iqtree_node_traversal:
            iqtree_diff.append(u)

    for v in iqtree_node_traversal:
        if v not in mcmctree_node_traversal:
            mcmctree_diff.append(v)
    iqtree_diff.reverse()
    mcmctree_node_traversal_v1 = copy.deepcopy(mcmctree_node_traversal)
    for idx_iter, item in enumerate(iqtree_diff):
        index_diff = mcmctree_node_traversal.index(item)
        mcmctree_node_traversal[index_diff] = mcmctree_diff[idx_iter]

    iqtree_idx = []
    mcmctree_node_traversal_post = mcmctree_node_traversal[1:]

    for val in iqtree_node_traversal[1:]:
        iqtree_idx.append(mcmctree_node_traversal_post.index(val))
    mcmc_blengths_v2 = []
    for idx in iqtree_idx:
        mcmc_blengths_v2.append(mcmctree_blengths[idx])

    iqtree_blengths_np = np.array(iqtree_blengths)
    mcmctree_blengths_np = np.array(mcmc_blengths_v2)

    blengths_diff = iqtree_blengths_np - mcmctree_blengths_np
    # pickle_data = [iqtree_idx, iqtree_node_traversal, mcmctree_node_traversal_v1, iqtree_diff, mcmctree_diff,
    #                blengths_diff, iqtree_blengths_np, mcmctree_blengths_np, iqtree_blengths, mcmctree_blengths]
    pickle_data = {"iqtree_idx": iqtree_idx,
                   "iqtree_node_traversal": iqtree_node_traversal,
                   "mcmctree_node_traversal": mcmctree_node_traversal_v1,  # initial mode traversal for MCMC-Tree
                   "iqtree_diff": iqtree_diff,
                   "mcmctree_diff": mcmctree_diff,
                   "blengths_diff": blengths_diff,
                   "iqtree_blengths": iqtree_blengths_np,
                   "mcmctree_blengths": mcmctree_blengths_np,  # reordered branch lengths
                   "iqtree_blengths_init": iqtree_blengths,
                   "mcmctree_blengths_init": mcmctree_blengths,

                   }

    with open(f'{file_path}/pickle/br_idx_mapping_{num_taxa}_{seq_len}.pkl', 'wb') as f:
        pickle.dump(pickle_data, f)

    plt.plot(iqtree_blengths_np, mcmctree_blengths_np, 'bo', )
    plt.plot([min(iqtree_blengths_np), max(iqtree_blengths_np)], [min(iqtree_blengths_np), max(iqtree_blengths_np)],
             'r--', label='x=y')
    plt.title(f'IQTree vs Baseml branch lengths for {num_taxa} taxa and {seq_len} seq length')
    plt.xlabel('IQtree branch lengths')
    plt.ylabel('Baseml branch lengths')
    plt.savefig(f'{file_path}/plots/br_plots/br_lengths_{num_taxa}_{seq_len}.png')
    plt.clf()


def generate_hessian_plots(file_path, num_taxa, seq_len):
    with open(f'{file_path}/pickle/br_idx_mapping_{num_taxa}_{seq_len}.pkl', 'rb') as f:
        pickle_data = pickle.load(f)
    index_list = pickle_data["iqtree_idx"]

    iqtree_gradients = f'{file_path}/iqtree_output/{num_taxa}/{seq_len}/output_gradient.gh'
    baseml_gradients = f'{file_path}/mcmctree_output/{num_taxa}/{seq_len}/gradient.txt'

    with open(iqtree_gradients) as f:
        data = f.readline()
        iqtree_gradient_vals = data.strip().split()
        iqtree_gradient_vals = [float(i) for i in iqtree_gradient_vals]

    with open(baseml_gradients) as f:
        data = f.readline()
        baseml_gradient_vals = data.strip().split()
        baseml_gradient_vals = [float(i) for i in baseml_gradient_vals]

    gradients = {"iqtree_gradients": iqtree_gradient_vals, "baseml_gradients": baseml_gradient_vals}

    with open(f'{file_path}/pickle/gradient_data_{num_taxa}_{seq_len}.pkl', 'wb') as f:
        pickle.dump(gradients, f)

    iqtree_hessian = f'{file_path}/iqtree_output/{num_taxa}/{seq_len}/output_hessian.gh'
    baseml_hessian = f'{file_path}/mcmctree_output/{num_taxa}/{seq_len}/hessian.txt'

    iqtree_hessian_values = []
    iqtree_h = open(iqtree_hessian, 'r')
    lines_iqtree = iqtree_h.readlines()
    iqtree_h.close()

    for line in lines_iqtree:
        values = line.strip().split()
        values = [float(i) for i in values]
        iqtree_hessian_values.append(values)

    baseml_hessian_values = []
    baseml_h = open(baseml_hessian, 'r')
    lines_baseml = baseml_h.readlines()
    baseml_h.close()

    for line in lines_baseml:
        values = line.strip().split()
        values = [float(i) for i in values]
        baseml_hessian_values.append(values)

    baseml_hessian_global = []
    for m in range(len(index_list)):
        baseml_hessian_local = []
        for n in range(len(index_list)):
            baseml_hessian_local.append(baseml_hessian_values[index_list[m]][index_list[n]])
        else:
            baseml_hessian_global.append(baseml_hessian_local)

    iqtree_hessian_diagonal = [element[idx] for idx, element in enumerate(iqtree_hessian_values)]
    baseml_hessian_diagonal = [element[idx] for idx, element in enumerate(baseml_hessian_global)]

    iqtree_hessian_off_diagonal = []
    for u, v in enumerate(iqtree_hessian_values):
        for k, l in enumerate(v):
            if k != u:
                iqtree_hessian_off_diagonal.append(l)

    baseml_hessian_off_diagonal = []
    for u, v in enumerate(baseml_hessian_global):
        for k, l in enumerate(v):
            if k != u:
                baseml_hessian_off_diagonal.append(l)

    hessians = {"iqtree_hessian": iqtree_hessian_values,
                "baseml_hessian": baseml_hessian_values,
                "iqtree_hessian_diagonal": iqtree_hessian_diagonal,
                "baseml_hessian_diagonal": baseml_hessian_diagonal,
                "iqtree_hessian_off_diagonal": iqtree_hessian_off_diagonal,
                "baseml_hessian_off_diagonal": baseml_hessian_off_diagonal}

    with open(f'{file_path}/pickle/hessian_data_{num_taxa}_{seq_len}.pkl', 'wb') as f:
        pickle.dump(hessians, f)

    plt.plot(iqtree_hessian_diagonal, baseml_hessian_diagonal, 'bo')
    plt.plot([min(iqtree_hessian_diagonal), max(iqtree_hessian_diagonal)],
             [min(iqtree_hessian_diagonal), max(iqtree_hessian_diagonal)], 'r--', label='x=y')
    plt.title(f'IQTree vs Baseml hessian diagonal values for {num_taxa} taxa and {seq_len} seq length')
    plt.xlabel('IQtree hessian')
    plt.ylabel('Baseml hessian')
    plt.savefig(f'{file_path}/plots/hessian_plots/hessian_diagonal_{num_taxa}_{seq_len}.png')
    plt.clf()

    plt.plot(iqtree_hessian_off_diagonal, baseml_hessian_off_diagonal, 'bo')
    plt.plot([min(iqtree_hessian_off_diagonal), max(iqtree_hessian_off_diagonal)],
             [min(iqtree_hessian_off_diagonal), max(iqtree_hessian_off_diagonal)], 'r--', label='x=y')
    plt.title(f'IQTree vs Baseml hessian off diagonal values for {num_taxa} taxa and {seq_len} seq length')
    plt.xlabel('IQtree hessian')
    plt.ylabel('Baseml hessian')
    plt.savefig(f'{file_path}/plots/hessian_plots/hessian_off_diagonal_{num_taxa}_{seq_len}.png')


def generate_mle_plots(file_path, num_taxa, seq_len):
    with open(f'{file_path}/pickle/br_idx_mapping_{num_taxa}_{seq_len}.pkl', 'rb') as f:
        pickle_data_br = pickle.load(f)
    iqtree_br = pickle_data_br["iqtree_blengths_init"]
    baseml_br = pickle_data_br["mcmctree_blengths_init"]
    iqtree_idx = pickle_data_br["iqtree_idx"]

    with open(f'{file_path}/pickle/gradient_data_{num_taxa}_{seq_len}.pkl', 'rb') as f:
        pickle_data_gr = pickle.load(f)
    iqtree_gradients = pickle_data_gr["iqtree_gradients"]
    baseml_gradients = pickle_data_gr["baseml_gradients"]

    with open(f'{file_path}/pickle/hessian_data_{num_taxa}_{seq_len}.pkl', 'rb') as f:
        pickle_data_v3 = pickle.load(f)
    iqtree_hessian = pickle_data_v3["iqtree_hessian"]
    baseml_hessian = pickle_data_v3["baseml_hessian"]

    iq_tree_mle, baseml_mle, iq_tree_mle_delta, baseml_mle_delta = calculate_approx_likelihood([iqtree_br, baseml_br],
                                                                                               [iqtree_gradients,
                                                                                                baseml_gradients],
                                                                                               [iqtree_hessian,
                                                                                                baseml_hessian],
                                                                                               iqtree_idx)
    mle_data = {
        "iqtree_mle": iq_tree_mle,
        "baseml_mle": baseml_mle
    }
    with open(f'{file_path}/pickle/likelihood_data_{num_taxa}_{seq_len}.pkl', 'wb') as f:
        pickle.dump(mle_data, f)

    plt.plot(iq_tree_mle, baseml_mle, 'bo')
    plt.plot([min(iq_tree_mle), max(iq_tree_mle)],
             [min(iq_tree_mle), max(iq_tree_mle)], 'r--', label='x=y')
    plt.title(f'IQTree vs Baseml likelihood values for {num_taxa} taxa and {seq_len} seq length')
    plt.xlabel('IQtree likelihood')
    plt.ylabel('Baseml likelihood')
    plt.savefig(f'{file_path}/plots/likelihood_plots/likelihood_{num_taxa}_{seq_len}.png')
    plt.clf()

    plt.plot(iq_tree_mle_delta, baseml_mle_delta, 'bo')
    plt.plot([min(iq_tree_mle_delta), max(iq_tree_mle_delta)],
             [min(iq_tree_mle_delta), max(iq_tree_mle_delta)], 'r--', label='x=y')
    plt.title(
        f'IQTree vs Baseml likelihood value for random branch lengths for {num_taxa} taxa and {seq_len} seq length')
    plt.xlabel('IQtree likelihood')
    plt.ylabel('Baseml likelihood')
    plt.savefig(f'{file_path}/plots/likelihood_plots/likelihood_random_{num_taxa}_{seq_len}.png')
    plt.clf()


def load_rmse_data(file_path, num_taxa, seq_len):
    print("----------------------------------- RMSE calculation --------------------------------------------")
    with open(f'{file_path}/pickle/br_idx_mapping_{num_taxa}_{seq_len}.pkl', 'rb') as f:
        pickle_data = pickle.load(f)
    iqtree_br = np.array(pickle_data["iqtree_blengths"])
    baseml_br = np.array(pickle_data["mcmctree_blengths"])
    error_br = rmse(iqtree_br, baseml_br)

    with open(f'{file_path}/pickle/hessian_data_{num_taxa}_{seq_len}.pkl', 'rb') as f:
        pickle_data_h = pickle.load(f)

    iqtree_diagonal_h = np.array(pickle_data_h["iqtree_hessian_diagonal"])
    baseml_diagonal_h = np.array(pickle_data_h["baseml_hessian_diagonal"])
    iqtree_off_diagonal_h = np.array(pickle_data_h["iqtree_hessian_off_diagonal"])
    baseml_off_diagonal_h = np.array(pickle_data_h["baseml_hessian_off_diagonal"])

    error_diagonal_h = rmse(iqtree_diagonal_h, baseml_diagonal_h)
    error_off_diagonal_h = rmse(iqtree_off_diagonal_h, baseml_off_diagonal_h)

    with open(f'{file_path}/pickle/likelihood_data_{num_taxa}_{seq_len}.pkl', 'rb') as f:
        pickle_data_l = pickle.load(f)
    iqtree_likelihood = pickle_data_l["iqtree_mle"]
    baseml_likelihood = pickle_data_l["baseml_mle"]

    error_likelihood = rmse(iqtree_likelihood, baseml_likelihood)

    return error_br, error_diagonal_h, error_off_diagonal_h, error_likelihood
