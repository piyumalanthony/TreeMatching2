import os
import pickle
import time

import ete3
import psutil
import shutil

from constants import IQTREE_PATH


def process_memory():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss


def simulate_tree_topology(min_taxa, max_taxa, gap, file_str):
    for i in range(min_taxa, max_taxa + gap, gap):
        if not (os.path.exists(f'{file_str}/simulated_data/{i}')):
            os.mkdir(f'{file_str}/simulated_data/{i}')
        # os.system(f'mkdir {file_str}/simulated_data/{i}')
        taxa = ["T" + str(i) for i in range(1, i + 1)]
        for j in range(1000, 11000, 1000):
            t = ete3.Tree()
            t.populate(i, taxa, random_branches=True, branch_range=[0.0005, 0.4])
            t.write(format=1,
                    outfile=f'{file_str}/simulated_data/{i}/{j}.nw')


def run_alisim(path, min_num_taxa, max_num_taxa, min_seq_len, max_seq_len, gap, model):
    for i in range(min_num_taxa, max_num_taxa + gap, gap):
        for j in range(min_seq_len, max_seq_len + min_seq_len, min_seq_len):
            # if not (os.path.exists(f'{path}/simulated_data/{i}/{j}')):
            #     os.mkdir(f'{path}/simulated_data/{i}/{j}')
            cmd_alisim = f'{IQTREE_PATH} --alisim {path}/simulated_data/{i}/{j} -m {model} -t {path}/simulated_data/{i}/{j}.nw -seed 1 --length {j}'
            os.system(cmd_alisim)


def run_alisim_aa(path, min_num_taxa, max_num_taxa, min_seq_len, max_seq_len, gap, model):
    for i in range(min_num_taxa, max_num_taxa + gap, gap):
        for j in range(min_seq_len, max_seq_len + min_seq_len, min_seq_len):
            cmd_alisim = f'{IQTREE_PATH} --alisim {path}/simulated_data/{i}/{j} -m {model} --seqtype AA -t {path}/simulated_data/{i}/{j}.nw -seed 1 --length {j}'
            os.system(cmd_alisim)


def run_iqtree(file_path, model, min_taxa, max_taxa, min_seq_len, max_seq_len, gap):
    if not (os.path.exists(f'{file_path}/time_pickle')):
        os.mkdir(f'{file_path}/time_pickle')
    time_global = []
    memory_global = []
    for num_taxa in range(min_taxa, max_taxa + gap, gap):
        if not (os.path.exists(f'{file_path}/iqtree_output/{num_taxa}')):
            os.mkdir(f'{file_path}/iqtree_output/{num_taxa}')
        time_local = []
        memory_local = []
        for seq_len in range(min_seq_len, max_seq_len + min_seq_len, min_seq_len):
            if not (os.path.exists(f'{file_path}/iqtree_output/{num_taxa}/{seq_len}')):
                os.mkdir(f'{file_path}/iqtree_output/{num_taxa}/{seq_len}')
            print(f'*********** Running IQ-TREE for MODEL:{model}, NUM TAXA:{num_taxa}, SEQ LEN:{seq_len} ***********')
            cmd_iqtree = f'/usr/bin/time -v -o {file_path}/iqtree_output/{num_taxa}/{seq_len}/time.txt {IQTREE_PATH} -s {file_path}/simulated_data/{num_taxa}/{seq_len}.phy --redo  -nt 1 -m {model} --dating mcmctree -seed 1 -te {file_path}/simulated_data/{num_taxa}/{seq_len}.nw --prefix {file_path}/iqtree_output/{num_taxa}/{seq_len}/output -keep-ident'
            mem_before = process_memory()
            start = time.time()
            os.system(cmd_iqtree)
            end = time.time()
            mem_after = process_memory()
            time_local.append(end - start)
            memory_local.append(mem_after - mem_before)
        time_global.append(time_local)
        memory_global.append(memory_local)
        with open(f'{file_path}/time_pickle/iqtree_{model}_{num_taxa}.pkl', 'wb') as f:
            pickle.dump({'time': time_local, 'memory': memory_local}, f)
    with open(f'{file_path}/time_pickle/iqtree_{model}_all.pkl', 'wb') as f:
        pickle.dump({'time': time_global, 'memory': memory_global}, f)


def generate_ctl_mcmctree(file_path, model, min_taxa, max_taxa, min_seq_len, max_seq_len, gap, data_type):
    model_index, gamma_rate, seq_data = 0, 0, 0
    for num_taxa in range(min_taxa, max_taxa + gap, gap):
        if not (os.path.exists(f'{file_path}/mcmctree_output/{num_taxa}')):
            os.mkdir(f'{file_path}/mcmctree_output/{num_taxa}')
        for seq_len in range(min_seq_len, max_seq_len + min_seq_len, min_seq_len):
            if not (os.path.exists(f'{file_path}/mcmctree_output/{num_taxa}/{seq_len}')):
                os.mkdir(f'{file_path}/mcmctree_output/{num_taxa}/{seq_len}')
            with open(f'{file_path}/simulated_data/{num_taxa}/{seq_len}.nw') as f:
                tree_str = f.readline()
                with open(f'{file_path}/simulated_data/{num_taxa}/{seq_len}_mcmc.treefile', "w+") as f:
                    f.write("{} {}\n".format(num_taxa, 1))
                    f.write(tree_str)
            if model == 'JC':
                model_index = 0
                gamma_rate = 0
            elif model == 'HKY':
                model_index = 4
                gamma_rate = 0
            if model == 'JC+G5{0.5}':
                model_index = 0
                gamma_rate = 0.5
            elif model == 'HKY+G5{0.5}':
                model_index = 4
                gamma_rate = 0.5
            elif model == 'WAG+G5{0.5}':
                model_index = 0
                gamma_rate = 0.5

            if data_type == 'DNA':
                seq_data = 0
            elif data_type == 'Codon':
                seq_data = 1
            elif data_type == 'AA':
                seq_data = 2

            mcmctree_file_str = [
                "seed = -1\n"
                f"seqfile = {file_path}/simulated_data/{num_taxa}/{seq_len}.phy\n",
                f"treefile = {file_path}/simulated_data/{num_taxa}/{seq_len}_mcmc.treefile\n",
                "outfile = out\n",
                "ndata = 1\n",
                f"seqtype = {seq_data}\n",
                "usedata = 3\n",
                "clock = 3\n",
                "RootAge = <1.0\n",
                f"model = {model_index}\n",
                f"alpha = {gamma_rate}\n",
                "ncatG = 5\n",
                "cleandata = 0\n",
                "BDparas = 1 1 0\n",
                "kappa_gamma = 6 2\n",
                "alpha_gamma = 1 1\n",
                "rgene_gamma = 2 2\n",
                "sigma2_gamma = 1 10\n",
                "finetune = 1: 0.1  0.1  0.1  0.01 .5\n",
                "print = 1\n",
                "burnin = 2000\n",
                "sampfreq = 2\n",
                "nsample = 20000\n"

            ]
            with open(f"{file_path}/simulated_data/{num_taxa}/{seq_len}_mcmctree.ctl", 'w+') as f:
                f.write("".join(mcmctree_file_str))


def generate_ctl_codeml(file_path, num_taxa, seq_len):
    os.remove(f'{file_path}/mcmctree_output/{num_taxa}/{seq_len}/tmp0001.ctl')
    os.remove(f'{file_path}/mcmctree_output/{num_taxa}/{seq_len}/out.BV')
    codeml_file_str = [
        'seqfile = tmp0001.txt\n',
        'treefile = tmp0001.trees\n',
        'outfile = tmp0001.out\n',
        'noisy = 3\n',
        'seqtype = 2\n',
        'model = 2\n'
        f'aaRatefile = {file_path}/wag.dat\n'
        'fix_alpha = 0\n',
        'alpha = .5\n',
        'ncatG = 5\n',
        'Small_Diff = 0.1e-6\n',
        'getSE = 2\n',
        'method = 1'
    ]

    with open(f"{file_path}/mcmctree_output/{num_taxa}/{seq_len}/tmp0001.ctl", 'w+') as f:
        f.write("".join(codeml_file_str))


def run_mcmctree(file_path, model, min_taxa, max_taxa, min_seq_len, max_seq_len, gap, data_type):
    time_global = []
    memory_global = []
    if model == 'WAG+G5{0.5}':
        shutil.copy('./dat/wag.dat', file_path)
    generate_ctl_mcmctree(file_path, model, min_taxa, max_taxa, min_seq_len, max_seq_len, gap, data_type)
    for num_taxa in range(min_taxa, max_taxa + gap, gap):
        time_local = []
        memory_local = []
        for seq_len in range(min_seq_len, max_seq_len + min_seq_len, min_seq_len):
            print(f'*********** Running BaseML for MODEL:{model}, NUM TAXA:{num_taxa}, SEQ LEN:{seq_len} ***********')
            os.chdir(f'{file_path}/mcmctree_output/{num_taxa}/{seq_len}')
            cmd_mcmctree = f'/usr/bin/time -v -o {file_path}/mcmctree_output/{num_taxa}/{seq_len}/time.txt mcmctree {file_path}/simulated_data/{num_taxa}/{seq_len}_mcmctree.ctl'
            mem_before = process_memory()
            start = time.time()
            os.system(cmd_mcmctree)
            if model == 'WAG+G5{0.5}':
                generate_ctl_codeml(file_path, num_taxa, seq_len)
                cmd_codeml = f'/usr/bin/time -v -o {file_path}/mcmctree_output/{num_taxa}/{seq_len}/time.txt codeml {file_path}/mcmctree_output/{num_taxa}/{seq_len}/tmp0001.ctl'
                os.system(cmd_codeml)
            end = time.time()
            mem_after = process_memory()
            time_local.append(end - start)
            memory_local.append(mem_after - mem_before)
        time_global.append(time_local)
        memory_global.append(memory_local)
        with open(f'{file_path}/time_pickle/mcmctree_{model}_{num_taxa}.pkl', 'wb') as f:
            pickle.dump({'time': time_local, 'memory': memory_local}, f)
    with open(f'{file_path}/time_pickle/mcmctree_{model}_all.pkl', 'wb') as f:
        pickle.dump({'time': time_global, 'memory': memory_global}, f)
