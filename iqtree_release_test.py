import os
import pickle
import time

import ete3

from constants import IQTREE_PATH
from data_generation import process_memory


def generate_tree_replicates_seq_len_fixed(num_replicates, min_taxa, max_taxa, gap, file_str):
    print('|| *********** Initialized tree topology generation for fix seq size *********** || \n\n')
    if not (os.path.exists(f'{file_str}')):
        os.mkdir(f'{file_str}')
    if not (os.path.exists(f'{file_str}/simulated_data')):
        os.mkdir(f'{file_str}/simulated_data')
    for i in range(min_taxa, max_taxa + gap, gap):
        if not (os.path.exists(f'{file_str}/simulated_data/{i}')):
            os.mkdir(f'{file_str}/simulated_data/{i}')
        for j in range(num_replicates):
            taxa = ["T" + str(i) for i in range(1, i + 1)]
            t = ete3.Tree()
            t.populate(size=i, names_library=taxa, random_branches=True, branch_range=[0.0005, 0.4])
            t.write(format=1,
                    outfile=f'{file_str}/simulated_data/{i}/{i}_{j}.nw')
        print(f'*********** Generated {num_replicates} Replicates for NUM TAXA: {i} *********** \n\n')


def generate_tree_replicates_num_taxa_fixed(num_replicates, num_taxa, min_seq_len, max_seq_len, gap_seq_len, file_str):
    # print('*********** Initialized tree topology generation for fix seq size *********** \n\n')
    if not (os.path.exists(f'{file_str}')):
        os.mkdir(f'{file_str}')
    if not (os.path.exists(f'{file_str}/simulated_data')):
        os.mkdir(f'{file_str}/simulated_data')
    for i in range(min_seq_len, max_seq_len + gap_seq_len, gap_seq_len):
        if not (os.path.exists(f'{file_str}/simulated_data/{i}')):
            os.mkdir(f'{file_str}/simulated_data/{i}')
        for j in range(num_replicates):
            taxa = ["T" + str(i) for i in range(1, num_taxa + 1)]
            t = ete3.Tree()
            t.populate(size=i, names_library=taxa, random_branches=True, branch_range=[0.0005, 0.4])
            t.write(format=1,
                    outfile=f'{file_str}/simulated_data/{i}/{i}_{j}.nw')


def generate_alignments_fixed_seq_len(path, n_replicates, seq_len, min_taxa, max_taxa, gap):
    print('|| *********** Initialized Alignment generation with alisim for fix seq size *********** || \n\n')
    for i in range(min_taxa, max_taxa + gap, gap):
        for j in range(n_replicates):
            cmd_alisim = f'{IQTREE_PATH} --alisim {path}/simulated_data/{i}/{i}_{j} -m GTR -t {path}/simulated_data/{i}/{i}_{j}.nw -seed 1 --length {seq_len}'
            os.system(cmd_alisim)
        print(f'*********** Generated {n_replicates} Replicates of Alignments for NUM TAXA: {i} *********** \n\n')


def run_iqtree_fixed_seq_len(file_path, iqtree_path, min_taxa, max_taxa, n_replicates, gap, version):
    if not (os.path.exists(f'{file_path}/time_pickle_{version}')):
        os.mkdir(f'{file_path}/time_pickle_{version}')
    if not (os.path.exists(f'{file_path}/iqtree_output_{version}')):
        os.mkdir(f'{file_path}/iqtree_output_{version}')

    time_global = []
    memory_global = []

    for num_taxa in range(min_taxa, max_taxa + gap, gap):
        if not (os.path.exists(f'{file_path}/iqtree_output_{version}/{num_taxa}')):
            os.mkdir(f'{file_path}/iqtree_output_{version}/{num_taxa}')
        time_local = []
        memory_local = []
        for replicate in range(n_replicates):
            if not (os.path.exists(f'{file_path}/iqtree_output_{version}/{num_taxa}/{num_taxa}_{replicate}')):
                os.mkdir(f'{file_path}/iqtree_output_{version}/{num_taxa}/{num_taxa}_{replicate}')
            print(f'*********** Running IQ-TREE for NUM TAXA:{num_taxa}, REPLICATE:{replicate} ***********')
            cmd_iqtree = f'{iqtree_path} -s {file_path}/simulated_data/{num_taxa}/{num_taxa}_{replicate}.phy --redo  -nt 1 -seed 1 -te {file_path}/simulated_data/{num_taxa}/{num_taxa}_{replicate}.nw --prefix {file_path}/iqtree_output_{version}/{num_taxa}/{num_taxa}_{replicate}/output'
            mem_before = process_memory()
            start = time.time()
            os.system(cmd_iqtree)
            end = time.time()
            mem_after = process_memory()
            time_local.append(end - start)
            memory_local.append(mem_after - mem_before)
        time_global.append(time_local)
        memory_global.append(memory_local)
        with open(f'{file_path}/time_pickle_{version}/iqtree_{version}_{num_taxa}.pkl', 'wb') as f:
            pickle.dump({'time': time_local, 'memory': memory_local}, f)
    with open(f'{file_path}/time_pickle_{version}/iqtree_{version}_all.pkl', 'wb') as f:
        pickle.dump({'time': time_global, 'memory': memory_global}, f)


def run_iqtree_emperical(file_path, iqtree_path, dataset_name, data_path, version):
    if not (os.path.exists(f'{file_path}/iqtree_output_{version}_{dataset_name}')):
        os.mkdir(f'{file_path}/iqtree_output_{version}_{dataset_name}')
    print(f'*********** Running IQ-TREE for {data_path} ***********')
    cmd_iqtree = f'{iqtree_path} -s {data_path} --redo  -nt 1 -seed 1 --prefix {file_path}/iqtree_output_{version}_{dataset_name}/output'
    mem_before = process_memory()
    start = time.time()
    os.system(cmd_iqtree)
    end = time.time()
    mem_after = process_memory()
    time_elapsed = end - start
    memory_footprint = mem_after - mem_before
    data = f'Time elapsed: {time_elapsed}\n Memory usage:{memory_footprint}'
    with open(f'{file_path}/iqtree_output_{version}_{dataset_name}/{dataset_name}.txt', 'w+') as f:
        f.write(data)
