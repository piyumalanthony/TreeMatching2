import pickle

path_1 = '/home/piyumal/PHD/IQTree_experiments/performace_test/seq_fixed_test_2/time_pickle_v1'
path_2 = '/home/piyumal/PHD/IQTree_experiments/performace_test/seq_fixed_test_2/time_pickle_v2'

with open(f'{path_1}/iqtree_v1_all.pkl', 'rb') as f:
    pickle_data = pickle.load(f)

print(pickle_data['time'])


with open(f'{path_2}/iqtree_v2_all.pkl', 'rb') as f:
    pickle_data_v2 = pickle.load(f)

print(pickle_data_v2['time'])