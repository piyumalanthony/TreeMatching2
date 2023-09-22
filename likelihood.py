import numpy as np


def calculate_approx_likelihood(br_lengths, gradients, hessians, iqtree_idx, replicates=100):
    iqtree_br = np.array(br_lengths[0])
    baseml_br = np.array(br_lengths[1])

    iqtree_gradients = np.array(gradients[0])
    baseml_gradients = np.array(gradients[1])

    iqtree_hessian = np.array(hessians[0])
    baseml_hessian = np.array(hessians[1])

    iq_tree_mle_estimates = []
    baseml_mle_estimates = []
    iq_tree_mle_estimates_delta = []
    baseml_mle_estimates_delta = []

    for i in range(replicates):
        noise = np.random.default_rng().uniform(0.0001, 0.01, len(iqtree_br))
        baseml_br_lengths = baseml_br + noise

        iqtree_noise = []
        noise_b = list(noise)
        for index in iqtree_idx:
            iqtree_noise.append(noise_b[index])
        iqtree_noise = np.array(iqtree_noise)
        iq_br_lengths = iqtree_br + iqtree_noise
        iqtree_mle = np.dot(iqtree_gradients, iq_br_lengths) + np.dot(iq_br_lengths,
                                                                      iqtree_hessian.dot(iq_br_lengths)) / 2
        iqtree_mle_noise = np.dot(iqtree_gradients, iqtree_noise) + np.dot(iqtree_noise,
                                                                           iqtree_hessian.dot(iqtree_noise)) / 2
        baseml_mle = np.dot(baseml_gradients, baseml_br_lengths) + np.dot(baseml_br_lengths,
                                                                          baseml_hessian.dot(baseml_br_lengths)) / 2
        baseml_mle_noise = np.dot(baseml_gradients, noise) + np.dot(noise, baseml_hessian.dot(noise)) / 2
        iq_tree_mle_estimates.append(iqtree_mle)
        iq_tree_mle_estimates_delta.append(iqtree_mle_noise)
        baseml_mle_estimates.append(baseml_mle)
        baseml_mle_estimates_delta.append(baseml_mle_noise)

    return iq_tree_mle_estimates, baseml_mle_estimates, iq_tree_mle_estimates_delta, baseml_mle_estimates_delta
