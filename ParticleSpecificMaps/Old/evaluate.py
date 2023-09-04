
from os.path import join

import matplotlib.pyplot as plt



def evaluate(inspectors, config):

    map_params = {}
    insp = inspectors[0]

    losses = insp.get_losses()

    # find minimum index
    best_trial = insp._trials_complete[index]

    best_cwd = insp.get_annotations_per_trial("cwd")[index]

    for i, insp in enumerate(inspectors):
        figure, _ = insp.plot_loss_feature_history(map_params=map_params, n_most_important=20)
        figure.tight_layout()
        figure.savefig(f"loss_feature_history_{i}.png")
        plt.close(figure)

        figure, _ = insp.plot_importance(map_params=map_params, n_most_important=50)
        figure.tight_layout()
        figure.savefig(f"importance_parameters_{i}.png")
        plt.close(figure)

        figure, _ = insp.plot_parallel_coordinates(map_params=map_params)
        figure.savefig(f"parallel_coordinates_{i}.png")
        plt.close(figure)

        figure, _ = insp.plot_slices(map_params=map_params)
        figure.savefig(f"slices_{i}.png")
        plt.close(figure)

        figure, _ = insp.plot_correlations(map_params=map_params)
        figure.savefig(f"parameter_correlations_{i}.png")
        plt.close(figure)

        figure, _ = insp.plot_pairwise_scatter(map_params=map_params)
        figure.savefig(f"pairwise_scatter_{i}.png")
        plt.close(figure)

    return True