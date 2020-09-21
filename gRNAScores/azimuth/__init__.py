# from . import features
# from . import metrics
# from . import models
# from .__main__ import (
#     score_spacers,
#     regenerate_stored_models,
#     model_comparison
# )
# # from .cli_run_model import main
# from .corrstats import rz_ci, rho_rxy_rxz, dependent_corr, independent_corr
# from .load_data import (
#     from_custom_file,
#     from_file,
#     set_V2_target_names,
#     combine_organisms,
#     read_V1_data,
#     read_V2_data,
#     merge_all,
#     mergeV1_V2,
#     get_V1_genes,
#     get_V2_genes,
#     get_V3_genes,
#     get_human_genes,
#     get_mouse_genes,
# )
# from .local_multiprocessing import configure
# from .metrics import (
#     mean_reciprocal_rank,
#     r_precision,
#     precision_at_k,
#     average_precision,
#     mean_average_precision,
#     dcg_at_k,
#     ndcg_at_k,
#     ndcg_at_k_ties,
#     dcg_helper,
#     dcg_at_k_ties,
#     get_discount_factors,
#     rank_data,
#     dcg_alt,
#     ndcg_alt,
#     ndcg_at_k_swap_perm_test,
# )
# from .model_comparison import (
#     set_target,
#     GP_setup,
#     SVC_setup,
#     L1_setup,
#     L2_setup,
#     mean_setup,
#     random_setup,
#     elasticnet_setup,
#     DNN_setup,
#     RF_setup,
#     doench_setup,
#     sgrna_from_doench_setup,
#     linreg_setup,
#     logregL1_setup,
#     LASSOs_ensemble_setup,
#     xu_et_al_setup,
#     adaboost_setup,
#     shared_setup,
#     setup,
#     run_models,
#     save_final_model_V3,
#     predict,
#     override_learn_options,
#     fill_learn_options,
#     write_results,
# )
# from .predict import (
#     fill_in_truth_and_predictions,
#     construct_filename,
#     extract_fpr_tpr_for_fold,
#     extract_NDCG_for_fold,
#     extract_spearman_for_fold,
#     get_train_test,
#     cross_validate,
# )
# from .util import (
#     get_thirty_one_mer_data,
#     convert_to_thirty_one,
#     concatenate_feature_sets,
#     spearmanr_nonan,
#     impute_gene_position,
#     get_gene_sequence,
#     get_ranks,
#     get_data,
#     extract_feature_from_model,
#     extract_feature_from_model_sum,
#     feature_importances,
# )
#
# __author__ = ("Nicolo Fusi", "Jennifer Listgarten", "Miles Smith")
# __email__ = (
#     "fusi@microsoft.com",
#     "jennl@microsoft.com",
#     "mileschristiansmith@gmail.com",
# )
