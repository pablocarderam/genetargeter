import os
from typing import Optional, Dict

import click
import numpy as np
import pandas as pd
from dill import load, dump
from sklearn.ensemble import GradientBoostingRegressor

from .model_comparison import run_models, save_final_model_V3, predict

#
# @click.group()
# def main():
#     pass
#
#
# @main.command()
# @click.option(
#     "--models",
#     help="model e.g. L1, or GP",
#     default="AdaBoost",
#     required=False,
#     type=str,
# )
# @click.option("--test", type=bool, default=False, is_flag=True)
# @click.option("--order", type=int, default=2)
# @click.option("--likelihood", type=str, default="gaussian")  # applies only to GP:
# @click.option("--weighted_degree", type=int, default=3)  # applies only to GP:
# @click.option("--adaboost_learning_rate", type=float, default=0.1)
# @click.option("--adaboost_max_depth", type=int, default=3)
# @click.option("--adaboost_num_estimators", type=int, default=100)
# @click.option("--adaboost_cv", type=bool, is_flag=True, default=False)
# @click.option("--output_model_file", type=str, default="model_file.pkl")
# @click.option("--exp_name", type=str, default=None)
# @click.help_option()
# def model_comparison(
#     models,
#     test,
#     order,
#     likelihood,
#     weighted_degree,
#     adaboost_learning_rate,
#     adaboost_max_depth,
#     adaboost_num_estimators,
#     adaboost_cv,
#     output_model_file,
# ):
#     """command-line version of model_comparison.py
#     \f
#     (see that file for more options?)
#     """
#
#     # store current directory
#     cur_dir = os.getcwd()
#     abspath = os.path.abspath(__file__)
#     dname = os.path.dirname(abspath)
#     # change directory to script directory so that relative paths work
#     os.chdir(dname)
#
#     with open(output_model_file, "rb") as f:
#         _, learn_options = load(f)
#
#     results, all_learn_options = run_models(
#         models=[models],
#         learn_options_set={"final": learn_options},
#         orders=[order],
#         test=test,
#         GP_likelihoods=[likelihood],
#         WD_kernel_degrees=[weighted_degree],
#         adaboost_num_estimators=[adaboost_num_estimators],
#         adaboost_max_depths=[adaboost_max_depth],
#         adaboost_learning_rates=[adaboost_learning_rate],
#         adaboost_CV=adaboost_cv,
#     )
#
#     # if exp_name is None:
#     #     exp_name = list(results.keys())[0]
#
#     os.chdir(cur_dir)
#
#     model = list(results.values())[0][3][0]
#     with open(output_model_file, "wb") as f:
#         dump((model, all_learn_options), f)
#
#
# @main.command()
# @click.option("--filename", help="model filename", required=False, type=str)
# @click.option(
#     "--include_position",
#     help="include position",
#     type=bool,
#     is_flag=True,
#     default=False,
# )
# @click.option(
#     "--learn_options",
#     help="training options set",
#     default=None,
#     required=False,
#     type=dict,
# )
# @click.option("--short_name", default=None, required=False, type=str)
# @click.option("--pam_audit", default=False, is_flag=True, type=bool)
# @click.option("--length_audit", default=False, is_flag=True, type=bool)
# @click.help_option()
# def regenerate_stored_models(
#     filename: str,
#     include_position: bool,
#     learn_options: Optional[Dict[str, str]],
#     short_name: str,
#     pam_audit: bool,
#     length_audit: bool,
# ):
#     save_final_model_V3(
#         filename, include_position, learn_options, short_name, pam_audit, length_audit
#     )
#
#
# @main.command()
# @click.option(
#     "--seq",
#     help="Sequence(s) to score.  If passing several sequences, enclose them "
#     "all in quotations.",
#     required=False,
#     type=str,
# )
# @click.option(
#     "--model",
#     help="model instance to use for prediction",
#     required=False,
#     type=GradientBoostingRegressor,
# )
# @click.option(
#     "--model_file",
#     help=" file name of pickled model to use for prediction",
#     required=False,
#     type=str,
# )
# @click.option(
#     "--pam_audit",
#     help="Enables checking the PAM of each sequence",
#     required=False,
#     is_flag=True,
#     type=bool,
# )
# @click.option(
#     "--length_audit",
#     help="Enables checking the length of each sequences",
#     required=False,
#     is_flag=True,
#     type=bool,
# )
# @click.option(
#     "--learn_options_override",
#     help="a dictionary indicating which learn_options to override",
#     required=False,
#     type=Dict[str, str],
# )
# @click.option(
#     "--print_scores", "-p", help="print scores", required=False, is_flag=True, type=bool
# )
# @click.help_option()
# def score_spacers(
#     seq: str,
#     model: Optional[GradientBoostingRegressor] = None,
#     model_file: Optional[str] = None,
#     pam_audit: bool = True,
#     length_audit: bool = False,
#     learn_options_override: Optional[Dict[str, str]] = None,
#     print_scores: bool = False,
# ) -> np.ndarray:
#
#     if isinstance(seq, str):
#         seq = np.array(seq.split())
#     scores = predict(
#         seq=seq,
#         aa_cut=None,
#         percent_peptide=None,
#         model=model,
#         model_file=model_file,
#         pam_audit=pam_audit,
#         length_audit=length_audit,
#         learn_options_override=learn_options_override,
#     )
#     if print_scores:
#         scores_df = pd.DataFrame({"sequence": seq, "score": scores})
#         print(scores_df)
#     return scores
#
#
# if __name__ == "__main__":
#     main()
