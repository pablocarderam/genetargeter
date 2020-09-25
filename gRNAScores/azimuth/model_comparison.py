from copy import deepcopy
from math import log10
from os import path
from typing import Optional, List, Union, Callable, Dict

import numpy as np
import pandas as pd
from dill import load, dump
from pkg_resources import resource_filename
from sklearn.ensemble import GradientBoostingRegressor

from .features.featurization import featurize_data
from .load_data import get_V3_genes, from_file
from .local_multiprocessing import configure
from .predict import cross_validate
from .util import concatenate_feature_sets, convert_to_thirty_one

DATA_PATH = resource_filename("gRNAScores.azimuth", "saved_models/")


def set_target(learn_options: Dict[str, str],
               classification: bool):
    """Set conditional target learn options

    :param learn_options:typing.Dict[str,str]
    :param classification:bool
    :return:typing.Dict[str,str]
    """


    if "target_name" in learn_options and learn_options["target_name"] is None:
        raise AssertionError("changed it to be automatically set here")
    if not classification:
        learn_options["target_name"] = learn_options["rank-transformed target name"]
        learn_options["training_metric"] = "spearmanr"
        learn_options["ground_truth_label"] = learn_options["target_name"]
    else:
        learn_options["target_name"] = learn_options["binary target name"]
        learn_options["training_metric"] = "AUC"
        learn_options["ground_truth_label"] = learn_options["binary target name"]

    if learn_options["V"] == 3:
        if (
            learn_options["target_name"] != "score_drug_gene_rank"
            and learn_options["target_name"] != "score_drug_gene_threshold"
        ):
            raise AssertionError("cannot use raw scores when mergind data")
        if (
            learn_options["ground_truth_label"] != "score_drug_gene_rank"
            and learn_options["ground_truth_label"] != "score_drug_gene_threshold"
        ):
            raise AssertionError("cannot use raw scores when mergind data")

    return learn_options


def GP_setup(learn_options: Dict[str,str],
             likelihood: str = "gaussian",
             degree: int = 3,
             set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    """

    :param learn_options: typing.Dict[str,str]
    :param likelihood: str
    :param degree:
    :param set_target_fn:
    :return: typing.Dict[str,str]
    """

    learn_options["method"] = "GPy"
    learn_options["kernel degree"] = degree

    if likelihood == "warped":
        learn_options["warpedGP"] = True
    else:
        learn_options["warpedGP"] = False
    learn_options = set_target_fn(learn_options, classification=False)

    return learn_options


def SVC_setup(learn_options: Dict[str, str],
              set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options["method"] = "SVC"
    learn_options = set_target_fn(learn_options, classification=True)

    return learn_options


def L1_setup(learn_options: Dict[str, str],
             set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=False)
    learn_options["method"] = "linreg"
    learn_options["penalty"] = "L1"
    learn_options["feature_select"] = False
    if "alpha" not in learn_options:
        learn_options["alpha"] = np.logspace(-5, log10(1.5e5), num=100)
    learn_options["loss"] = "squared"

    return learn_options


def L2_setup(learn_options: Dict[str, str],
             set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=False)
    learn_options["method"] = "linreg"
    learn_options["penalty"] = "L2"
    learn_options["feature_select"] = False
    if "alpha" not in learn_options:
        learn_options["alpha"] = np.logspace(-5, log10(1.5e5), num=100)
    learn_options["loss"] = "squared"

    return learn_options


def mean_setup(learn_options: Dict[str, str],
               set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=False)
    learn_options["method"] = "mean"
    return learn_options


def random_setup(learn_options,
                 set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=False)
    learn_options["method"] = "random"
    return learn_options


def elasticnet_setup(learn_options: Dict[str, str],
                     set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=False)
    learn_options["method"] = "linreg"
    learn_options["penalty"] = "EN"
    learn_options["feature_select"] = False
    learn_options["loss"] = "squared"
    if "alpha" not in learn_options:
        learn_options["alpha"] = np.array([1e-5 * pow(2, x) for x in range(0, 30)])
    return learn_options


def DNN_setup(learn_options: Dict[str, str],
              set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=False)
    learn_options["method"] = "DNN"
    learn_options["DNN target variable"] = "score"  # 'score_drug_gene_quantized'
    # learn_options['DNN architecture'] = (119, 10, 10, 10, 2)
    return learn_options


def RF_setup(learn_options: Dict[str, str],
             set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=False)
    learn_options["method"] = "RandomForestRegressor"
    return learn_options


def doench_setup(learn_options: Dict[str, str],
                 set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=True)
    learn_options["method"] = "doench"
    return learn_options


def sgrna_from_doench_setup(learn_options: Dict[str, str],
                            set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=False)
    learn_options["method"] = "sgrna_from_doench"
    return learn_options


def linreg_setup(learn_options: Dict[str, str],
                 set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options["method"] = "linreg"
    learn_options["penalty"] = None
    learn_options["feature_select"] = False
    if "alpha" not in learn_options:
        learn_options["alpha"] = np.array([0.0])
    learn_options["loss"] = "squared"
    learn_options = set_target_fn(learn_options, classification=False)

    return learn_options


def logregL1_setup(learn_options: Dict[str, str],
                   set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=True)
    learn_options["method"] = "logregL1"
    learn_options["penalty"] = "L1"
    learn_options["feature_select"] = False
    if "alpha" not in learn_options:
        learn_options["alpha"] = np.logspace(-5, log10(1.5e5), num=100)
    if "fit_intercept" not in learn_options:
        learn_options["fit_intercept"] = True
    return learn_options


def LASSOs_ensemble_setup(learn_options: Dict[str, str],
                          set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=False)
    learn_options["method"] = "lasso_ensemble"
    learn_options["penalty"] = "L1"
    learn_options["feature_select"] = False
    if "alpha" not in learn_options:
        learn_options["alpha"] = np.logspace(-5, log10(1.5e5), num=100)
    learn_options["loss"] = "squared"

    return learn_options


def xu_et_al_setup(learn_options: Dict[str, str],
                   set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target):
    learn_options = set_target_fn(learn_options, classification=True)
    learn_options["method"] = "xu_et_al"

    return learn_options


def adaboost_setup(
    learn_options: Dict[str, str],
    num_estimators=100,
    max_depth=3,
    learning_rate=0.1,
    set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target,
    model="AdaBoost",
):
    """
    """
    learn_options = set_target_fn(
        learn_options, classification=(model == "AdaBoostClassifier")
    )
    if model == "AdaBoost":
        learn_options["method"] = "AdaBoostRegressor"
    elif model == "AdaBoostClassifier":
        learn_options["method"] = "AdaBoostClassifier"
    else:
        raise Exception("model must be either AdaBoost or AdaBoost Classifier")
    learn_options["adaboost_version"] = "python"  # "R" or "python"

    if "adaboost_loss" not in learn_options and model == "AdaBoostRegressor":
        learn_options[
            "adaboost_loss"
        ] = (
            "ls"
        )  # alternatives: "lad", "huber", "quantile", see scikit docs for details
    if "adaboost_alpha" not in learn_options:
        learn_options[
            "adaboost_alpha"
        ] = 0.5  # this parameter is only used by the huber and quantile loss functions.

    if not learn_options["adaboost_CV"]:
        learn_options["adaboost_learning_rate"] = learning_rate
        learn_options["adaboost_n_estimators"] = num_estimators
        learn_options["adaboost_max_depth"] = max_depth
    else:
        learn_options["adaboost_n_estimators"] = num_estimators

    return learn_options


def shared_setup(learn_options: Dict[str, str],
                 order,
                 test):
    if "num_proc" not in learn_options:
        learn_options["num_proc"] = None
    if "num_thread_per_proc" not in learn_options:
        learn_options["num_thread_per_proc"] = None

    num_proc = configure(
        test=test,
        num_proc=learn_options["num_proc"],
        num_thread_per_proc=learn_options["num_thread_per_proc"],
    )
    if num_proc > 1:
        learn_options["num_proc"] = num_proc - 1
    else:
        learn_options["num_proc"] = num_proc

    learn_options["order"] = order  # gets used many places in code, not just here

    if "cv" not in learn_options:
        # if no CV preference is specified, use leave-one-gene-out
        learn_options["cv"] = "gene"

    if "normalize_features" not in learn_options:
        # if no CV preference is specified, use leave-one-gene-out
        learn_options["normalize_features"] = True

    if "weighted" not in learn_options:
        learn_options["weighted"] = None

    if "all pairs" not in learn_options:
        learn_options["all pairs"] = False

    if "include_known_pairs" not in learn_options:
        learn_options["include_known_pairs"] = False

    if "include_gene_guide_feature" not in learn_options:
        learn_options[
            "include_gene_guide_feature"
        ] = 0  # used as window size, so 0 is none

    # these should default to true to match experiments before they were options:
    if "gc_features" not in learn_options:
        learn_options["gc_features"] = True
    if "nuc_features" not in learn_options:
        learn_options["nuc_features"] = True

    if "train_genes" not in learn_options:
        learn_options["train_genes"] = None
    if "test_genes" not in learn_options:
        learn_options["test_genes"] = None

    if "num_proc" not in learn_options:
        learn_options["num_proc"] = None
    if "num_thread_per_proc" not in learn_options:
        learn_options["num_thread_per_proc"] = None

    if "seed" not in learn_options:
        learn_options["seed"] = 1

    if "flipV1target" not in learn_options:
        learn_options["flipV1target"] = False

    if "num_genes_remove_train" not in learn_options:
        learn_options["num_genes_remove_train"] = None

    if "include_microhomology" not in learn_options:
        learn_options["include_microhomology"] = False

    if "algorithm_hyperparam_search" not in learn_options:
        learn_options[
            "algorithm_hyperparam_search"
        ] = "grid"  # other options is bo for bayesian optimization

    return num_proc


def setup(
    test=False,
    order=1,
    learn_options: Optional[Dict[str, str]] = None,
    data_file=None,
    pam_audit=True,
    length_audit=True,
):
    print(f"learn_options: {learn_options}")
    num_proc = shared_setup(learn_options, order, test)

    if "testing_non_binary_target_name" not in learn_options:
        raise AssertionError(
            "need this in order to get metrics though used to be not needed, so you may newly see this error"
        )
    if learn_options["testing_non_binary_target_name"] not in ["ranks", "raw", "thrs"]:
        raise Exception(
            'learn_options["testing_non_binary_target_name"] must be in ["ranks", "raw", "thrs"]'
        )

    x_df, Y, gene_position, target_genes = from_file(
        data_file=data_file, data_file2=None, learn_options=learn_options
    )
    learn_options["all_genes"] = target_genes

    if test:
        learn_options["order"] = 1

    if (
        "convert_30mer_to_31mer" in learn_options
        and learn_options["convert_30mer_to_31mer"] is True
    ):
        print(
            "WARNING!!! converting 30 mer to 31 mer (and then cutting off first nucleotide to go "
            "back to 30mer with a right shift)"
        )
        for i in range(x_df.shape[0]):
            x_df["30mer"].iloc[i] = convert_to_thirty_one(
                x_df.iloc[i]["30mer"], x_df.index.values[i][1], x_df.iloc[i]["Strand"]
            )

        x_df["30mer"] = x_df["30mer"].apply(
            lambda x: x[1:]
        )  # chop the first nucleotide

    if (
        "left_right_guide_ind" in learn_options
        and learn_options["left_right_guide_ind"] is not None
    ):
        seq_start, seq_end, expected_length = learn_options["left_right_guide_ind"]
        if len(x_df["30mer"].values[0]) != expected_length:
            raise AssertionError("incorrect spacer length")
        x_df["30mer"] = x_df["30mer"].apply(lambda seq: seq[seq_start:seq_end])

    feature_sets = featurize_data(
        x_df,
        learn_options,
        Y,
        gene_position,
        pam_audit=pam_audit,
        length_audit=length_audit,
    )
    np.random.seed(learn_options["seed"])

    return Y, feature_sets, target_genes, learn_options, num_proc


def run_models(
    models,
    orders,
    GP_likelihoods: Optional[List[str]] = None,
    WD_kernel_degrees: Union[List[int], int] = 3,
    adaboost_learning_rates: Union[List[float], float] = 0.1,
    adaboost_num_estimators: Union[List[float], float] = 100,
    adaboost_max_depths: Union[List[int], int] = 3,
    learn_options_set: Optional[dict] = None,
    test: bool = False,
    adaboost_CV: bool = True,
    setup_function: Callable = setup,
    set_target_fn: Callable[[Dict[str, str], bool], Dict[str, str]] = set_target,
    pam_audit: bool = True,
    length_audit: bool = True,
):
    """
    CV is set to false if want to train a final model and not cross-validate, but it goes in to what
    looks like cv code
    """

    if isinstance(WD_kernel_degrees, int):
        WD_kernel_degrees = [WD_kernel_degrees]
    if isinstance(adaboost_learning_rates, float):
        adaboost_learning_rates = [adaboost_learning_rates]
    if isinstance(adaboost_num_estimators, float):
        adaboost_num_estimators = [adaboost_num_estimators]
    if isinstance(adaboost_max_depths, int):
        adaboost_max_depths = [adaboost_max_depths]

    if GP_likelihoods is None:
        GP_likelihoods = ["gaussian", "warped"]

    results = {}
    if learn_options_set is None:
        raise AssertionError("need to specify learn_options_set")
    all_learn_options = {}

    # shorten so easier to display on graphs
    feat_models_short = {
        "L1": "L1",
        "L2": "L2",
        "elasticnet": "EN",
        "linreg": "LR",
        "RandomForest": "RF",
        "AdaBoost": "AB",
        "AdaBoostClassifier": "ABClass",
        "doench": "doench",
        "logregL1": "logregL1",
        "sgrna_from_doench": "sgrna_from_doench",
        "SVC": "SVC",
        "xu_et_al": "xu_et_al",
    }

    if not adaboost_CV:
        print(
            "Received option adaboost_CV=False, so I'm training using all of the data"
        )
        if len(learn_options_set) != 1:
            raise AssertionError(
                f"When CV is False, only 1 set of learn options is allowed.  Instead, "
                f"{len(learn_options_set)} were given.\n"
            )
        if len(models) != 1:
            raise AssertionError("when CV is False, only 1 model is allowed")

    for learn_options_str in learn_options_set:
        # these options get augmented in setup
        partial_learn_opt = learn_options_set[learn_options_str]
        # if the model requires encoded features
        for model in models:
            # models requiring explicit featurization
            if model in feat_models_short:
                for order in orders:
                    print(f"running {model}, order {order} for {learn_options_str}")

                    Y, feature_sets, _, learn_options, _ = setup_function(
                        test=test,
                        order=order,
                        learn_options=partial_learn_opt,
                        pam_audit=pam_audit,
                        length_audit=length_audit,
                    )  # TODO precompute features for all orders, as this is repeated for each model

                    if model == "L1":
                        learn_options_model = L1_setup(
                            deepcopy(learn_options), set_target_fn=set_target_fn
                        )
                    elif model == "L2":
                        learn_options_model = L2_setup(
                            deepcopy(learn_options), set_target_fn=set_target_fn
                        )
                    elif model == "elasticnet":
                        learn_options_model = elasticnet_setup(
                            deepcopy(learn_options), set_target_fn=set_target_fn
                        )
                    elif model == "linreg":
                        learn_options_model = linreg_setup(
                            deepcopy(learn_options), set_target_fn=set_target_fn
                        )
                    elif model == "logregL1":
                        learn_options_model = logregL1_setup(
                            deepcopy(learn_options), set_target_fn=set_target_fn
                        )
                    elif model == "RandomForest":
                        learn_options_model = RF_setup(
                            deepcopy(learn_options), set_target_fn=set_target_fn
                        )
                    elif model == "SVC":
                        learn_options_model = SVC_setup(
                            deepcopy(learn_options), set_target_fn=set_target_fn
                        )
                    elif model == "doench":
                        learn_options_model = doench_setup(
                            deepcopy(learn_options), set_target_fn=set_target_fn
                        )
                    elif model == "sgrna_from_doench":
                        learn_options_model = sgrna_from_doench_setup(
                            deepcopy(learn_options), set_target_fn=set_target_fn
                        )
                    elif model == "xu_et_al":
                        learn_options_model = xu_et_al_setup(
                            deepcopy(learn_options), set_target_fn=set_target_fn
                        )
                    elif model == "AdaBoost" or "AdaBoostClassifier":
                        for learning_rate in adaboost_learning_rates:
                            for num_estimators in adaboost_num_estimators:
                                for max_depth in adaboost_max_depths:
                                    learn_options_model = adaboost_setup(
                                        deepcopy(learn_options),
                                        learning_rate=learning_rate,
                                        num_estimators=num_estimators,
                                        max_depth=max_depth,
                                        set_target_fn=set_target_fn,
                                        model=model,
                                    )
                        model_string = (
                            feat_models_short[model]
                            + f"_or{learn_options_set[learn_options_str]['order']}_"
                            f"md{max_depth}_lr{learning_rate:.2f}_"
                            f"n{num_estimators}_{learn_options_str}"
                        )
                    if model != "AdaBoost":
                        model_string = (
                            feat_models_short[model]
                            + f"_ord{learn_options_set[learn_options_str]['order']}_"
                            f"{learn_options_str}"
                        )

                    results[model_string] = cross_validate(
                        Y,
                        feature_sets,
                        learn_options=learn_options_model,
                        TEST=test,
                        CV=adaboost_CV,
                    )

                    all_learn_options[model_string] = learn_options_model
                    print(f"finished computing {model}")
            # if the model doesn't require explicit featurization
            else:
                if setup_function != setup:
                    raise AssertionError("not yet modified to handle this")
                print(f"running {model} for {learn_options_str}")
                Y, feature_sets, _, learn_options, _ = setup(
                    test=test,
                    order=1,
                    learn_options=partial_learn_opt,
                    pam_audit=pam_audit,
                    length_audit=length_audit,
                )
                if model == "mean":
                    learn_options_model = mean_setup(deepcopy(learn_options))
                elif model == "random":
                    learn_options_model = random_setup(deepcopy(learn_options))
                elif model == "DNN":
                    learn_options_model = DNN_setup(deepcopy(learn_options))
                elif model == "GP":
                    for likelihood in GP_likelihoods:
                        for degree in WD_kernel_degrees:
                            learn_options_model = GP_setup(
                                deepcopy(learn_options),
                                likelihood=likelihood,
                                degree=degree,
                            )
                            model_string = f"{model}_{likelihood}_degree{degree}_{learn_options_str}"
                            results[model_string] = cross_validate(
                                Y,
                                feature_sets,
                                learn_options=learn_options_model,
                                TEST=test,
                                CV=adaboost_CV,
                            )

                else:
                    raise NotImplementedError(f"model {model} not supported")

                # "GP" already calls cross_validate() and has its own model_string, so skip this.
                if model != "GP":
                    model_string = model + f"_{learn_options_str}"
                    results[model_string] = cross_validate(
                        Y,
                        feature_sets,
                        learn_options=learn_options_model,
                        TEST=test,
                        CV=adaboost_CV,
                    )

            all_learn_options[model_string] = learn_options_model

    return results, all_learn_options


def save_final_model_V3(
    filename: str = None,
    include_position: bool = True,
    learn_options: Optional[Dict[str, str]] = None,
    short_name: str = "final",
    pam_audit: bool = True,
    length_audit: bool = True,
):
    """
    run_models(produce_final_model=True) is what saves the model
    """
    test = False
    if filename is None:
        raise AssertionError("need to provide filename to save final model")

    if learn_options is None:
        learn_options = {
            "V": 3,
            "train_genes": get_V3_genes(),
            "test_genes": get_V3_genes(),
            "testing_non_binary_target_name": "ranks",
            "include_pi_nuc_feat": True,
            "gc_features": True,
            "nuc_features": True,
            "include_gene_position": False,
            "include_NGGX_interaction": True,
            "include_Tm": True,
            "include_strand": False,
            "include_gene_feature": False,
            "include_gene_guide_feature": 0,
            "extra pairs": False,
            "weighted": None,
            "training_metric": "spearmanr",
            "NDGC_k": 10,
            "cv": "gene",
            "include_gene_effect": False,
            "include_drug": False,
            "include_sgRNAscore": False,
            "adaboost_loss": "ls",  # main 'ls', alternatives: 'lad', 'huber', 'quantile', see scikit docs for details
            "adaboost_alpha": 0.5,  # this parameter is only used by the huber and quantile loss functions.
            "normalize_features": False,
            "adaboost_CV": False,
        }
        if include_position:
            learn_options["include_gene_position"] = True

    learn_options_set = {short_name: learn_options}
    results, _ = run_models(
        ["AdaBoost"],
        orders=[2],
        adaboost_learning_rates=[0.1],
        adaboost_max_depths=[3],
        adaboost_num_estimators=[100],
        learn_options_set=learn_options_set,
        test=test,
        adaboost_CV=False,
        pam_audit=pam_audit,
        length_audit=length_audit,
    )
    model = list(results.values())[0][3][0]

    with open(filename, "wb") as f:
        dump((model, learn_options), f, -1)

    return model


def predict(
    seq: np.ndarray,
    aa_cut: Optional[np.ndarray] = None,
    percent_peptide: Optional[np.ndarray] = None,
    model: Optional[GradientBoostingRegressor] = None,
    model_file: Optional[str] = None,
    pam_audit: bool = True,
    length_audit: bool = False,
    learn_options_override: Optional[Dict[str, str]] = None,
    verbose: bool = False,
):
    """

    Parameters
    ----------
    seq : :class:np.ndarray
        numpy array of 30 nt sequences.
    aa_cut : numpy array of amino acid cut positions (optional).
    percent_peptide : numpy array of percent peptide (optional).
    model : model instance to use for prediction (optional).
    model_file : file name of pickled model to use for prediction (optional).
    pam_audit : check PAM of each sequence.
    length_audit : check length of each sequence.
    learn_options_override : a dictionary indicating which learn_options to override (optional).
    verbose : bool
        display extra information

    Return
    ------
    :class:`~np.array`
    """

    if not isinstance(seq, np.ndarray):
        raise AssertionError("Please ensure seq is a numpy array")
    if len(seq[0]) <= 0:
        raise AssertionError("Make sure that seq is not empty")
    if not isinstance(seq[0], str):
        raise AssertionError(
            f"Please ensure input sequences are in string format, i.e. 'AGAG' "
            f"rather than ['A' 'G' 'A' 'G'] or alternate representations"
        )

    if aa_cut is not None:
        if len(aa_cut) <= 0:
            raise AssertionError("Make sure that aa_cut is not empty")
        if not isinstance(aa_cut, np.ndarray):
            raise AssertionError("Please ensure aa_cut is a numpy array")
        if not np.all(np.isreal(aa_cut)):
            raise AssertionError("amino-acid cut position needs to be a real number")

    if percent_peptide is not None:
        if len(percent_peptide) <= 0:
            raise AssertionError("Make sure that percent_peptide is not empty")
        if not isinstance(percent_peptide, np.ndarray):
            raise AssertionError("Please ensure percent_peptide is a numpy array")
        if not np.all(np.isreal(percent_peptide)):
            raise AssertionError("percent_peptide needs to be a real number")

    if model_file is None:
        if np.any(percent_peptide == -1) or (
            percent_peptide is None and aa_cut is None
        ):
            if verbose:
                print("No model file specified, using V3_model_nopos")
            model_name = "V3_model_nopos.pickle"
        else:
            if verbose:
                print("No model file specified, using V3_model_full")
            model_name = "V3_model_full.pickle"

        model_file = path.join(DATA_PATH, model_name)

    if model is None:
        with open(model_file, "rb") as f:
            model, learn_options = load(f)
    else:
        model, learn_options = model

    learn_options["V"] = 2

    learn_options = override_learn_options(learn_options_override, learn_options)

    x_df = pd.DataFrame(
        columns=["30mer", "Strand"],
        data=list(zip(seq, ["NA" for x in range(len(seq))])),
    )

    if np.all(percent_peptide != -1) and (
        percent_peptide is not None and aa_cut is not None
    ):
        gene_position = pd.DataFrame(
            columns=["Percent Peptide", "Amino Acid Cut position"],
            data=list(zip(percent_peptide, aa_cut)),
        )
    else:
        gene_position = pd.DataFrame(
            columns=["Percent Peptide", "Amino Acid Cut position"],
            data=list(zip(np.ones(seq.shape[0]) * -1, np.ones(seq.shape[0]) * -1)),
        )

    feature_sets = featurize_data(
        x_df,
        learn_options,
        pd.DataFrame(),
        gene_position,
        pam_audit=pam_audit,
        length_audit=length_audit,
    )
    inputs, *_ = concatenate_feature_sets(feature_sets)

    # call to scikit-learn, returns a vector of predicted values
    preds = model.predict(inputs)

    # also check that predictions are not 0/1 from a classifier.predict()
    # (instead of predict_proba() or decision_function())

    if np.all([True if pr in (0, 1) else False for pr in np.unique(preds)]):
        raise AssertionError("model returned only 0s and 1s")
    return preds


def override_learn_options(learn_options_override: Optional[Dict[str, str]],
                           learn_options: Dict[str, str]):
    """
    override all keys seen in learn_options_override to alter learn_options
    """
    if learn_options_override is not None:
        for k in learn_options_override:
            learn_options[k] = learn_options_override[k]
    return learn_options


def fill_learn_options(learn_options_used_to_fill: Optional[Dict[str, str]],
                       learn_options_with_possible_missing: Dict[str, str]):
    """
    only fill in keys that are missing from learn_options from learn_options_fill
    """
    if learn_options_used_to_fill is not None:
        for k in learn_options_used_to_fill:
            if k not in learn_options_with_possible_missing:
                learn_options_with_possible_missing[k] = learn_options_used_to_fill[k]
    return learn_options_with_possible_missing


def write_results(predictions,
                  file_to_predict: str):
    newfile = file_to_predict.replace(".csv", ".pred.csv")
    data = pd.read_csv(file_to_predict)
    data["predictions"] = predictions
    data.to_csv(newfile)
    print(f"wrote results to {newfile}")
    return data, newfile


if __name__ == "__main__":

    save_final_model_V3(
        filename="saved_models/V3_model_nopos.pickle", include_position=False
    )
    save_final_model_V3(
        filename="saved_models/V3_model_full.pickle", include_position=True
    )

    learn_options = {
        "V": 3,
        "train_genes": get_V3_genes(),
        "test_genes": get_V3_genes(),
        "target_name": "score_drug_gene_rank",
        "testing_non_binary_target_name": "ranks",
        "include_pi_nuc_feat": True,
        "gc_features": True,
        "nuc_features": True,
        "include_gene_position": True,
        "include_NGGX_interaction": True,
        "include_Tm": True,
        "include_strand": False,
        "include_gene_feature": False,
        "include_gene_guide_feature": 0,
        "extra pairs": False,
        "weighted": None,
        "training_metric": "spearmanr",
        "NDGC_k": 10,
        "cv": "gene",
        "include_gene_effect": False,
        "include_drug": False,
        "include_sgRNAscore": False,
        "adaboost_loss": "ls",
        # main "ls", alternatives: "lad", "huber", "quantile", see scikit docs for details
        "adaboost_alpha": 0.5,  # this parameter is only used by the huber and quantile loss functions.
        "adaboost_CV": False,
    }

    learn_options_set = {"post bug fix": learn_options}
