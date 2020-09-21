from copy import deepcopy
from multiprocessing import Pool
from time import time

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder

from .metrics import ndcg_at_k_ties
from .models import DNN, GP, baselines, ensembles, regression
from .util import spearmanr_nonan, concatenate_feature_sets


def fill_in_truth_and_predictions(
    truth, predictions, fold, y_all, y_pred, learn_options, test
):
    truth[fold]["ranks"] = np.hstack(
        (
            truth[fold]["ranks"],
            y_all[learn_options["rank-transformed target name"]].values[test].flatten(),
        )
    )

    truth[fold]["thrs"] = np.hstack(
        (
            truth[fold]["thrs"],
            y_all[learn_options["binary target name"]].values[test].flatten(),
        )
    )

    if "raw_target_name" in learn_options:
        truth[fold]["raw"] = np.hstack(
            (
                truth[fold]["raw"],
                y_all[learn_options["raw target name"]].values[test].flatten(),
            )
        )

    predictions[fold] = np.hstack((predictions[fold], y_pred.flatten()))

    return truth, predictions


def construct_filename(learn_options, TEST):
    if "V" in learn_options:
        filename = f"V{learn_options['V']}"
    else:
        filename = "offV1"

    if TEST:
        filename = "TEST."

    filename += learn_options["method"]
    filename += f'.order{learn_options["order"]}'

    filename += learn_options["target_name"]
    if learn_options["method"] == "linreg" and learn_options["penalty"] is not None:
        filename += f".{learn_options['penalty']}"
    filename += "." + learn_options["cv"]

    if learn_options["training_metric"] == "NDCG":
        filename += f".NDGC_{learn_options['NDGC_k']}"
    elif learn_options["training_metric"] == "AUC":
        filename += ".AUC"
    elif learn_options["training_metric"] == "spearmanr":
        filename += ".spearman"

    print(f"filename = {filename}")
    return filename


def extract_fpr_tpr_for_fold(aucs, y_binary, test, y_pred):
    if len(np.unique(y_binary)) > 2:
        raise AssertionError("if using AUC need binary targets")
    fpr, tpr, _ = roc_curve(y_binary[test], y_pred)
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)


def extract_NDCG_for_fold(metrics, y_ground_truth, test, y_pred, learn_options):
    NDCG_fold = ndcg_at_k_ties(
        y_ground_truth[test].flatten(), y_pred.flatten(), learn_options["NDGC_k"]
    )
    metrics.append(NDCG_fold)


def extract_spearman_for_fold(metrics, y_ground_truth, test, y_pred):
    spearman = np.nan_to_num(
        spearmanr(y_ground_truth[test].flatten(), y_pred.flatten())[0]
    )
    if np.isnan(spearman):
        raise AssertionError("found nan spearman")
    metrics.append(spearman)


def get_train_test(test_gene, y_all, train_genes=None):
    # this is a bit convoluted because the train_genes+test_genes may not add up to all genes
    # for e.g. when we load up V3, but then use only V2, etc.

    not_test = y_all.index.get_level_values("Target gene").values != test_gene

    if train_genes is not None:
        in_train_genes = np.zeros(not_test.shape, dtype=bool)
        for t_gene in train_genes:
            in_train_genes = np.logical_or(
                in_train_genes,
                (y_all.index.get_level_values("Target gene").values == t_gene),
            )
        train = np.logical_and(not_test, in_train_genes)
    else:
        train = not_test
    # y_all['test'] as to do with extra pairs in V2
    if test_gene == "dummy":
        test = train
    else:
        test = y_all.index.get_level_values("Target gene").values == test_gene

    # convert to indices
    test = np.where(test)[0]
    train = np.where(train)[0]
    return train, test


def cross_validate(y_all, feature_sets, learn_options=None, TEST=False, CV=True):
    # feature_sets is a dictionary of "set name" to pandas.DataFrame
    # one set might be single-nucleotide, position-independent features of order X, for e.g.
    # Method: "GPy" or "linreg"
    # Metric: NDCG (learning to rank metric, Normalized Discounted Cumulative Gain); AUC
    # Output: cv_score_median, gene_rocs
    # When CV=False, it trains on everything (and tests on everything, just to fit the code)

    # print(f"range of y_all is [{np.min(y_all[learn_options['target_name']].values)}, "
    #       f"{np.max(y_all[learn_options['target_name']].values)}]")

    allowed_methods = [
        "GPy",
        "linreg",
        "AdaBoostRegressor",
        "AdaBoostClassifier",
        "DecisionTreeRegressor",
        "RandomForestRegressor",
        "ARDRegression",
        "mean",
        "random",
        "DNN",
        "lasso_ensemble",
        "doench",
        "logregL1",
        "sgrna_from_doench",
        "SVC",
        "xu_et_al",
    ]

    if learn_options["method"] not in allowed_methods:
        raise AssertionError("invalid method: {learn_options['method']}")
    if (
        learn_options["method"] != "linreg" or learn_options["penalty"] != "L2"
    ) and learn_options["weighted"] is not None:
        raise AssertionError(
            f"{learn_options['method']} {learn_options['weighted']} weighted "
            f"only works with linreg L2 right now"
        )

    # construct filename from options
    filename = construct_filename(learn_options, TEST)

    print("Cross-validating genes...")
    t2 = time()

    y = np.array(y_all[learn_options["target_name"]].values[:, None], dtype=np.float64)

    # concatenate feature sets in to one nparray, and get dimension of each
    inputs, _, dimsum, feature_names = concatenate_feature_sets(feature_sets)

    if not CV:
        if learn_options["cv"] != "gene":
            raise AssertionError(
                "Must use gene-CV when CV is False (I need to use all of the "
                "genes and stratified complicates that)"
            )

    # set-up for cross-validation
    # for outer loop, the one Doench et al use genes for
    if learn_options["cv"] == "stratified":
        if "extra_pairs" in learn_options or learn_options["extra pairs"]:
            raise AssertionError(
                "can't use extra pairs with stratified CV need to figure out how to properly "
                "account for genes affected by two drugs"
            )
        label_encoder = LabelEncoder()
        label_encoder.fit(y_all["Target gene"].values)
        gene_classes = label_encoder.transform(y_all["Target gene"].values)
        if "n_folds" in learn_options:
            n_splits = learn_options["n_folds"]
        elif (
            learn_options["train_genes"] is not None
            and learn_options["test_genes"] is not None
        ):
            n_splits = len(learn_options["test_genes"])
        else:
            n_splits = len(learn_options["all_genes"])

        skf = StratifiedKFold(n_splits=n_splits, shuffle=True)
        cv = skf.split(np.zeros(len(gene_classes), dtype=np.bool), gene_classes)
        fold_labels = [f"fold{i:d}" for i in range(1, n_splits + 1)]
        if learn_options["num_genes_remove_train"] is not None:
            raise NotImplementedError
    elif learn_options["cv"] == "gene":
        cv = []

        if not CV:
            train_test_tmp = get_train_test(
                "dummy", y_all
            )  # get train, test split using a dummy gene
            # train_tmp, test_tmp = train_test_tmp
            # not a typo, using training set to test on as well, just for this case.
            # Test set is not used for internal cross-val, etc. anyway.
            # train_test_tmp = (train_tmp, train_tmp)
            cv.append(train_test_tmp)
            fold_labels = ["dummy_for_no_cv"]  # learn_options['all_genes']

        elif (
            learn_options["train_genes"] is not None
            and learn_options["test_genes"] is not None
        ):
            if (
                learn_options["train_genes"] is None
                or learn_options["test_genes"] is None
            ):
                raise AssertionError("use both or neither")
            for i, gene in enumerate(learn_options["test_genes"]):
                cv.append(get_train_test(gene, y_all, learn_options["train_genes"]))
            fold_labels = learn_options["test_genes"]
            # if train and test genes are seperate, there should be only one fold
            # train_test_disjoint = set.isdisjoint(set(learn_options["train_genes"].tolist()),
            #                                      set(learn_options["test_genes"].tolist()))

        else:
            for i, gene in enumerate(learn_options["all_genes"]):
                train_test_tmp = get_train_test(gene, y_all)
                cv.append(train_test_tmp)
            fold_labels = learn_options["all_genes"]

        if learn_options["num_genes_remove_train"] is not None:
            for i, (train, test) in enumerate(cv):
                unique_genes = np.random.permutation(
                    np.unique(np.unique(y_all["Target gene"][train]))
                )
                genes_to_keep = unique_genes[
                    0 : len(unique_genes) - learn_options["num_genes_remove_train"]
                ]
                filtered_train = []
                for j, gene in enumerate(y_all["Target gene"]):
                    if j in train and gene in genes_to_keep:
                        filtered_train.append(j)
                cv_i_orig = deepcopy(cv[i])
                cv[i] = (filtered_train, test)
                if learn_options["num_genes_remove_train"] == 0:
                    if np.any(cv_i_orig[0] != cv[i][0]):
                        raise AssertionError()
                    if np.any(cv_i_orig[1] != cv[i][1]):
                        raise AssertionError()
                print(
                    f"# train/train after/before is {len(cv[i][0])}, {len(cv_i_orig[0])}"
                )
                print(
                    f"# test/test after/before is {len(cv[i][1])}, {len(cv_i_orig[1])}"
                )
    else:
        raise Exception(f"invalid cv options given: {learn_options['cv']}")

    cv = [c for c in cv]  # make list from generator, so can subset for TEST case
    if TEST:
        ind_to_use = [0]  # [0,1]
        cv = [cv[i] for i in ind_to_use]
        fold_labels = [fold_labels[i] for i in ind_to_use]

    truth = dict(
        [
            (t, dict([(m, np.array([])) for m in ["raw", "ranks", "thrs"]]))
            for t in fold_labels
        ]
    )
    predictions = dict([(t, np.array([])) for t in fold_labels])

    m = {}
    metrics = []

    # do the cross-validation
    num_proc = learn_options["num_proc"]
    X = inputs
    if num_proc > 1:
        num_proc = np.min([num_proc, len(cv)])
        print(f"using multiprocessing with {num_proc} procs -- one for each fold")
        jobs = []
        pool = Pool(processes=num_proc)
        for i, fold in enumerate(cv):
            train, test = fold
            print(
                f"working on fold {i+1} of {len(cv)}, with {len(train)} train and {len(test)} test"
            )
            if learn_options["method"] == "GPy":
                job = pool.apply_async(
                    GP.gp_on_fold,
                    args=(feature_sets, train, test, y, y_all, learn_options),
                )
            elif learn_options["method"] == "linreg":
                job = pool.apply_async(
                    regression.linreg_on_fold,
                    args=(train, test, y, y_all, X, learn_options, fold),
                )
            elif learn_options["method"] == "logregL1":
                job = pool.apply_async(
                    regression.logreg_on_fold,
                    args=(train, test, y, y_all, X, learn_options),
                )
            elif learn_options["method"] == "AdaBoostRegressor":
                job = pool.apply_async(
                    ensembles.adaboost_on_fold,
                    args=(train, test, y, y_all, X, learn_options, False),
                )
            elif learn_options["method"] == "AdaBoostClassifier":
                job = pool.apply_async(
                    ensembles.adaboost_on_fold,
                    args=(train, test, y, y_all, X, learn_options, True),
                )
            elif learn_options["method"] == "DecisionTreeRegressor":
                job = pool.apply_async(
                    ensembles.decisiontree_on_fold, args=(train, test, y, X)
                )
            elif learn_options["method"] == "RandomForestRegressor":
                job = pool.apply_async(
                    ensembles.randomforest_on_fold, args=(train, test, y, X)
                )
            elif learn_options["method"] == "ARDRegression":
                job = pool.apply_async(
                    regression.ARDRegression_on_fold, args=(train, test, y, X)
                )
            elif learn_options["method"] == "random":
                job = pool.apply_async(baselines.random_on_fold, args=(test))
            elif learn_options["method"] == "mean":
                job = pool.apply_async(baselines.mean_on_fold, args=(train, test, y))
            elif learn_options["method"] == "SVC":
                job = pool.apply_async(
                    baselines.SVC_on_fold, args=(train, test, y_all, X, learn_options)
                )
            elif learn_options["method"] == "DNN":
                job = pool.apply_async(
                    DNN.DNN_on_fold, args=(train, test, y_all, X, learn_options)
                )
            elif learn_options["method"] == "lasso_ensemble":
                job = pool.apply_async(
                    ensembles.LASSOs_ensemble_on_fold,
                    args=(feature_sets, train, test, y, y_all, X, learn_options),
                )
            elif learn_options["method"] == "doench":
                job = pool.apply_async(
                    baselines.doench_on_fold,
                    args=(train, test, y, y_all, X, learn_options),
                )
            elif learn_options["method"] == "sgrna_from_doench":
                job = pool.apply_async(
                    baselines.sgrna_from_doench_on_fold, args=(feature_sets, test, X)
                )
            elif learn_options["method"] == "xu_et_al":
                job = pool.apply_async(
                    baselines.xu_et_al_on_fold, args=(test, X, learn_options)
                )
            else:
                raise Exception(f"did not find method={learn_options['method']}")
            jobs.append(job)
        pool.close()
        pool.join()
        print(f"finished fold {i + 1}")
        for i, fold in enumerate(cv):  # i in range(0,len(jobs)):
            y_pred, m[i] = jobs[i].get()
            train, test = fold

            if learn_options["training_metric"] == "AUC":
                extract_fpr_tpr_for_fold(
                    aucs=metrics,
                    y_binary=y_all[learn_options["ground_truth_label"]].values,
                    test=test,
                    y_pred=y_pred,
                )
            elif learn_options["training_metric"] == "NDCG":
                extract_NDCG_for_fold(
                    metrics=metrics,
                    y_ground_truth=y_all[learn_options["ground_truth_label"]].values,
                    test=test,
                    y_pred=y_pred,
                    learn_options=learn_options,
                )
            elif learn_options["training_metric"] == "spearmanr":
                extract_spearman_for_fold(
                    metrics=metrics,
                    y_ground_truth=y_all[learn_options["ground_truth_label"]].values,
                    test=test,
                    y_pred=y_pred,
                )
            else:
                raise Exception(
                    f"invalid 'training_metric' in learn_options: {learn_options['training_metric']}"
                )

            truth, predictions = fill_in_truth_and_predictions(
                truth, predictions, fold_labels[i], y_all, y_pred, learn_options, test
            )

        pool.terminate()

    else:
        # non parallel version
        for i, fold in enumerate(cv):
            train, test = fold
            if learn_options["method"] == "GPy":
                y_pred, m[i] = GP.gp_on_fold(
                    feature_sets, train, test, y, y_all, learn_options
                )
            elif learn_options["method"] == "linreg":
                y_pred, m[i] = regression.linreg_on_fold(
                    train, test, y, y_all, X, learn_options
                )
            elif learn_options["method"] == "logregL1":
                y_pred, m[i] = regression.logreg_on_fold(
                    train, test, y, y_all, X, learn_options
                )
            elif learn_options["method"] == "AdaBoostRegressor":
                y_pred, m[i] = ensembles.adaboost_on_fold(
                    train, test, y, y_all, X, learn_options, classification=False
                )
            elif learn_options["method"] == "AdaBoostClassifier":
                y_pred, m[i] = ensembles.adaboost_on_fold(
                    train, test, y, y_all, X, learn_options, classification=True
                )
            elif learn_options["method"] == "DecisionTreeRegressor":
                y_pred, m[i] = ensembles.decisiontree_on_fold(train, test, y, X)
            elif learn_options["method"] == "RandomForestRegressor":
                y_pred, m[i] = ensembles.randomforest_on_fold(train, test, y, X)
            elif learn_options["method"] == "ARDRegression":
                y_pred, m[i] = regression.ARDRegression_on_fold(train, test, y, X)
            elif learn_options["method"] == "random":
                y_pred, m[i] = baselines.random_on_fold(test)
            elif learn_options["method"] == "mean":
                y_pred, m[i] = baselines.mean_on_fold(train, test, y)
            elif learn_options["method"] == "SVC":
                y_pred, m[i] = baselines.SVC_on_fold(
                    train, test, y_all, X, learn_options
                )
            elif learn_options["method"] == "DNN":
                y_pred, m[i] = DNN.DNN_on_fold(train, test, y_all, X, learn_options)
            elif learn_options["method"] == "lasso_ensemble":
                y_pred, m[i] = ensembles.LASSOs_ensemble_on_fold(
                    feature_sets, train, test, y, y_all, X, learn_options
                )
            elif learn_options["method"] == "doench":
                y_pred, m[i] = baselines.doench_on_fold(
                    train, test, y, y_all, X, learn_options
                )
            elif learn_options["method"] == "sgrna_from_doench":
                y_pred, m[i] = baselines.sgrna_from_doench_on_fold(
                    feature_sets, test, X
                )
            elif learn_options["method"] == "xu_et_al":
                y_pred, m[i] = baselines.xu_et_al_on_fold(test, X, learn_options)
            else:
                raise Exception(f"invalid method found: {learn_options['method']}")

            if learn_options["training_metric"] == "AUC":
                # fills in truth and predictions
                extract_fpr_tpr_for_fold(
                    aucs=metrics,
                    y_binary=y_all[learn_options["ground_truth_label"]].values,
                    test=test,
                    y_pred=y_pred,
                )
            elif learn_options["training_metric"] == "NDCG":
                extract_NDCG_for_fold(
                    metrics=metrics,
                    y_ground_truth=y_all[learn_options["ground_truth_label"]].values,
                    test=test,
                    y_pred=y_pred,
                    learn_options=learn_options,
                )
            elif learn_options["training_metric"] == "spearmanr":
                extract_spearman_for_fold(
                    metrics=metrics,
                    y_ground_truth=y_all[learn_options["ground_truth_label"]].values,
                    test=test,
                    y_pred=y_pred,
                )

            truth, predictions = fill_in_truth_and_predictions(
                truth, predictions, fold_labels[i], y_all, y_pred, learn_options, test
            )

            print(f"\t\tRMSE: {np.sqrt(((y_pred - y[test]) ** 2).mean())}")
            print(
                f"\t\tSpearman correlation: {np.nan_to_num(spearmanr_nonan(y[test], y_pred)[0])}"
            )
            print(f"\t\tfinished fold/gene {i + 1} of {len(fold_labels)}")

    cv_median_metric = [np.median(metrics)]
    gene_pred = [(truth, predictions)]
    print(
        f"\t\tmedian {learn_options['training_metric']} across gene folds: {cv_median_metric[-1]:.3f}"
    )

    t3 = time()
    print(f"\t\tElapsed time for cv is {(t3 - t2):.{2}} seconds")
    return metrics, gene_pred, fold_labels, m, dimsum, filename, feature_names
