from numbers import Number

import numpy as np
from numpy.core.multiarray import ndarray
from scipy.stats import spearmanr
from sklearn import linear_model
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder

from .. import predict
from ..metrics import ndcg_at_k_ties


def ARDRegression_on_fold(train, test, y, X):
    """
    """
    clf = linear_model.ARDRegression()
    clf.fit(X[train], y[train][:, 0])
    y_pred = clf.predict(X[test])[:, None]
    return y_pred, clf


def train_linreg_model(alpha, l1r, learn_options, fold, X, y, y_all):
    """
    fold is something like train_inner (boolean array specifying what is in the fold)
    """
    if learn_options["penalty"] == "L2":
        clf = linear_model.Ridge(
            alpha=alpha,
            fit_intercept=learn_options["fit_intercept"],
            normalize=learn_options["normalize_features"],
            copy_X=True,
            max_iter=None,
            tol=0.001,
            solver="auto",
        )
        weights = get_weights(learn_options, fold, y, y_all)
        clf.fit(X[fold], y[fold], sample_weight=weights)
    elif learn_options["penalty"] == "EN" or learn_options["penalty"] == "L1":
        if learn_options["loss"] == "squared":
            clf = linear_model.ElasticNet(
                alpha=alpha,
                l1_ratio=l1r,
                fit_intercept=learn_options["fit_intercept"],
                normalize=learn_options["normalize_features"],
                # max_iter=30000,
                warm_start=True,
            )
        elif learn_options["loss"] == "huber":
            clf = linear_model.SGDRegressor(
                "huber",
                alpha=alpha,
                l1_ratio=l1r,
                fit_intercept=learn_options["fit_intercept"],
                n_iter=10,
                penalty="elasticnet",
                shuffle=True,
                max_iter=30000,
            )
        # print(f"Fitting X:{X} and y:{y}")  # TODO: remove after DEBUG
        clf.fit(X[fold], y[fold])
    elif learn_options["penalty"] is None:
        clf = linear_model.LinearRegression(
            fit_intercept=learn_options["fit_intercept"],
            normalize=learn_options["normalize_features"],
            copy_X=True,
        )
        # print(f"Fitting X:{X} and y:{y}")  # TODO: remove after DEBUG
        clf.fit(X[fold], y[fold])
    return clf


def logreg_on_fold(train, test, y, y_all, X, learn_options):
    """
    (L1/L2 penalized) logistic reggresion using scikitlearn
    """

    if len(np.unique(y)) > 2:
        raise AssertionError("if using logreg need binary targets")
    if learn_options["weighted"] is not None:
        raise AssertionError("cannot do weighted Log reg")
    if learn_options["feature_select"] is True:
        raise AssertionError(
            "cannot do feature selection yet in logistic regression--see "
            "linreg_on_fold to implement"
        )

    cv, n_folds = set_up_inner_folds(learn_options, y_all.iloc[train])

    if learn_options["penalty"] != "L1" and learn_options["penalty"] != "L2":
        raise AssertionError("can only use L1 or L2 with logistic regression")

    tol = 0.00001

    performance: ndarray = np.zeros((len(learn_options["alpha"]), 1))
    for train_inner, test_inner in cv:
        for i, alpha in enumerate(learn_options["alpha"]):
            clf: LogisticRegression = linear_model.LogisticRegression(
                penalty=learn_options["penalty"].lower(),
                dual=False,
                fit_intercept=learn_options["fit_intercept"],
                class_weight=(
                    learn_options["class_weight"]
                    if "class_weight" in learn_options
                    else None
                ),
                tol=tol,
                C=1.0 / alpha,
            )

            clf.fit(X[train][train_inner], y[train][train_inner].flatten())
            tmp_pred = clf.predict_proba(X[train][test_inner])[:, 1]

            if learn_options["training_metric"] == "AUC":
                fpr, tpr, _ = roc_curve(
                    y_all[learn_options["ground_truth_label"]][train][test_inner],
                    tmp_pred,
                )
                if np.any(np.isnan(fpr)):
                    raise AssertionError("found nan fpr")
                if np.any(np.isnan(tpr)):
                    raise AssertionError("found nan tpr")

                tmp_auc = auc(fpr, tpr)
                performance[i] += tmp_auc
            else:
                raise Exception("can only use AUC metric for cv with classification")

    performance /= n_folds

    max_score_ind = np.where(performance == np.nanmax(performance))
    if max_score_ind == len(performance):
        raise AssertionError("enlarge alpha range as hitting max boundary")

    # in the unlikely event of tied scores, take the first one.
    if len(max_score_ind[0]) > 1:
        max_score_ind = [max_score_ind[0][0], max_score_ind[1][0]]

    best_alpha = learn_options["alpha"][max_score_ind[0]]

    best_alpha = best_alpha[0]
    if not isinstance(best_alpha, Number):
        raise Exception(f"best_alpha must be a number but is {type(best_alpha)}")

    print("\tbest alpha is {best_alpha} from range={learn_options['alpha'][[0, -1]]}")
    max_perf = np.nanmax(performance)

    if max_perf < 0.0:
        raise Exception("performance is negative")

    print(f"\t\tbest performance is {np.nanmax(performance)}")

    clf = linear_model.LogisticRegression(
        penalty=learn_options["penalty"],
        dual=False,
        fit_intercept=learn_options["fit_intercept"],
        class_weight=learn_options["class_weight"],
        tol=tol,
        C=1.0 / best_alpha,
    )
    clf.fit(X[train], y[train].flatten())

    y_pred = clf.predict_proba(X[test])[:, 1]
    y_pred = y_pred[:, None]

    return y_pred, clf


def linreg_on_fold(train, test, y, y_all, X, learn_options):
    """
    linreg using scikitlearn, using more standard regression models with penalization requiring
    nested-cross-validation
    """

    if learn_options["weighted"] is not None and (
        learn_options["penalty"] != "L2" or learn_options["method"] != "linreg"
    ):
        raise NotImplementedError(
            "weighted prediction not implemented for any methods by L2 at the moment"
        )

    if "fit_intercept" not in learn_options:
        learn_options["fit_intercept"] = True
    if "normalize_features" not in learn_options:
        learn_options["normalize_features"] = True

    cv, n_folds = set_up_inner_folds(learn_options, y_all.iloc[train])

    if learn_options["penalty"] == "L1" or learn_options["penalty"] is None:
        l1_ratio = [1.0]
    elif learn_options["penalty"] == "L2":
        l1_ratio = [0.0]
    elif learn_options["penalty"] == "EN":  # elastic net
        l1_ratio = np.linspace(0.0, 1.0, 20)

    performance = np.zeros((len(learn_options["alpha"]), len(l1_ratio)))
    degenerate_pred = np.zeros((len(learn_options["alpha"])))
    num_tests = len(l1_ratio) * len(learn_options["alpha"]) * len(cv)
    test_position = 0
    for train_inner, test_inner in cv:
        # print(f"testing {train_inner} and {test_inner}")
        for i, alpha in enumerate(learn_options["alpha"]):
            for j, l1r in enumerate(l1_ratio):
                print(f"Test {test_position} of {num_tests} for fold {fold_number}")
                test_position += 1
                clf = train_linreg_model(
                    alpha,
                    l1r,
                    learn_options,
                    train_inner,
                    X[train],
                    y[train],
                    y_all.iloc[train],
                )
                if learn_options["feature_select"]:
                    clf, tmp_pred = feature_select(
                        clf, learn_options, test_inner, train_inner, X[train], y[train]
                    )
                else:
                    tmp_pred = clf.predict(X[train][test_inner])

                if learn_options["training_metric"] == "AUC":
                    fpr, tpr, _ = roc_curve(
                        y_all[learn_options["ground_truth_label"]][train][test_inner],
                        tmp_pred,
                    )
                    if np.any(np.isnan(fpr)):
                        raise AssertionError("found nan fpr")
                    if np.any(np.isnan(tpr)):
                        raise AssertionError("found nan tpr")
                    tmp_auc = auc(fpr, tpr)
                    performance[i, j] += tmp_auc

                elif learn_options["training_metric"] == "spearmanr":
                    spearman = np.nan_to_num(
                        spearmanr(
                            y_all[learn_options["ground_truth_label"]][train][
                                test_inner
                            ],
                            tmp_pred.flatten(),
                        )[0]
                    )
                    performance[i, j] += spearman

                elif learn_options["training_metric"] == "score":
                    performance[i, j] += clf.score(
                        X[test_inner],
                        y_all[learn_options["ground_truth_label"]][train][test_inner],
                    )

                elif learn_options["training_metric"] == "NDCG":
                    if "thresh" in learn_options["ground_truth_label"]:
                        raise AssertionError(
                            "for NDCG must not use thresholded ranks, but pure " "ranks"
                        )

                    tmp_truth = (
                        y_all[learn_options["ground_truth_label"]]
                        .values[train][test_inner]
                        .flatten()
                    )
                    tmp_perf = ndcg_at_k_ties(
                        tmp_truth, tmp_pred.flatten(), learn_options["NDGC_k"]
                    )
                    performance[i, j] += tmp_perf

                    degenerate_pred_tmp = len(np.unique(tmp_pred)) < len(tmp_pred) / 2.0
                    degenerate_pred[i] += degenerate_pred_tmp
                    # tmp_pred_r, tmp_truth_r = rank_data(tmp_pred, tmp_truth)

    performance /= n_folds

    max_score_ind = np.where(performance == np.nanmax(performance))
    if max_score_ind == len(performance):
        raise AssertionError("enlarge alpha range as hitting max boundary")

    # in the unlikely event of tied scores, take the first one.
    # we take the first one regardless because a change in NumPy made it an
    # error to use a single-value ndarray as an array index:
    # https://stackoverflow.com/questions/42128830/typeerror-only-integer-scalar-arrays-can-be-converted-to-a-scalar-index/42444003#42444003
    max_score_ind = [max_score_ind[0][0], max_score_ind[1][0]]

    best_alpha, best_l1r = (
        learn_options["alpha"][max_score_ind[0]],
        l1_ratio[max_score_ind[1]],
    )

    # try:
    #     print(f"\tbest alpha is {best_alpha} from range={learn_options['alpha'][[0, -1]]}")
    # except:
    #     raise Exception("Uh... something went wrong with figuring out 'best_alpha' at line 227 in regression.py.  Fix it.")

    if learn_options["penalty"] == "EN":
        print(f"\t\tbest l1_ratio is {best_l1r} from range={l1_ratio[[0, -1]]}")
    max_perf = np.nanmax(performance)

    if max_perf < 0.0:
        raise Exception("performance is negative")

    print(f"\t\tbest performance is {max_perf}")

    clf = train_linreg_model(best_alpha, l1r, learn_options, train, X, y, y_all)
    if learn_options["feature_select"]:
        # raise Exception("untested in a long time, should double check")
        clf, y_pred = feature_select(clf, learn_options, test, train, X, y)
    else:
        y_pred = clf.predict(X[test])

    if learn_options["penalty"] != "L2" and learn_options["penalty"] is not None:
        y_pred = y_pred[:, None]

    return y_pred, clf


def feature_select(clf, learn_options, test_inner, train_inner, X, y):
    if learn_options["weighted"] is not None:
        raise AssertionError(
            "cannot currently do feature selection with weighted regression"
        )
    if learn_options["loss"] == "huber":
        raise AssertionError("won't use huber loss function with feature selection")
    non_zero_coeff = clf.coef_ != 0.0
    if non_zero_coeff.sum() > 0:
        clf = linear_model.LinearRegression()
        clf.fit(X[train_inner][:, non_zero_coeff.flatten()], y[train_inner])
        tmp_pred = clf.predict(X[test_inner][:, non_zero_coeff.flatten()])
    else:
        tmp_pred = np.ones_like(test_inner)
    return clf, tmp_pred


def get_weights(learn_options, fold, y, y_all):
    """
    fold is an object like train_inner which is boolean for which indexes are in the fold
    """
    weights = None
    if learn_options["weighted"] == "variance":
        weights = 1.0 / y_all["variance"].values[fold]
    elif learn_options["weighted"] == "ndcg":
        # DCG: r[0] + np.sum(r[1:] / np.log2(np.arange(2, r.size + 1)))
        N = len(fold)
        r = np.ones(N)
        discount = np.concatenate(
            (np.array([r[0]]), r[1:] / np.log2(np.arange(2, r.size + 1)))
        )[::1]
        ind = np.argsort(y[fold], axis=0).flatten()
        weights = np.ones(len(ind))
        weights[ind] = discount
    elif learn_options["weighted"] == "rank":
        N = len(y[fold])
        inverse_ranks = (np.arange(N) + 1.0)[::-1]
        ind = np.argsort(y[fold], axis=0).flatten()
        weights = np.ones(len(ind))
        weights[ind] = inverse_ranks
    elif learn_options["weighted"] == "score":
        N = len(y[fold])
        score = y[fold] + np.abs(np.min(y[fold]))
        ind = np.argsort(y[fold], axis=0).flatten()
        weights = np.ones(len(ind))
        weights[ind] = score
    elif learn_options["weighted"] == "random":
        N = len(y[fold])
        weights = np.random.rand(N)
    elif learn_options["weighted"] is not None:
        raise Exception(f"invalid weighted type, {learn_options['weighted']}")
    # plt.plot(weights, y[train_inner],'.')
    return weights


def set_up_inner_folds(learn_options: dict, y):
    label_encoder = LabelEncoder()
    label_encoder.fit(y["Target gene"].values)
    gene_classes = label_encoder.transform(y["Target gene"].values)
    n_genes = len(np.unique(gene_classes))
    cv = []
    if (
        (
            "ignore_gene_level_for_inner_loop" in learn_options
            and learn_options["ignore_gene_level_for_inner_loop"]
        )
        or learn_options["cv"] == "stratified"
        or n_genes == 1
    ):
        if "n_folds" not in learn_options:
            n_splits = len(np.unique(gene_classes))
        else:
            n_splits = learn_options["n_folds"]
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True)
        cv = skf.split(np.zeros(len(gene_classes), dtype=np.bool), gene_classes)
    elif learn_options["cv"] == "gene":
        gene_list = np.unique(y["Target gene"].values)
        for gene in gene_list:
            cv.append(predict.get_train_test(gene, y))
    n_folds = len(cv)
    return cv, n_folds
