"""
from https://gist.github.com/bwhite/3726239

Information Retrieval metrics

Useful Resources:
http://www.cs.utexas.edu/~mooney/ir-course/slides/Evaluation.ppt
http://www.nii.ac.jp/TechReports/05-014E.pdf
http://www.stanford.edu/class/cs276/handouts/EvaluationNew-handout-6-per.pdf
http://hal.archives-ouvertes.fr/docs/00/72/67/60/PDF/07-busa-fekete.pdf
Learning to Rank for Information Retrieval (Tie-Yan Liu)
"""
from time import time

import numpy as np
from scipy.stats.mstats import rankdata

from .elevation.metrics import spearman_weighted_swap_perm_test


def mean_reciprocal_rank(relevance_scores: list) -> np.ndarray:
    """Score is reciprocal of the rank of the first relevant item

    First element is 'rank 1'.  Relevance is binary (nonzero is relevant).

    Example from http://en.wikipedia.org/wiki/Mean_reciprocal_rank
    > rs = [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
    > mean_reciprocal_rank(rs)
    0.61111111111111105
    > rs = np.array([[0, 0, 0], [0, 1, 0], [1, 0, 0]])
    > mean_reciprocal_rank(rs)
    0.5
    > rs = [[0, 0, 0, 1], [1, 0, 0], [1, 0, 0]]
    > mean_reciprocal_rank(rs)
    0.75

    Args:
        relevance_scores: Iterator of relevance scores (list or numpy) in rank order
            (first element is the first item)

    Returns:
        Mean reciprocal rank
    """
    relevance_scores = (np.asarray(r).nonzero()[0] for r in relevance_scores)
    return np.mean([1.0 / (r[0] + 1) if r.size else 0.0 for r in relevance_scores])


def r_precision(relevance: list) -> np.ndarray:
    """Score is precision after all relevant documents have been retrieved

    Relevance is binary (nonzero is relevant).

    > r = [0, 0, 1]
    > r_precision(r)
    0.33333333333333331
    > r = [0, 1, 0]
    > r_precision(r)
    0.5
    > r = [1, 0, 0]
    > r_precision(r)
    1.0

    Args:
        r: Relevance scores (list or numpy) in rank order
            (first element is the first item)

    Returns:
        R Precision
    """
    relevance = np.asarray(relevance) != 0
    z = relevance.nonzero()[0]
    if not z.size:
        return 0.0
    return np.mean(relevance[: z[-1] + 1])


def precision_at_k(r, k):
    """Score is precision @ k

    Relevance is binary (nonzero is relevant).

    > r = [0, 0, 1]
    > precision_at_k(r, 1)
    0.0
    > precision_at_k(r, 2)
    0.0
    > precision_at_k(r, 3)
    0.33333333333333331
    > precision_at_k(r, 4)
    Traceback (most recent call last):
        File "<stdin>", line 1, in ?
    ValueError: Relevance score length < k


    Args:
        r: Relevance scores (list or numpy) in rank order
            (first element is the first item)

    Returns:
        Precision @ k

    Raises:
        ValueError: len(r) must be >= k
        :param k:
    """
    if k < 1:
        raise AssertionError()
    r = np.asarray(r)[:k] != 0
    if r.size != k:
        raise ValueError("Relevance score length < k")
    return np.mean(r)


def average_precision(r):
    """Score is average precision (area under PR curve)

    Relevance is binary (nonzero is relevant).

    > r = [1, 1, 0, 1, 0, 1, 0, 0, 0, 1]
    > delta_r = 1. / sum(r)
    > sum([sum(r[:x + 1]) / (x + 1.) * delta_r for x, y in enumerate(r) if y])
    0.7833333333333333
    > average_precision(r)
    0.78333333333333333

    Args:
        r: Relevance scores (list or numpy) in rank order
            (first element is the first item)

    Returns:
        Average precision
    """
    r = np.asarray(r) != 0
    out = [precision_at_k(r, k + 1) for k in range(r.size) if r[k]]
    if not out:
        return 0.0
    return np.mean(out)


def mean_average_precision(rs):
    """Score is mean average precision

    Relevance is binary (nonzero is relevant).

    >>> rs = [[1, 1, 0, 1, 0, 1, 0, 0, 0, 1]]
    >>> mean_average_precision(rs)
    0.78333333333333333
    >>> rs = [[1, 1, 0, 1, 0, 1, 0, 0, 0, 1], [0]]
    >>> mean_average_precision(rs)
    0.39166666666666666

    Args:
        rs: Iterator of relevance scores (list or numpy) in rank order
            (first element is the first item)

    Returns:
        Mean average precision
    """
    return np.mean([average_precision(r) for r in rs])


def dcg_at_k(r, k, method=0):
    """Score is discounted cumulative gain (dcg)

    Relevance is positive real values.  Can use binary
    as the previous methods.

    Example from
    http://www.stanford.edu/class/cs276/handouts/EvaluationNew-handout-6-per.pdf
    > r = [3, 2, 3, 0, 0, 1, 2, 2, 3, 0]
    > dcg_at_k(r, 1)
    3.0
    > dcg_at_k(r, 1, method=1)
    3.0
    > dcg_at_k(r, 2)
    5.0
    > dcg_at_k(r, 2, method=1)
    4.2618595071429155
    > dcg_at_k(r, 10)
    9.6051177391888114
    > dcg_at_k(r, 11)
    9.6051177391888114

    Args:
        r: Relevance scores (list or numpy) in rank order
            (first element is the first item)
        k: Number of results to consider
        method: If 0 then weights are [1.0, 1.0, 0.6309, 0.5, 0.4307, ...]
                If 1 then weights are [1.0, 0.6309, 0.5, 0.4307, ...]

    Returns:
        Discounted cumulative gain
    """
    r = np.asfarray(r)[:k]
    if r.size:
        if method == 0:
            return r[0] + np.sum(r[1:] / np.log2(np.arange(2, r.size + 1)))
        elif method == 1:
            return np.sum(r / np.log2(np.arange(2, r.size + 2)))
        else:
            raise ValueError("method must be 0 or 1.")
    return 0.0


def ndcg_at_k(r, k, method=0):
    """Score is normalized discounted cumulative gain (ndcg)

    Relevance is positive real values.  Can use binary
    as the previous methods.

    Example from
    http://www.stanford.edu/class/cs276/handouts/EvaluationNew-handout-6-per.pdf
    >>> r = [3, 2, 3, 0, 0, 1, 2, 2, 3, 0]
    >>> ndcg_at_k(r, 1)
    1.0
    >>> r = [2, 1, 2, 0]
    >>> ndcg_at_k(r, 4)
    0.9203032077642922
    >>> ndcg_at_k(r, 4, method=1)
    0.96519546960144276
    >>> ndcg_at_k([0], 1)
    0.0
    >>> ndcg_at_k([1], 2)
    1.0

    Args:
        r: Relevance scores (list or numpy) in rank order
            (first element is the first item)
        k: Number of results to consider
        method: If 0 then weights are [1.0, 1.0, 0.6309, 0.5, 0.4307, ...]
                If 1 then weights are [1.0, 0.6309, 0.5, 0.4307, ...]

    Returns:
        Normalized discounted cumulative gain
    """
    dcg_max = dcg_at_k(sorted(r, reverse=True), k, method)
    if not dcg_max:
        return 0.0
    return dcg_at_k(r, k, method) / dcg_max


# ------------------------------------------------------------------------------------
# custom stuff from us to avoid problem with ties


def ndcg_at_k_ties(
    labels: list,
    predictions: list,
    k: int,
    method: int = 0,
    normalize_from_below_too: bool = False,
    theta=None,
) -> float:
    """
    See 2008 McSherry et al on how to efficiently compute NDCG with ties
    labels are ground truth

    if k=None then k gets set to len(labels)

    labels and predictions get flattened here

    set normalize_from_below_too=False for conventional
    ndcg_at_k_ties, but note this will only
    ensure the max is 1, not that the min is zero.
    to get that added guarantee, set this argument to True
    """

    if isinstance(labels, list):
        labels = np.array(labels)
    if isinstance(predictions, list):
        predictions = np.array(predictions)

    if len(labels.shape) != 1 and np.min(labels.shape) != 1:
        raise AssertionError("should be 1D array or equivalent")
    if len(predictions.shape) != 1 and np.min(predictions.shape) != 1:
        raise AssertionError("should be 1D array or equivalent")

    labels = labels.flatten()
    predictions = predictions.flatten()

    if np.any(labels.shape != predictions.shape):
        raise AssertionError("labels and predictions should have the same shape")

    if k is None:
        k = len(labels)

    labels = labels.copy()

    dcg = dcg_at_k_ties(labels, predictions, k, method=method, theta=theta)

    dcg_max = dcg_at_k_ties(labels, labels, k, method, theta=theta)
    # NOTE: I have checked that dcg_at_k_ties and dcg_at_k match when there are no ties,
    # or ties in the labels

    if normalize_from_below_too:
        dcg_min = dcg_at_k_ties(
            np.sort(labels)[::-1], np.sort(predictions), k, method, theta=theta
        )
    else:
        dcg_min = 0
    numerator = dcg - dcg_min
    if numerator <= -1e-5:
        raise AssertionError()
    numerator = np.max((0, numerator))
    ndcg = numerator / (dcg_max - dcg_min)
    if not 1.0 >= ndcg >= 0.0:
        raise AssertionError(f"ndcg={ndcg} should be in [0,1]")
    if not dcg_max:
        ndcg = 0.0
    return ndcg


def dcg_helper(discount_factors, gain, k, labels, method, predictions):
    # step through, in current order (of decreasing predictions), accumulating tied gains
    # (which may be singletons)
    ii = 0
    dcg = 0.0
    while ii < k:
        current_pred = predictions[ii]
        current_gain = gain(labels[ii], method)
        # intializing the tied cumulative variables
        cum_tied_gain = current_gain
        cum_tied_disc = discount_factors[ii]
        num_ties = 1
        ii += 1
        # count number of ties in predictions
        while ii < len(predictions) and predictions[ii] == current_pred:  # while tied
            num_ties += 1.0
            cum_tied_gain += gain(labels[ii], method)
            if ii < k:
                cum_tied_disc += discount_factors[ii]
            ii += 1
        avg_gain = cum_tied_gain / num_ties
        dcg += avg_gain * cum_tied_disc
        if np.isnan(dcg):
            raise AssertionError("found nan dcg")
    return dcg


def dcg_at_k_ties(labels, predictions, k, method=0, theta=None):
    """
    See 2008 McSherry et al on how to efficiently compute NDCG (method=0 here) with ties
    (in the predictions)
    'labels' are what the "ground truth" judges assign
    'predictions' are the algorithm predictions corresponding to each label
    Also, http://en.wikipedia.org/wiki/Discounted_cumulative_gain for basic defns
    """
    if not isinstance(predictions, np.ndarray):
        raise AssertionError()
    if len(labels) != len(predictions):
        raise AssertionError("labels and predictions should be of same length")
    if k > len(labels):
        raise AssertionError("k should be <= len(labels)")

    # order both labels and preds so that they are in order of decreasing predictive score
    sorted_ind = np.argsort(predictions)[::-1]
    predictions = predictions[sorted_ind]
    labels = labels[sorted_ind]

    def gain(label, method):
        if method == 0:
            return label
        elif method == 1:
            return 2 ** label - 1.0
        elif method == 2 or method == 3 or method == 4:
            return label
        else:
            raise NotImplementedError()

    if method == 0:
        discount_factors = get_discount_factors(len(labels), discount="log2")
    elif method == 1:
        raise Exception("need to implement: log_2(i+1)")
    elif method == 2:
        discount_factors = get_discount_factors(len(labels), discount="linear")
    elif method == 3:
        discount_factors = get_discount_factors(len(labels), discount="combination")
    elif method == 4:
        if theta is None:
            raise AssertionError("need to specify theta or theta")
        discount_factors = get_discount_factors(
            len(labels), discount="1/rtheta", theta=theta
        )

    else:
        raise NotImplementedError()

    if len(discount_factors) != len(labels):
        raise AssertionError("discount factors has wrong length")

    dcg = dcg_helper(discount_factors, gain, k, labels, method, predictions)
    if np.isnan(dcg):
        raise AssertionError("found nan dcg")

    return dcg


def get_discount_factors(num_labels, discount="log2", theta=None):
    ii_range = np.arange(num_labels) + 1

    if discount == "log2":
        discount_factors = np.concatenate(
            (np.array([1.0]), 1.0 / np.log2(ii_range[1:]))
        )
    elif discount == "linear":
        discount_factors = -ii_range / float(num_labels) + 1.0
    elif discount == "combination":
        l2 = np.concatenate((np.array([1.0]), 1.0 / np.log2(ii_range[1:])))
        linear = -ii_range / float(num_labels) + 1.0
        discount_factors = np.max((l2, linear), axis=0)
    elif discount == "1/rtheta":
        discount_factors = 1.0 / (ii_range ** theta)
    else:
        raise NotImplementedError

    return discount_factors


def rank_data(r, rground):
    # we checked this heavily, and is correct, e.g. rground will go from largest rank to smallest
    r = rankdata(r)
    rground = rankdata(rground)
    if np.sum(r) != np.sum(rground):
        raise AssertionError("ranks should add up to the same")
    return r, rground


def dcg_alt(relevances, rank=20):
    relevances = np.asarray(relevances)[:rank]
    n_relevances = len(relevances)
    if n_relevances == 0:
        return 0.0
    discounts = np.log2(np.arange(n_relevances) + 2)
    return np.sum(relevances / discounts)


def ndcg_alt(relevances, rank=20):
    best_dcg = dcg_alt(sorted(relevances, reverse=True), rank)
    if best_dcg == 0:
        return 0.0
    return dcg_alt(relevances, rank) / best_dcg


def ndcg_at_k_swap_perm_test(
    preds1, preds2, true_labels, nperm, method, k, normalize_from_below_too, theta=None
):
    # pVal is the probability that we would observe as big an AUC diff as we
    # did if the ROC curves were drawn from the null hypothesis (which is that
    # one model does not perform better than the other)
    #
    # null hypothesis is that the prediction ranking are the same, so we exchange a random
    # number of them with each other.
    #
    # see ndcg_at_k_ties for all but the first four parameters
    #
    # balance_zeros = True means that when we swap a zero for a non-zero value, we will also do
    # a reverse swap
    #
    # this is a two-sided test, but since it is a symmetric null distribution, one should
    # be able to divide the p-value by 2 to get the one-sided version (but think this through
    # before using)

    if isinstance(preds1, list):
        preds1 = np.array(preds1)
    else:
        preds1 = preds1.flatten()

    if isinstance(preds2, list):
        preds2 = np.array(preds2)
    else:
        preds2 = preds2.flatten()

    if isinstance(true_labels, list):
        true_labels = np.array(true_labels)
    else:
        true_labels = true_labels.flatten()

    if len(preds1) != len(preds2):
        raise AssertionError("need same number of preditions from each model")
    if len(preds1) != len(true_labels):
        raise AssertionError("need same number of preditions in truth and predictions")
    N = len(preds1)

    # re-sort all by truth ordering so that when swap they are aligned
    sorted_ind = np.argsort(true_labels)[::-1]
    true_labels = true_labels[sorted_ind]
    preds1 = preds1[sorted_ind]
    preds2 = preds2[sorted_ind]

    ranks1 = rankdata(preds1)
    ranks2 = rankdata(preds2)

    ndcg1 = ndcg_at_k_ties(
        true_labels,
        ranks1,
        k=k,
        method=method,
        normalize_from_below_too=normalize_from_below_too,
        theta=theta,
    )
    ndcg2 = ndcg_at_k_ties(
        true_labels,
        ranks2,
        k=k,
        method=method,
        normalize_from_below_too=normalize_from_below_too,
        theta=theta,
    )

    real_ndcg_diff = np.abs(ndcg1 - ndcg2)
    perm_ndcg_diff = np.nan * np.zeros(nperm)

    if np.all(preds1 == preds2):
        pval = 1.0
    else:
        zero_ind = true_labels == 0
        if np.sum(zero_ind) >= len(zero_ind):
            raise AssertionError("balancing assumes there are more zeros than ones")

        for _ in range(nperm):
            pair_ind_to_swap = np.random.rand(N) < 0.5

            ranks1_perm = ranks1.copy()
            ranks1_perm[pair_ind_to_swap] = ranks2[pair_ind_to_swap]

            ranks2_perm = ranks2.copy()
            ranks2_perm[pair_ind_to_swap] = ranks1[pair_ind_to_swap]

            ndcg1_perm = ndcg_at_k_ties(
                true_labels,
                ranks1_perm,
                k=k,
                method=method,
                normalize_from_below_too=normalize_from_below_too,
                theta=theta,
            )
            ndcg2_perm = ndcg_at_k_ties(
                true_labels,
                ranks2_perm,
                k=k,
                method=method,
                normalize_from_below_too=normalize_from_below_too,
                theta=theta,
            )

            for thing in theta:
                tmp_diff = np.abs(ndcg1_perm[thing] - ndcg2_perm[thing])
                perm_ndcg_diff[thing][_] = tmp_diff

        num_stat_greater = np.max((((perm_ndcg_diff > real_ndcg_diff).sum() + 1), 1.0))
        pval = num_stat_greater / nperm

    return pval, real_ndcg_diff, perm_ndcg_diff, ndcg1, ndcg2


if __name__ == "__main__":
    simulated_data = True
    permute_real_data = True

    T = 1000

    nperm = 100

    weights = np.array([0.001])
    theta_range = weights  # just to make life easier

    # only for simulated data
    N = 100
    frac_zeros = 0

    k = None

    allp = np.nan * np.zeros((len(theta_range) + 1, T))

    if not simulated_data:
        # print(
        #     "loading up saved data..."
        # )  # two-fold CV data from CRISPR off-target GUIDE-SEQ
        # with open(r"\\nerds5\kevin\from_nicolo\gs.pickle", "rb") as f:
        #     predictions, truth_all = pickle.load(f)
        # print("done.")
        # N = len(truth_all[0])
        pass  # that gs.pickle file was not in the source repo

    for t in range(T):

        # totally simulated
        if simulated_data:
            truth = np.random.rand(N)
            zero_ind = np.random.rand(N) < frac_zeros
            truth[zero_ind] = 0
            pred1 = np.random.rand(N)
            pred2 = np.random.rand(N)
        # this all refers to stuff from that unavailable gs.pickle from above
        # else:
        #     fold = 0
        #     truth = truth_all[fold]
        #     pred1 = predictions["CFD"][fold]
        #     pred2 = predictions["product"][fold]

        #     if permute_real_data:
        #         truth = np.random.permutation(truth)

        t0 = time()
        for i, w in enumerate(weights):
            weights_array = truth.copy()
            weights_array += w

            pvaltmp, real_corr_diff, perm_corr_diff, corr1, corr2 = spearman_weighted_swap_perm_test(
                pred1, pred2, truth, nperm, weights_array
            )

            allp[i, t] = pvaltmp
            t1 = time()

    truth = np.array([3, 4, 2, 1, 0, 0, 0])
    pred1 = np.array([3, 4, 2, 1, 0, 0, 0])
    pred2 = np.array([2, 1, 3, 4, 5, 6, 7])

    truth3 = np.array([3, 4, 2, 1, 0, 0, 0])
    truth4 = np.zeros(7)
    truth4[0] = 1
    pred3 = np.array([2, 1, 3, 4, 5, 6, 7]) * 10
    pred4 = np.array([4, 3, 2, 1, 0, 0, 0])
    pred5 = np.array([4, 3, 1, 2, 0, 0, 0])

    nperm = 1000
    method = 4
    theta = 0.5
    normalize_from_below_too = True
    k = len(pred3)

    pval, real_ndcg_diff, perm_ndcg_diff, ndcg1, ndcg2 = ndcg_at_k_swap_perm_test(
        pred1, pred2, truth, nperm, method, k, normalize_from_below_too, theta=theta
    )
    print(f"ndcg1={ndcg1}, ndcg2={ndcg2}, ndcg_diff={real_ndcg_diff}, p={pval}")

    pval, real_ndcg_diff, perm_ndcg_diff, ndcg1, ndcg2 = ndcg_at_k_swap_perm_test(
        pred1, pred1, truth, nperm, method, k, normalize_from_below_too, theta=theta
    )
    print(f"ndcg1={ndcg1}, ndcg2={ndcg2}, ndcg_diff={real_ndcg_diff}, p={pval}")

    pval, real_ndcg_diff, perm_ndcg_diff, ndcg1, ndcg2 = ndcg_at_k_swap_perm_test(
        pred1, pred4, truth, nperm, method, k, normalize_from_below_too, theta=theta
    )
    print(f"ndcg1={ndcg1}, ndcg2={ndcg2}, ndcg_diff={real_ndcg_diff}, p={pval}")

    pval, real_ndcg_diff, perm_ndcg_diff, ndcg1, ndcg2 = ndcg_at_k_swap_perm_test(
        pred1, pred5, truth, nperm, method, k, normalize_from_below_too, theta=theta
    )
    print(f"ndcg1={ndcg1}, ndcg2={ndcg2}, ndcg_diff={real_ndcg_diff}, p={pval}")

    print(ndcg_at_k_ties(truth4, pred2, k, method=3, normalize_from_below_too=True))

    print(ndcg_alt(truth[np.argsort(pred2)[::-1]], 5))
    print(ndcg_at_k(truth[np.argsort(pred2)[::-1]], 5, method=1))
    print(ndcg_at_k(truth[np.argsort(pred2)[::-1]], 5, method=0))

    print(ndcg_at_k_ties(truth, pred2, 5, method=1))
    print(ndcg_at_k_ties(truth, pred2, 5, method=0))
