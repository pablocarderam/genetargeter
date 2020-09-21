import Bio.Seq as Seq
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.stats.mstats import rankdata
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model.coordinate_descent import ElasticNet

from .load_data import combine_organisms


def get_thirty_one_mer_data():
    """
    Load up our processed data file for all of V1 and V2, make a 31mer so that
    we can use the SSC trained model to compare to
    Assumes we call this from the analysis subdirectory
    """
    myfile = r"..\data\FC_plus_RES_withPredictions.csv"
    # was originally FC_RES_5304.csv but that isn't present
    newfile = r"..\data\FC_plus_RES_withPredictions_w_31mer.csv"
    data = pd.read_csv(myfile)
    thirty_one_mer = []
    for i in range(data.shape[0]):
        thirty_one_mer.append(
            convert_to_thirty_one(
                data.iloc[i]["30mer"], data.iloc[i]["Target"], data.iloc[i]["Strand"]
            )
        )
    data["31mer"] = thirty_one_mer
    data.to_csv(newfile)


def convert_to_thirty_one(guide_seq, gene, strand):
    """
    Given a guide sequence, a gene name, and strand (e.g. "sense"),
    return a 31mer string which is our 30mer,
    plus one more at the end.
    """
    guide_seq = Seq.Seq(guide_seq)
    gene_seq = Seq.Seq(get_gene_sequence(gene)).reverse_complement()
    if strand == "sense":
        guide_seq = guide_seq.reverse_complement()
    ind = gene_seq.find(guide_seq)
    if ind == -1:
        print(
            f"returning sequence+'A', could not find guide {guide_seq} in gene {gene}"
        )
        return gene_seq + "A"
    if gene_seq[ind : (ind + len(guide_seq))] != guide_seq:
        raise AssertionError("match not right")
    new_mer = gene_seq[(ind - 1) : (ind + len(guide_seq))]
    # this actually tacks on an extra one at the end for some reason
    if strand == "sense":
        new_mer = new_mer.reverse_complement()
    return str(new_mer)


def concatenate_feature_sets(feature_sets, keys=None):
    """
    Given a dictionary of sets of features, each in a pd.DataFrame,
    concatenate them together to form one big np.array, and get the dimension
    of each set
    Returns: inputs, dim
    """
    if feature_sets == {}:
        raise AssertionError("no feature sets present")
    if keys is None:
        keys = list(feature_sets.keys())

    F = feature_sets[keys[0]].shape[0]
    for assemblage in feature_sets:
        F2 = feature_sets[assemblage].shape[0]
        if F != F2:
            raise AssertionError(
                f"not same # individuals for features {keys[0]} and {assemblage}"
            )

    N = feature_sets[keys[0]].shape[0]
    inputs = np.zeros((N, 0))
    feature_names = []
    dim = {}
    dimsum = 0
    for assemblage in keys:
        inputs_set = feature_sets[assemblage].values
        dim[assemblage] = inputs_set.shape[1]
        dimsum = dimsum + dim[assemblage]
        inputs = np.hstack((inputs, inputs_set))
        feature_names.extend(feature_sets[assemblage].columns.tolist())

    return inputs, dim, dimsum, feature_names


def spearmanr_nonan(x, y):
    """
    same as scipy.stats.spearmanr, but if all values are equal, returns 0 instead of nan
    (Output: rho, pval)
    """
    r, p = spearmanr(x, y)
    r = np.nan_to_num(r)
    p = np.nan_to_num(p)
    return r, p


def impute_gene_position(gene_position):
    """
    Some amino acid cut position and percent peptide are blank because of stop codons, but
    we still want a number for these, so just set them to 101 as a proxy
    """

    gene_position["Percent Peptide"] = gene_position["Percent Peptide"].fillna(101.00)

    if "Amino Acid Cut position" in gene_position.columns:
        gene_position["Amino Acid Cut position"] = gene_position[
            "Amino Acid Cut position"
        ].fillna(gene_position["Amino Acid Cut position"].mean())

    return gene_position


def get_gene_sequence(gene_name):
    try:
        gene_file = f"../../gene_sequences/{gene_name}_sequence.txt"
        with open(gene_file, "rb") as f:
            seq = f.read()
            seq = seq.replace("\r\n", "")
    except ValueError:
        print(
            f"could not find gene sequence file {gene_file}, "
            f"please see examples and generate one for your gene "
            f"as needed, with this filename"
        )

    return seq


def get_ranks(y, thresh=0.8, prefix="", flip=False):
    """
    y should be a DataFrame with one column
    thresh is the threshold at which to call it a knock-down or not
    col_name = 'score' is only for V2 data
    flip should be FALSE for both V1 and V2!
    """

    if prefix is not None:
        prefix = prefix + "_"

    # y_rank = y.apply(ranktrafo)
    y_rank = y.apply(rankdata)
    y_rank /= y_rank.max()

    if flip:
        y_rank = (
            1.0 - y_rank
        )  # before this line, 1-labels where associated with low ranks, this flips it around
        # (hence the y_rank > thresh below)
        # we should NOT flip (V2), see README.txt in ./data

    y_rank.columns = [prefix + "rank"]
    y_threshold = (y_rank > thresh) * 1

    y_threshold.columns = [prefix + "threshold"]

    # JL: undo the log2 transform (not sure this matters?)
    y_rank_raw = (2 ** y).apply(rankdata)
    y_rank_raw /= y_rank_raw.max()
    if flip:
        y_rank_raw = 1.0 - y_rank_raw
    y_rank_raw.columns = [prefix + "rank raw"]
    if np.any(np.isnan(y_rank)):
        raise AssertionError("found NaN in ranks")

    y_quantized = y_threshold.copy()
    y_quantized.columns = [prefix + "quantized"]

    return y_rank, y_rank_raw, y_threshold, y_quantized


def get_data(data, y_names, organism="human", target_gene=None):
    """
    this is called once for each gene (aggregating across cell types)
    y_names are cell types
    e.g. call: X_CD13, Y_CD13 = get_data(cd13, y_names=['NB4 CD13', 'TF1 CD13'])
    """
    outputs = pd.DataFrame()
    # generate ranks for each cell type before aggregating to match what is in Doench et al
    thresh = 0.8
    for y_name in y_names:  # for each cell type
        y = pd.DataFrame(data[y_name])
        # these thresholds/quantils are not used:
        y_rank, y_rank_raw, y_threshold, _ = get_ranks(y, thresh=thresh, flip=False)
        y_rank.columns = [y_name + " rank"]
        y_rank_raw.columns = [y_name + " rank raw"]
        y_threshold.columns = [y_name + " threshold"]

        outputs = pd.concat([outputs, y, y_rank, y_threshold, y_rank_raw], axis=1)

    # aggregated rank across cell types
    average_activity = pd.DataFrame(outputs[[y_name for y_name in y_names]].mean(1))
    average_activity.columns = ["average activity"]

    average_rank_from_avg_activity = get_ranks(
        average_activity, thresh=thresh, flip=False
    )[0]
    average_rank_from_avg_activity.columns = ["average_rank_from_avg_activity"]
    average_threshold_from_avg_activity = (average_rank_from_avg_activity > thresh) * 1
    average_threshold_from_avg_activity.columns = [
        "average_threshold_from_avg_activity"
    ]

    average_rank = pd.DataFrame(
        outputs[[y_name + " rank" for y_name in y_names]].mean(1)
    )
    average_rank.columns = ["average rank"]
    # higher ranks are better (when flip=False as it should be)
    average_threshold = (average_rank > thresh) * 1
    average_threshold.columns = ["average threshold"]

    # undo the log2 trafo on the reads per million, apply rank trafo right away
    average_rank_raw = pd.DataFrame(
        outputs[[y_name + " rank raw" for y_name in y_names]].mean(1)
    )
    average_rank_raw.columns = ["average rank raw"]
    outputs = pd.concat(
        [
            outputs,
            average_rank,
            average_threshold,
            average_activity,
            average_rank_raw,
            average_rank_from_avg_activity,
            average_threshold_from_avg_activity,
        ],
        axis=1,
    )

    # import pdb; pdb.set_trace()

    # sequence-specific computations
    # features = featurize_data(data)
    # strip out featurization to later
    features = pd.DataFrame(data["30mer"])

    if organism == "human":
        target_gene = y_names[0].split(" ")[1]

    outputs["Target gene"] = target_gene
    outputs["Organism"] = organism

    features["Target gene"] = target_gene
    features["Organism"] = organism
    features["Strand"] = pd.DataFrame(data["Strand"])

    return features, outputs


def extract_feature_from_model(method, results, split):
    model_type = results[method][3][split]
    if isinstance(model_type, ElasticNet):
        tmp_imp = results[method][3][split].coef_[:, None]
    elif isinstance(model_type, GradientBoostingRegressor):
        tmp_imp = results[method][3][split].feature_importances_[:, None]
    else:
        raise Exception(f"need to add model {model_type} to feature extraction")
    return tmp_imp


def extract_feature_from_model_sum(method, results, split, indexes):
    model_type = results[method][3][split]
    if isinstance(model_type, ElasticNet):
        tmp_imp = np.sum(results[method][3][split].coef_[indexes])
    elif isinstance(model_type, GradientBoostingRegressor):
        tmp_imp = np.sum(results[method][3][split].feature_importances_[indexes])
    else:
        raise Exception(f"need to add model {model_type} to feature extraction")
    return tmp_imp


def feature_importances(results):
    for method in results:
        feature_names = results[method][6]

        seen = set()
        uniq = []
        for ft in feature_names:
            if ft not in seen:
                uniq.append(ft)
            else:
                seen.add(ft)
        if seen:
            raise Exception(f"feature name appears more than once: {seen}")

        pd_order1, pi_order1, pd_order2, pi_order2, nggx = [], [], [], [], []
        for i, s in enumerate(feature_names):
            if "False" in s:
                continue
            elif "_" in s:
                nucl, _ = s.split("_")
                if len(nucl) == 1:
                    pd_order1.append(i)
                elif len(nucl) == 2:
                    pd_order2.append(i)
            elif "NGGX_pd.Order2" in s:
                nggx.append(i)
            else:
                nucl = s
                if len(nucl) == 1:
                    pi_order1.append(i)
                elif len(nucl) == 2:
                    pi_order2.append(i)

        grouped_feat = {
            "pd_order2": pd_order2,
            "pi_order2": pi_order2,
            "pd_order1": pd_order1,
            "pi_order1": pi_order1,
            "NGGX_pd.Order2": nggx,
        }

        grouped_feat_ind = [grouped_feat[a] for a in grouped_feat]
        remaining_features_ind = set.difference(
            set(range(len(feature_names))), set(grouped_feat_ind)
        )

        for i in remaining_features_ind:
            grouped_feat[feature_names[i]] = [i]

        feature_importances_grouped = {}
        for k in grouped_feat:
            if not grouped_feat[k]:
                continue
            else:
                for split in results[method][3]:
                    split_feat_importance = extract_feature_from_model_sum(
                        method, results, split, grouped_feat[k]
                    )
                    if k not in feature_importances_grouped:
                        feature_importances_grouped[k] = [split_feat_importance]
                    else:
                        feature_importances_grouped[k].append(split_feat_importance)

        all_split_importances = None
        for split in results[method][3]:

            split_feat_importance = extract_feature_from_model(method, results, split)

            if all_split_importances is None:
                all_split_importances = split_feat_importance.copy()
            else:
                all_split_importances = np.append(
                    all_split_importances, split_feat_importance, axis=1
                )

        avg_importance = np.mean(all_split_importances, axis=1)[:, None]
        std_importance = np.std(all_split_importances, axis=1)[:, None]
        imp_array = np.concatenate(
            (np.array(feature_names)[:, None], avg_importance, std_importance), axis=1
        )

        df = pd.DataFrame(
            data=imp_array,
            columns=["Feature name", "Mean feature importance", "Std. Dev."],
        )
        df = df.convert_objects(convert_numeric=True)

        feature_dictionary = {
            "pd_order2": "position dep. order 2 ",
            "pd_order1": "position dep. order 1 ",
            "pi_order1": "position ind. order 1 ",
            "pi_order2": "position ind. order 2 ",
            "5mer_end_False": "Tm (5mer end)",
            "5mer_start_False": "Tm (5mer start)",
            "Amino Acid Cut position": "amino acid cut position ",
            "8mer_middle_False": "Tm (8mer middle)",
            "NGGX_pd.Order2": "NGGN interaction ",
            "Tm global_False": "Tm (30mer)",
            "Percent Peptide": "percent peptide ",
        }

        for i in range(df.shape[0]):
            thisfeat = df["Feature name"].iloc[i]
            if thisfeat in feature_dictionary:
                df["Feature name"].iloc[i] = feature_dictionary[thisfeat]

        return df


if __name__ == "__main__":
    # get_thirty_one_mer_data()

    V = "1"
    if V == "1":
        HUMAN_DATA = pd.read_excel("data/V1_data.xlsx", sheetname=0, index_col=[0, 1])
        MOUSE_DATA = pd.read_excel("data/V1_data.xlsx", sheetname=1, index_col=[0, 1])
        X, Y = combine_organisms(HUMAN_DATA, MOUSE_DATA)
        X.to_pickle("../data/X.pd")  # sequence features (i.e. inputs to prediction)
        Y.to_pickle(
            "../data/Y.pd"
        )  # cell-averaged ranks, plus more (i.e. possible targets for prediction)
        print("done writing to file")
    elif V == "2":
        # this is now all in predict.py
        pass
    elif V == "0":
        pass
