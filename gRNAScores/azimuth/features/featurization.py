from itertools import product
from sys import maxsize
from time import time
from typing import List

import Bio.Seq as Seq
import Bio.SeqUtils as SeqUtil
import Bio.SeqUtils.MeltingTemp as Tm
import numpy as np
import pandas as pd
from dill import dump
from sklearn.preprocessing import OneHotEncoder, LabelEncoder

from gRNAScores.azimuth import util
from gRNAScores.azimuth.features.microhomology import compute_score


def featurize_data(
    data, learn_options, Y, gene_position, pam_audit=True, length_audit=True, quiet=True
):
    """
    assumes that data contains the 30mer
    returns set of features from which one can make a kernel for each one
    """

    if np.any(data["30mer"].str.len() != 30):
        raise AssertionError(f"should only have sequences 30 nt long")

    if not quiet:
        print("Constructing features...")
    time_0 = time()

    feature_sets = {}

    if learn_options["nuc_features"]:
        # spectrum kernels (position-independent) and weighted degree kernels (position-dependent)
        get_all_order_nuc_features(
            data["30mer"],
            feature_sets,
            learn_options,
            learn_options["order"],
            max_index_to_use=30,
            quiet=quiet,
        )

    check_feature_set(feature_sets)

    if learn_options["gc_features"]:
        gc_above_10, gc_below_10, gc_count = gc_features(data, length_audit)
        feature_sets["gc_above_10"] = pd.DataFrame(gc_above_10)
        feature_sets["gc_below_10"] = pd.DataFrame(gc_below_10)
        feature_sets["gc_count"] = pd.DataFrame(gc_count)

    if learn_options["include_gene_position"]:
        for item in gene_position.columns:
            set_name = item
            feature_sets[set_name] = pd.DataFrame(gene_position[item])
        feature_sets["Percent Peptide <50%"] = feature_sets["Percent Peptide"] < 50
        feature_sets["Percent Peptide <50%"]["Percent Peptide <50%"] = feature_sets[
            "Percent Peptide <50%"
        ].pop("Percent Peptide")

    if learn_options["include_gene_effect"]:
        print("including gene effect")
        gene_names = Y["Target gene"]
        enc = OneHotEncoder()
        label_encoder = LabelEncoder()
        label_encoder.fit(gene_names)
        one_hot_genes = np.array(
            enc.fit_transform(label_encoder.transform(gene_names)[:, None]).todense()
        )
        feature_sets["gene effect"] = pd.DataFrame(
            one_hot_genes,
            columns=[f"gene_{i}" for i in range(one_hot_genes.shape[1])],
            index=gene_names.index,
        )

    if learn_options["include_known_pairs"]:
        feature_sets["known pairs"] = pd.DataFrame(Y["test"])

    if learn_options["include_NGGX_interaction"]:
        feature_sets["NGGX"] = NGGX_interaction_feature(data, pam_audit)

    if learn_options["include_Tm"]:
        feature_sets["Tm"] = Tm_feature(data, pam_audit, learn_options=None)

    if learn_options["include_sgRNAscore"]:
        feature_sets["sgRNA Score"] = pd.DataFrame(data["sgRNA Score"])

    if learn_options["include_drug"]:
        drug_names = Y.index.get_level_values("drug").tolist()
        enc = OneHotEncoder()
        label_encoder = LabelEncoder()
        label_encoder.fit(drug_names)
        one_hot_drugs = np.array(
            enc.fit_transform(label_encoder.transform(drug_names)[:, None]).todense()
        )
        feature_sets["drug"] = pd.DataFrame(
            one_hot_drugs,
            columns=[f"drug_{i:d}" for i in range(one_hot_drugs.shape[1])],
            index=drug_names,
        )

    if learn_options["include_strand"]:
        feature_sets["Strand effect"] = (pd.DataFrame(data["Strand"]) == "sense") * 1

    if learn_options["include_gene_feature"]:
        feature_sets["gene features"] = gene_feature(Y)
        # feature_sets["gene features"] = gene_feature(Y, data, learn_options)

    if learn_options["include_gene_guide_feature"] > 0:
        tmp_feature_sets = gene_guide_feature(Y, data, learn_options)
        for key in tmp_feature_sets:
            feature_sets[key] = tmp_feature_sets[key]

    if learn_options["include_microhomology"]:
        feature_sets["microhomology"] = get_micro_homology_features(
            Y["Target gene"], data
        )

    time_1 = time()
    if not quiet:
        print(
            f"\t\tElapsed time for constructing features is {time_1 - time_0:.2f} seconds"
        )

    check_feature_set(feature_sets)

    if learn_options["normalize_features"]:
        print(
            f"should not be here as doesn't make sense when we make one-off predictions, "
            f"but could make sense for internal model comparisons when using regularized models"
        )
        feature_sets = normalize_feature_sets(feature_sets)
        check_feature_set(feature_sets)

    return feature_sets


def check_feature_set(feature_sets):
    """
    Ensure the # of people is the same in each feature set
    """
    if feature_sets == {}:
        raise AssertionError("no feature sets present")

    N = None
    for ft in feature_sets:
        N2 = feature_sets[ft].shape[0]
        if N is None:
            N = N2
        else:
            if N < 1:
                raise AssertionError("should be at least one individual")
            if N != N2:
                raise AssertionError(
                    "# of individuals do not match up across feature sets"
                )

    for item in feature_sets:
        if np.any(np.isnan(feature_sets[item])):
            raise Exception(f"found Nan in set {item}")


def NGGX_interaction_feature(data, pam_audit=True):
    """
    assuming 30-mer, grab the NGGX _ _ positions, and make a one-hot
    encoding of the NX nucleotides yielding 4x4=16 features
    """
    sequence = data["30mer"].values
    feat_NX = pd.DataFrame()
    # check that GG is where we think
    for seq in sequence:
        if pam_audit and seq[25:27] != "GG":
            raise Exception(f"expected GG but found {seq[25 :27]}")
        NX = seq[24] + seq[27]
        NX_onehot = nucleotide_features(
            NX, order=2, feature_type="pos_dependent", max_index_to_use=2, prefix="NGGX"
        )
        feat_NX = pd.concat([feat_NX, NX_onehot], axis=1, sort=True)
    return feat_NX.T


def get_all_order_nuc_features(
    data,
    feature_sets,
    learn_options,
    maxorder,
    max_index_to_use,
    prefix="",
    quiet=False,
):
    for order in range(1, maxorder + 1):
        if not quiet:
            print(f"\t\tconstructing order {order} features")
        nuc_features_pd, nuc_features_pi = apply_nucleotide_features(
            data,
            order,
            include_pos_independent=True,
            max_index_to_use=max_index_to_use,
            prefix=prefix,
        )
        feature_sets[f"{prefix}_nuc_pd_Order{order:d}"] = nuc_features_pd
        if learn_options["include_pi_nuc_feat"]:
            feature_sets[f"{prefix}_nuc_pi_Order{order:d}"] = nuc_features_pi
        check_feature_set(feature_sets)

        if not quiet:
            print("\t\t\t\t\t\t\tdone")


def countGC(s, length_audit=True):
    """
    GC content for only the 20mer, as per the Doench paper/code
    """
    if length_audit:
        if len(s) != 30:
            raise AssertionError("seems to assume 30mer")
    return len(s[4:24].replace("A", "").replace("T", ""))


def SeqUtilFeatures(data):
    """
    assuming '30-mer'is a key
    get melting temperature features from:
        0-the 30-mer ("global Tm")
        1-the Tm (melting temperature) of the DNA:RNA hybrid from positions 16 - 20 of the sgRNA,
        i.e. the 5nts immediately proximal of the NGG PAM
        2-the Tm of the DNA:RNA hybrid from position 8 - 15 (i.e. 8 nt)
        3-the Tm of the DNA:RNA hybrid from position 3 - 7  (i.e. 5 nt)
    """
    sequence = data["30mer"].values
    num_features = 1
    featarray = np.ones((sequence.shape[0], num_features))
    for i, seq in enumerate(sequence):
        if len(seq) != 30:
            raise AssertionError("seems to assume 30mer")
        featarray[i, 0] = SeqUtil.molecular_weight(str(seq))

    feat = pd.DataFrame(pd.DataFrame(featarray))
    return feat


def organism_feature(data):
    """
    Human vs. mouse
    """
    organism = np.array(data["Organism"].values)
    feat = pd.DataFrame(pd.DataFrame(organism))
    import pdb

    pdb.set_trace()
    return feat


def get_micro_homology_features(gene_names, X):
    # originally was flipping the guide itself as necessary, but now flipping the gene instead

    print("building microhomology features")
    feat = pd.DataFrame(index=X.index)
    feat["mh_score"] = ""
    feat["oof_score"] = ""

    # number of nulceotides to take to the left and right of the guide
    K_MER_LENGTH_LEFT = 9
    K_MER_LENGTH_RIGHT = 21
    for gene in gene_names.unique():
        gene_seq = Seq.Seq(util.get_gene_sequence(gene)).reverse_complement()
        guide_inds = np.where(gene_names.values == gene)[0]
        print(f"getting microhomology for all {len(guide_inds)} guides in gene {gene}")
        for ps in guide_inds:
            guide_seq = Seq.Seq(X["30mer"][ps])
            strand = X["Strand"][ps]
            if strand == "sense":
                gene_seq = gene_seq.reverse_complement()
            # figure out the sequence to the left and right of this guide, in the gene
            ind = gene_seq.find(guide_seq)
            if ind == -1:
                gene_seq = gene_seq.reverse_complement()
                ind = gene_seq.find(guide_seq)
            else:
                pass
            if ind == -1:
                mh_score = 0
                oof_score = 0
            else:

                if gene_seq[ind : (ind + len(guide_seq))] != guide_seq:
                    raise AssertionError("match not right")

                left_win = gene_seq[(ind - K_MER_LENGTH_LEFT) : ind]
                right_win = gene_seq[
                    (ind + len(guide_seq)) : (ind + len(guide_seq) + K_MER_LENGTH_RIGHT)
                ]
                if len(left_win.tostring()) != K_MER_LENGTH_LEFT:
                    raise AssertionError()
                if len(right_win.tostring()) != K_MER_LENGTH_RIGHT:
                    raise AssertionError()
                sixtymer = str(left_win) + str(guide_seq) + str(right_win)
                if len(sixtymer) != 60:
                    raise AssertionError("should be of length 60")
                mh_score, oof_score = compute_score(sixtymer)

            feat.ix[ps, "mh_score"] = mh_score
            feat.ix[ps, "oof_score"] = oof_score
        print(f"computed microhomology of {str(gene)}")

    return pd.DataFrame(feat, dtype="float")


def local_gene_seq_features(gene_names, learn_options, X):
    print(f"building local gene sequence features")
    feat = pd.DataFrame(index=X.index)
    feat["gene_left_win"] = ""
    feat["gene_right_win"] = ""

    # number of nulceotides to take to the left and right of the guide
    k_mer_length = learn_options["include_gene_guide_feature"]
    for gene in gene_names.unique():
        gene_seq = Seq.Seq(util.get_gene_sequence(gene)).reverse_complement()
        for ps in np.where(gene_names.values == gene)[0]:
            guide_seq = Seq.Seq(X["30mer"][ps])
            strand = X["Strand"][ps]
            if strand == "sense":
                guide_seq = guide_seq.reverse_complement()
            # figure out the sequence to the left and right of this guide, in the gene
            ind = gene_seq.find(guide_seq)
            if ind == -1:
                if ind == -1:
                    raise AssertionError("could not find guide in gene")
            if gene_seq[ind : (ind + len(guide_seq))] != guide_seq:
                raise AssertionError("match not right")
            left_win = gene_seq[(ind - k_mer_length) : ind]
            right_win = gene_seq[
                (ind + len(guide_seq)) : (ind + len(guide_seq) + k_mer_length)
            ]

            if strand == "antisense":
                # it's arbitrary which of sense and anti-sense we flip, we just want
                # to keep them in the same relative alphabet/direction
                left_win = left_win.reverse_complement()
                right_win = right_win.reverse_complement()
            if left_win.tostring() == "":
                raise AssertionError(f"k_mer_context, {k_mer_length}, is too large")
            if len(left_win) != len(right_win):
                raise AssertionError(f"k_mer_context, {k_mer_length}, is too large")
            feat.ix[ps, "gene_left_win"] = left_win.tostring()
            feat.ix[ps, "gene_right_win"] = right_win.tostring()
        print(f"featurizing local context of {gene}")

    feature_sets = {}
    get_all_order_nuc_features(
        feat["gene_left_win"],
        feature_sets,
        learn_options,
        learn_options["order"],
        max_index_to_use=maxsize,
        prefix="gene_left_win",
    )
    get_all_order_nuc_features(
        feat["gene_right_win"],
        feature_sets,
        learn_options,
        learn_options["order"],
        max_index_to_use=maxsize,
        prefix="gene_right_win",
    )
    return feature_sets


def gene_feature(Y):
    """
    Things like the sequence of the gene, the DNA Tm of the gene, etc.
    """

    gene_names = Y["Target gene"]

    gene_length = np.zeros((gene_names.values.shape[0], 1))
    gc_content = np.zeros((gene_names.shape[0], 1))
    temperature = np.zeros((gene_names.shape[0], 1))
    molecular_weight = np.zeros((gene_names.shape[0], 1))

    for gene in gene_names.unique():
        seq = util.get_gene_sequence(gene)
        gene_length[gene_names.values == gene] = len(seq)
        gc_content[gene_names.values == gene] = SeqUtil.GC(seq)
        temperature[gene_names.values == gene] = Tm.Tm_staluc(seq, rna=False)
        molecular_weight[gene_names.values == gene] = SeqUtil.molecular_weight(
            seq, "DNA"
        )

    everything = np.concatenate(
        (gene_length, gc_content, temperature, molecular_weight), axis=1
    )
    df = pd.DataFrame(
        data=everything,
        index=gene_names.index,
        columns=[
            "gene length",
            "gene GC content",
            "gene temperature",
            "gene molecular weight",
        ],
    )
    return df


def gene_guide_feature(Y, X, learn_options):
    # features, which are related to parts of the gene-local to the guide, and
    # possibly incorporating the guide or interactions with it

    # expensive, so pickle if necessary
    gene_file = (
        f"../data/gene_seq_feat_V{learn_options['V']}_km"
        f"{learn_options['include_gene_guide_feature']}.ord{learn_options['order']}.pickle"
    )

    # if False:  # os.path.isfile(gene_file): #while debugging, comment out
    #     print(f"loading local gene seq feats from file {gene_file}")
    #     with open(gene_file, "rb") as file_to_open:
    #         feature_sets = load(file_to_open)
    # else:
    feature_sets = local_gene_seq_features(Y["Target gene"], learn_options, X)
    print(f"writing local gene seq feats to file {gene_file}")
    with open(gene_file, "wb") as file_to_open:
        dump(feature_sets, file_to_open)

    return feature_sets


def gc_cont(seq):
    return (seq.count("G") + seq.count("C")) / float(len(seq))


def Tm_feature(data, pam_audit=True, learn_options=None):
    """
    assuming '30-mer'is a key
    get melting temperature features from:
        0-the 30-mer ("global Tm")
        1-the Tm (melting temperature) of the DNA:RNA hybrid from positions 16 - 20 of the sgRNA,
        i.e. the 5nts immediately proximal of the NGG PAM
        2-the Tm of the DNA:RNA hybrid from position 8 - 15 (i.e. 8 nt)
        3-the Tm of the DNA:RNA hybrid from position 3 - 7  (i.e. 5 nt)
    """

    if learn_options is None or "Tm segments" not in learn_options:
        segments = [(19, 24), (11, 19), (6, 11)]
    else:
        segments = learn_options["Tm segments"]

    sequence = data["30mer"].values
    featarray = np.ones((sequence.shape[0], 4))

    for i, seq in enumerate(sequence):
        if pam_audit and seq[25:27] != "GG":
            raise Exception(f"expected GG but found {seq[25:27]}")
        rna = False
        featarray[i, 0] = Tm.Tm_staluc(seq, rna=rna)  # 30mer Tm
        featarray[i, 1] = Tm.Tm_staluc(
            seq[segments[0][0] : segments[0][1]], rna=rna
        )  # 5nts immediately proximal of the NGG PAM
        featarray[i, 2] = Tm.Tm_staluc(
            seq[segments[1][0] : segments[1][1]], rna=rna
        )  # 8-mer
        featarray[i, 3] = Tm.Tm_staluc(
            seq[segments[2][0] : segments[2][1]], rna=rna
        )  # 5-mer

    feat = pd.DataFrame(
        featarray,
        index=data.index,
        columns=[
            f"Tm global_{rna}",
            f"5mer_end_{rna}",
            f"8mer_middle_{rna}",
            f"5mer_start_{rna}",
        ],
    )

    return feat


def gc_features(data, audit=True):
    gc_count = data["30mer"].apply(lambda seq: countGC(seq, audit))
    gc_count.name = "GC count"
    gc_above_10 = (gc_count > 10) * 1
    gc_above_10.name = "GC > 10"
    gc_below_10 = (gc_count < 10) * 1
    gc_below_10.name = "GC < 10"
    return gc_above_10, gc_below_10, gc_count


def normalize_features(data, axis):
    """
    input: pd.DataFrame of dtype=np.float64 array, of dimensions
    mean-center, and unit variance each feature
    """
    data -= data.mean(axis)
    data /= data.std(axis)
    # remove rows with NaNs
    data = data.dropna(1)
    if np.any(np.isnan(data.values)):
        raise Exception("found NaN in normalized features")
    return data


def apply_nucleotide_features(
    seq_data_frame, order, include_pos_independent, max_index_to_use, prefix=""
):
    if include_pos_independent:
        feat_pd = seq_data_frame.apply(
            nucleotide_features, args=(order, max_index_to_use, prefix, "pos_dependent")
        )
        feat_pi = seq_data_frame.apply(
            nucleotide_features,
            args=(order, max_index_to_use, prefix, "pos_independent"),
        )
        if np.any(np.isnan(feat_pd)):
            raise AssertionError(
                "nans here can arise from sequences of different lengths"
            )
        if np.any(np.isnan(feat_pi)):
            raise AssertionError(
                "nans here can arise from sequences of different lengths"
            )
    else:
        feat_pd = seq_data_frame.apply(
            nucleotide_features, args=(order, max_index_to_use, prefix, "pos_dependent")
        )
        if np.any(np.isnan(feat_pd)):
            raise AssertionError("found nan in feat_pd")
        feat_pi = None
    return feat_pd, feat_pi


def get_alphabet(order, raw_alphabet=("A", "T", "C", "G")):
    alphabet = ["".join(i) for i in product(raw_alphabet, repeat=order)]
    return alphabet


def nucleotide_features(
    s,
    order,
    max_index_to_use,
    prefix="",
    feature_type="all",
    raw_alphabet=("A", "T", "C", "G"),
):
    """
    compute position-specific order-mer features for the 4-letter alphabet
    (e.g. for a sequence of length 30, there are 30*4 single nucleotide features
          and (30-1)*4^2=464 double nucleotide features
    """
    if not feature_type in ["all", "pos_independent", "pos_dependent"]:
        raise AssertionError("Unknown feature_type")
    if max_index_to_use <= len(s):
        max_index_to_use = len(s)

    if max_index_to_use is not None:
        s = s[:max_index_to_use]
    # s = s[:30] #cut-off at thirty to clean up extra data that they accidentally left in,
    # nd were instructed to ignore in this way
    alphabet: List[str] = get_alphabet(order, raw_alphabet=raw_alphabet)
    features_pos_dependent = np.zeros(len(alphabet) * (len(s) - (order - 1)))
    features_pos_independent = np.zeros(np.power(len(raw_alphabet), order))

    index_dependent: List[str] = []
    index_independent: List[str] = []

    for position in range(0, len(s) - order + 1, 1):
        for l in alphabet:
            index_dependent.append(f"{prefix}{l}_{position:d}")

    for l in alphabet:
        index_independent.append(f"{prefix}{l}")

    for position in range(0, len(s) - order + 1, 1):
        nucl: object = s[position : position + order]
        features_pos_dependent[alphabet.index(nucl) + (position * len(alphabet))] = 1.0
        features_pos_independent[alphabet.index(nucl)] += 1.0

        # this is to check that the labels in the pd df actually match the nucl and position
        if (
            index_dependent[alphabet.index(nucl) + (position * len(alphabet))]
            != f"{prefix}{nucl}_{position:d}"
        ):
            raise AssertionError()
        if index_independent[alphabet.index(nucl)] != f"{prefix}{nucl}":
            raise AssertionError()

    if np.any(np.isnan(features_pos_dependent)):
        raise Exception("found nan features in features_pos_dependent")
    if np.any(np.isnan(features_pos_independent)):
        raise Exception("found nan features in features_pos_independent")

    if feature_type in ("all", "pos_independent"):
        if feature_type == "all":
            res = pd.Series(features_pos_dependent, index=index_dependent).append(
                pd.Series(features_pos_independent, index=index_independent)
            )
            if np.any(np.isnan(res.values)):
                raise AssertionError()
        else:
            res = pd.Series(features_pos_independent, index=index_independent)
            if np.any(np.isnan(res.values)):
                raise AssertionError()
    else:
        res = pd.Series(features_pos_dependent, index=index_dependent)

    if np.any(np.isnan(res.values)):
        raise AssertionError()
    return res


def nucleotide_features_dictionary(prefix=""):
    seqname = ["-4", "-3", "-2", "-1"]
    seqname.extend([str(i) for i in range(1, 21)])
    seqname.extend(["N", "G", "G", "+1", "+2", "+3"])

    orders = [1, 2, 3]
    sequence = 30
    feature_names_dep = []
    feature_names_indep = []
    index_dependent = []
    index_independent = []

    for order in orders:
        raw_alphabet = ["A", "T", "C", "G"]
        alphabet = ["".join(i) for i in product(raw_alphabet, repeat=order)]
        features_pos_dependent = np.zeros(len(alphabet) * (sequence - (order - 1)))
        features_pos_independent = np.zeros(np.power(len(raw_alphabet), order))

        index_dependent.extend(
            [
                f"{prefix}_pd.Order{order}_P{i}"
                for i in range(len(features_pos_dependent))
            ]
        )
        index_independent.extend(
            [
                f"{prefix}_pi.Order{order}_P{i}"
                for i in range(len(features_pos_independent))
            ]
        )

        for pos in range(sequence - (order - 1)):
            for letter in alphabet:
                feature_names_dep.append(f"{letter}_{seqname[pos]}")

        for letter in alphabet:
            feature_names_indep.append(f"{letter}")

        if len(feature_names_indep) != len(index_independent):
            raise AssertionError()
        if len(feature_names_dep) != len(index_dependent):
            raise AssertionError()

    index_all = index_dependent + index_independent
    feature_all = feature_names_dep + feature_names_indep

    return dict(zip(index_all, feature_all))


def normalize_feature_sets(feature_sets):
    """
    zero-mean, unit-variance each feature within each set
    :type feature_sets: List[np.ndarray]
    """

    print("Normalizing features...")
    t1 = time()

    new_feature_sets = {}
    for fset in feature_sets:
        new_feature_sets[fset] = normalize_features(feature_sets[fset], axis=0)
        if np.any(np.isnan(new_feature_sets[fset].values)):
            raise Exception(f"found Nan feature values in set={fset}")
        if new_feature_sets[fset].shape[1] <= 0:
            raise AssertionError("0 columns of features")
    t2 = time()
    print(f"\t\tElapsed time for normalizing features is {(t2 - t1):.2f} seconds")

    return new_feature_sets
