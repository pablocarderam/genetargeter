from os import path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from . import util

CUR_DIR = path.dirname(path.abspath(__file__))


def from_custom_file(data_file: str) -> tuple:
    # use semantics of when we load V2 data
    print(f"Loading inputs to predict from {data_file}")
    data = pd.read_csv(data_file)

    mandatory_columns = [
        "30mer",
        "Target gene",
        "Percent Peptide",
        "Amino Acid Cut position",
    ]
    for col in mandatory_columns:
        if col not in data.columns:
            raise AssertionError(
                f"inputs for prediction must include these columns: "
                f"{mandatory_columns}"
            )

    x_df = pd.DataFrame(data)
    x_df["30mercopy"] = x_df["30mer"]
    x_df = x_df.set_index(["30mer", "Target gene"])
    x_df["30mer"] = x_df["30mercopy"]
    x_df.index.names = ["Sequence", "Target"]
    x_df["drug"] = [f"dummydrug{i}" for i in range(x_df.shape[0])]
    x_df = x_df.set_index("drug", append=True)

    y_df = None
    gene_position = x_df[["Percent Peptide", "Amino Acid Cut position"]]
    target_genes = np.unique(x_df.index.levels[1])

    return x_df, y_df, gene_position, target_genes


def from_file(data_file: str, data_file2: str, learn_options: dict) -> tuple:
    if learn_options["V"] == 1:  # from Nature Biotech paper

        print(f"loading V{learn_options['V']} data")

        if learn_options["weighted"] is not None:
            raise AssertionError("not supported for V1 data")
        _, gene_position, target_genes, x_df, y_df = read_V1_data(
            data_file, learn_options
        )

        learn_options["binary target name"] = "average threshold"
        learn_options["rank-transformed target name"] = "average rank"
        learn_options["raw target name"] = "average activity"

    elif learn_options["V"] == 2:  # from Nov 2014, hot off the machines
        x_df, _, target_genes, y_df, gene_position = read_V2_data(
            data_file, learn_options
        )

        # check that data is consistent with sgRNA score
        xx = x_df["sgRNA Score"].values
        yy = y_df["score_drug_gene_rank"].values
        rr, _ = pearsonr(xx, yy)
        if rr <= 0:
            raise AssertionError(
                "data processing has gone wrong as correlation with previous "
                "predictions is negative"
            )

    elif (
        learn_options["V"] == 3
    ):  # merge of V1 and V2--this is what is used for the final model
        # these are relative to the V2 data, and V1 will be made to automatically match
        learn_options["binary target name"] = "score_drug_gene_threshold"
        learn_options["rank-transformed target name"] = "score_drug_gene_rank"
        learn_options["raw target name"] = None

        x_df, y_df, gene_position, target_genes = mergeV1_V2(
            data_file, data_file2, learn_options
        )

    elif learn_options["V"] == 4:  # merge of V1 and V2 and the Xu et al data
        # these are relative to the V2 data, and V1 and Xu et al. will be made
        # to automatically match
        learn_options["binary target name"] = "score_drug_gene_threshold"
        learn_options["rank-transformed target name"] = "score_drug_gene_rank"
        learn_options["raw target name"] = None

        x_df, y_df, gene_position, target_genes = merge_all(
            data_file, data_file2, learn_options
        )

    elif learn_options["V"] == 5:
        raise Exception(
            "The from_file() function is attempting to learn using the xu_et_al data.  "
            "This data is no longer available with Azimuth."
        )

    # truncate down to 30--some data sets gave us more.
    x_df["30mer"] = x_df["30mer"].apply(lambda x: x[0:30])

    return x_df, y_df, gene_position, target_genes


def set_V2_target_names(learn_options: dict) -> dict:
    if "binary target name" not in learn_options:
        learn_options["binary target name"] = "score_drug_gene_threshold"
    if "rank-transformed target name" not in learn_options:
        learn_options["rank-transformed target name"] = "score_drug_gene_rank"
    learn_options["raw target name"] = "score"
    return learn_options


def combine_organisms(human_data: pd.DataFrame, mouse_data: pd.DataFrame) -> tuple:
    # 'Target' is the column name, 'CD13' are some rows in that column
    # xs slices through the pandas data frame to return another one
    cd13 = human_data.xs("CD13", level="Target", drop_level=False)
    # y_names are column names, cd13 is a pd object
    x_cd13, y_cd13 = util.get_data(cd13, y_names=["NB4 CD13", "TF1 CD13"])
    cd33 = human_data.xs("CD33", level="Target", drop_level=False)
    x_cd33, y_cd33 = util.get_data(
        cd33, y_names=["MOLM13 CD33", "TF1 CD33", "NB4 CD33"]
    )
    cd15 = human_data.xs("CD15", level="Target", drop_level=False)
    x_cd15, y_cd15 = util.get_data(cd15, y_names=["MOLM13 CD15"])

    mouse_x = pd.DataFrame()
    mouse_y = pd.DataFrame()
    for k in mouse_data.index.levels[1]:
        # is k the gene
        x_df, y_df = util.get_data(
            mouse_data.xs(k, level="Target", drop_level=False),
            ["On-target Gene"],
            target_gene=k,
            organism="mouse",
        )
        mouse_x = pd.concat([mouse_x, x_df], axis=0)
        mouse_y = pd.concat([mouse_y, y_df], axis=0)

    x_df = pd.concat([x_cd13, x_cd15, x_cd33, mouse_x], axis=0, sort=True)
    y_df = pd.concat([y_cd13, y_cd15, y_cd33, mouse_y], axis=0, sort=True)

    return x_df, y_df


def read_V1_data(
    data_file: str = None,
    learn_options: dict = None,
    aml_file: str = f"{CUR_DIR}/data/V1_suppl_data.txt",
) -> tuple:
    if data_file is None:
        data_file = CUR_DIR + "/data/V1_data.xlsx"
    human_data = pd.read_excel(data_file, sheetname=0, index_col=[0, 1])
    mouse_data = pd.read_excel(data_file, sheetname=1, index_col=[0, 1])
    x_df, y_df = combine_organisms(human_data, mouse_data)

    # get position within each gene, then join and re-order
    # note that 11 missing guides we were told to ignore
    annotations = pd.read_csv(aml_file, delimiter="\t", index_col=[0, 4])
    annotations.index.names = x_df.index.names
    gene_position = pd.merge(
        x_df, annotations, how="inner", left_index=True, right_index=True
    )
    gene_position = util.impute_gene_position(gene_position)
    gene_position = gene_position[
        ["Amino Acid Cut position", "Nucleotide cut position", "Percent Peptide"]
    ]
    y_df = y_df.loc[gene_position.index]
    x_df = x_df.loc[gene_position.index]

    y_df[
        "test"
    ] = (
        1
    )  # for bookkeeping to keep consistent with V2 which uses this for "extra pairs"

    target_genes = y_df["Target gene"].unique()

    y_df.index.names = ["Sequence", "Target gene"]

    if not x_df.index.equals(y_df.index):
        raise AssertionError(
            "The index of x_df is different from the index of y_df "
            "(this can cause inconsistencies/random performance later on)"
        )

    if learn_options is not None and learn_options["flipV1target"]:
        print(
            "************************************************************************\n"
            "*****************MATCHING DOENCH CODE (DEBUG MODE)**********************\n"
            "************************************************************************"
        )
        # normally it is:y_df['average threshold'] =y_df['average rank'] > 0.8, where
        # 1s are good guides, 0s are not
        y_df["average threshold"] = y_df["average rank"] < 0.2  # 1s are bad guides
        print("press c to continue")
        import pdb

        pdb.set_trace()

    return annotations, gene_position, target_genes, x_df, y_df


def read_V2_data(
    data_file: str = None, learn_options: dict = None, verbose: bool = True
) -> tuple:
    if data_file is None:
        data_file = CUR_DIR + "/data/V2_data.xlsx"

    data = pd.read_excel(
        data_file,
        sheetname="ResultsFiltered",
        skiprows=range(0, 6 + 1),
        index_col=[0, 4],
    )
    # grab data relevant to each of three drugs, which exludes some genes
    # note gene MED12 has two drugs, all others have at most one
    x_df = pd.DataFrame()

    # This comes from the "Pairs" tab in their excel sheet,
    # note HPRT/HPRT1 are same thing, and also PLX_2uM/PLcX_2uM
    known_pairs = {
        "AZD_200nM": ["CCDC101", "MED12", "TADA2B", "TADA1"],
        "6TG_2ug/mL": ["HPRT1"],
        "PLX_2uM": ["CUL3", "NF1", "NF2", "MED12"],
    }

    drugs_to_genes = {
        "AZD_200nM": ["CCDC101", "MED12", "TADA2B", "TADA1"],
        "6TG_2ug/mL": ["HPRT1"],
        "PLX_2uM": ["CUL3", "NF1", "NF2", "MED12"],
    }

    if learn_options is not None:
        if learn_options["extra pairs"] or learn_options["all pairs"]:
            raise AssertionError(
                "extra pairs and all pairs options (in learn_options) can't be "
                "active simultaneously."
            )

        if learn_options["extra pairs"]:
            drugs_to_genes["AZD_200nM"].extend(["CUL3", "NF1", "NF2"])
        elif learn_options["all pairs"]:
            drugs_to_genes["AZD_200nM"].extend(["HPRT1", "CUL3", "NF1", "NF2"])
            drugs_to_genes["PLX_2uM"].extend(["HPRT1", "CCDC101", "TADA2B", "TADA1"])
            drugs_to_genes["6TG_2ug/mL"].extend(
                ["CCDC101", "MED12", "TADA2B", "TADA1", "CUL3", "NF1", "NF2"]
            )

    count = 0
    for drug in drugs_to_genes:
        genes = drugs_to_genes[drug]
        for gene in genes:
            xtmp = data.copy().xs(gene, level="Target gene", drop_level=False)
            xtmp["drug"] = drug
            xtmp["score"] = xtmp[
                drug
            ].copy()  # grab the drug results that are relevant for this gene

            if gene in known_pairs[drug]:
                xtmp["test"] = 1.0
            else:
                xtmp["test"] = 0.0

            count = count + xtmp.shape[0]
            x_df = pd.concat([x_df, xtmp], axis=0)
            if verbose:
                print(
                    f"Loaded {xtmp.shape[0]} samples for gene {gene} "
                    f"\ttotal number of samples: {count}"
                )

    # create new index that includes the drug
    x_df = x_df.set_index("drug", append=True)

    y_df = pd.DataFrame(x_df.pop("score"))
    y_df.columns.names = ["score"]

    test_gene = pd.DataFrame(x_df.pop("test"))
    target = pd.DataFrame(
        x_df.index.get_level_values("Target gene").values,
        index=y_df.index,
        columns=["Target gene"],
    )
    y_df = pd.concat((y_df, target, test_gene), axis=1)
    target_genes = y_df["Target gene"].unique()
    gene_position = x_df[["Percent Peptide", "Amino Acid Cut position"]].copy()

    # convert to ranks for each (gene, drug combo)
    # flip = True
    y_rank = pd.DataFrame()
    y_threshold = pd.DataFrame()
    y_quant = pd.DataFrame()
    for drug in drugs_to_genes:
        gene_list = drugs_to_genes[drug]
        for gene in gene_list:
            ytmp = pd.DataFrame(
                y_df.xs((gene, drug), level=["Target gene", "drug"], drop_level=False)[
                    "score"
                ]
            )
            y_ranktmp, _, y_thresholdtmp, y_quanttmp = util.get_ranks(
                ytmp, thresh=0.8, prefix="score_drug_gene", flip=False
            )
            # np.unique(y_rank.values-y_rank_raw.values)
            y_rank = pd.concat((y_rank, y_ranktmp), axis=0)
            y_threshold = pd.concat((y_threshold, y_thresholdtmp), axis=0)
            y_quant = pd.concat((y_quant, y_quanttmp), axis=0)

    yall = pd.concat((y_rank, y_threshold, y_quant), axis=1)
    y_df = pd.merge(y_df, yall, how="inner", left_index=True, right_index=True)

    # convert also by drug only, irrespective of gene
    y_rank = pd.DataFrame()
    y_threshold = pd.DataFrame()
    y_quant = pd.DataFrame()
    for drug in drugs_to_genes:
        ytmp = pd.DataFrame(y_df.xs(drug, level="drug", drop_level=False)["score"])
        y_ranktmp, _, y_thresholdtmp, y_quanttmp = util.get_ranks(
            ytmp, thresh=0.8, prefix="score_drug", flip=False
        )
        # np.unique(y_rank.values-y_rank_raw.values)
        y_rank = pd.concat((y_rank, y_ranktmp), axis=0)
        y_threshold = pd.concat((y_threshold, y_thresholdtmp), axis=0)
        y_quant = pd.concat((y_quant, y_quanttmp), axis=0)

    yall = pd.concat((y_rank, y_threshold, y_quant), axis=1)
    y_df = pd.merge(y_df, yall, how="inner", left_index=True, right_index=True)

    gene_position = util.impute_gene_position(gene_position)

    if learn_options is not None and learn_options["weighted"] == "variance":
        print("computing weights from replicate variance...")
        # compute the variance across replicates so can use it as a weight
        data = pd.read_excel(
            data_file,
            sheetname="Normalized",
            skiprows=range(0, 6 + 1),
            index_col=[0, 4],
        )
        data.index.names = ["Sequence", "Target gene"]

        experiments = {
            "AZD_200nM": ["Deep 25", "Deep 27", "Deep 29 ", "Deep 31"],
            "6TG_2ug/mL": ["Deep 33", "Deep 35", "Deep 37", "Deep 39"],
            "PLX_2uM": ["Deep 49", "Deep 51", "Deep 53", "Deep 55"],
        }

        variance = None
        for drug in drugs_to_genes:
            data_tmp = data.iloc[
                data.index.get_level_values("Target gene").isin(drugs_to_genes[drug])
            ][experiments[drug]]
            data_tmp["drug"] = drug
            data_tmp = data_tmp.set_index("drug", append=True)
            data_tmp["variance"] = np.var(data_tmp.values, axis=1)
            if variance is None:
                variance = data_tmp["variance"].copy()
            else:
                variance = pd.concat((variance, data_tmp["variance"]), axis=0)

        orig_index = y_df.index.copy()
        y_df = pd.merge(
            y_df, pd.DataFrame(variance), how="inner", left_index=True, right_index=True
        )
        y_df = y_df.ix[orig_index]
        print("done.")

    # Make sure to keep this check last in this function
    if not x_df.index.equals(y_df.index):
        raise AssertionError(
            "The index of x_df is different from the index of y_df "
            "(this can cause inconsistencies/random performance later on)"
        )

    return x_df, drugs_to_genes, target_genes, y_df, gene_position


def merge_all(data_file: str, data_file2: str, learn_options: dict) -> tuple:
    x_df, y_df, gene_position, target_genes = mergeV1_V2(
        data_file, data_file2, learn_options
    )
    return x_df, y_df, gene_position, target_genes


def mergeV1_V2(data_file: str, data_file2: str, learn_options: dict) -> tuple:
    """
    ground_truth_label, etc. are taken to correspond to the V2 data,
    and then the V1 is appropriately matched based on semantics
    """
    if learn_options["include_strand"]:
        raise AssertionError("don't currently have 'Strand' column in V1 data")

    _, gene_position1, target_genes1, x_df1, y_df1 = read_V1_data(
        data_file, learn_options
    )
    x_df2, _, target_genes2, y_df2, gene_position2 = read_V2_data(data_file2)

    y_df1.rename(
        columns={"average rank": learn_options["rank-transformed target name"]},
        inplace=True,
    )
    y_df1.rename(
        columns={"average threshold": learn_options["binary target name"]}, inplace=True
    )

    # rename columns, and add a dummy "drug" to V1 so can join the data sets
    y_df1["drug"] = ["nodrug" for _ in range(y_df1.shape[0])]
    y_df1 = y_df1.set_index("drug", append=True)
    y_df1.index.names = ["Sequence", "Target gene", "drug"]

    y_cols_to_keep = np.unique(
        ["Target gene", "test", "score_drug_gene_rank", "score_drug_gene_threshold"]
    )

    y_df1 = y_df1[y_cols_to_keep]
    y_df2 = y_df2[y_cols_to_keep]

    x_df1["drug"] = ["nodrug" for _ in range(x_df1.shape[0])]
    x_df1 = x_df1.set_index("drug", append=True)

    x_cols_to_keep = ["30mer", "Strand"]
    x_df1 = x_df1[x_cols_to_keep]
    x_df2 = x_df2[x_cols_to_keep]

    gene_position1["drug"] = ["nodrug" for _ in range(gene_position1.shape[0])]
    gene_position1 = gene_position1.set_index("drug", append=True)
    gene_position1.index.names = ["Sequence", "Target gene", "drug"]
    cols_to_keep = ["Percent Peptide", "Amino Acid Cut position"]
    gene_position1 = gene_position1[cols_to_keep]
    gene_position2 = gene_position2[cols_to_keep]

    y_df = pd.concat((y_df1, y_df2), axis=0)
    x_df = pd.concat((x_df1, x_df2), axis=0)
    gene_position = pd.concat((gene_position1, gene_position2))

    # target_genes = target_genes1 + target_genes2
    target_genes = np.concatenate((target_genes1, target_genes2))

    save_to_file = False

    if save_to_file:
        y_df.index.names = ["Sequence", "Target", "drug"]
        if np.any(x_df.index.values != y_df.index.values):
            raise AssertionError("rows don't match up")

        onedupind = np.where(y_df.index.duplicated())[0][0]
        alldupind = np.where(
            y_df.index.get_level_values(0).values == y_df.index[onedupind][0]
        )[0]

        # arbitrarily set one of these to have "nodrug2" as the third level index
        # so that they are not repeated, and the joints therefore do not augment the data set
        if len(alldupind) != 2:
            raise AssertionError("expected only duplicates")
        newindex = y_df.index.tolist()
        newindex[onedupind] = (
            newindex[onedupind][0],
            newindex[onedupind][1],
            "nodrug2",
        )
        y_df.index = pd.MultiIndex.from_tuples(newindex, names=y_df.index.names)
        x_df.index = pd.MultiIndex.from_tuples(newindex, names=y_df.index.names)

        # there seems to be a duplicate index, and thus this increases the data set size,
        # so doing it the hacky way...
        x_and_y = pd.merge(x_df, y_df, how="inner", left_index=True, right_index=True)
        gene_position_tmp = gene_position.copy()
        gene_position_tmp.index.names = ["Sequence", "Target", "drug"]
        gene_position_tmp.index = pd.MultiIndex.from_tuples(
            newindex, names=y_df.index.names
        )
        x_and_y = pd.merge(
            x_and_y, gene_position_tmp, how="inner", left_index=True, right_index=True
        )

        # truncate to 30mers
        x_and_y["30mer"] = x_and_y["30mer"].apply(lambda x: x[0:30])
        x_and_y.to_csv(r"D:\Source\CRISPR\data\tmp\V3.csv")

    return x_df, y_df, gene_position, target_genes


def get_V1_genes(data_file=None, learn_options: dict = None) -> np.ndarray:
    _, _, target_genes, _, _ = read_V1_data(data_file, learn_options)
    return target_genes


def get_V2_genes(data_file: str = None) -> np.ndarray:
    _, _, target_genes, _, _ = read_V2_data(data_file, verbose=False)
    return target_genes


def get_V3_genes(data_fileV1: str = None, data_fileV2: str = None) -> np.ndarray:
    target_genes = np.concatenate(
        (get_V1_genes(data_fileV1), get_V2_genes(data_fileV2))
    )
    return target_genes


def get_mouse_genes(data_file: str = None) -> np.ndarray:
    _, _, _, x_df, _ = read_V1_data(data_file, learn_options=None)
    return x_df[x_df["Organism"] == "mouse"]["Target gene"].unique()


def get_human_genes(data_file: str = None) -> np.ndarray:
    _, _, _, x_df, _ = read_V1_data(data_file, learn_options=None)
    mouse_genes = x_df[x_df["Organism"] == "mouse"]["Target gene"].unique()
    all_genes = get_V3_genes(None, None)
    return np.setdiff1d(all_genes, mouse_genes)
