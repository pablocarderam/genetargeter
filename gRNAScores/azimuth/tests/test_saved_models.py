import unittest
from os import path

from numpy import array, allclose
from pandas import read_csv

from ..model_comparison import predict

dirname, filename = path.split(path.abspath(__file__))


class SavedModelTests(unittest.TestCase):
    """
    This unit test checks that the predictions for 1000 guides match the predictions we expected in Nov 2016.
    This unit test can fail due to randomness in the model (e.g. random seed, feature reordering).
    """

    def test_predictions_nopos(self):
        df = read_csv(path.join(dirname, "1000guides.csv"), index_col=0)
        predictions = predict(array(df["guide"].values), None, None)
        self.assertTrue(allclose(predictions, df["truth nopos"].values, atol=1e-1))

    def test_predictions_pos(self):
        df = read_csv(path.join(dirname, "1000guides.csv"), index_col=0)
        predictions = predict(
            array(df["guide"].values),
            array(df["AA cut"].values),
            array(df["Percent peptide"].values),
        )
        self.assertTrue(allclose(predictions, df["truth pos"].values, atol=1e-1))


if __name__ == "__main__":
    unittest.main()
