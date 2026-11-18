import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from paired import score, input_counts


class PairedTests(unittest.TestCase):
    def test_score(self):
        for correct, wrong, expected in [
            (4, 0, "correct"), (4, 1, "wrong"), (3, 0, "no_solve"),
            (0, 0, "no_solve"), (0, 1, "wrong"),
        ]:
            with self.subTest(correct=correct, wrong=wrong):
                self.assertEqual(score(correct, wrong), expected)

    def test_input_counts(self):
        self.assertEqual(input_counts({
            "num_input_centroids": "2", "input_centroid_0_id": "123",
        }), (1, 1))

    def test_evaluation_distinguishes_failed_and_invalid_output(self):
        from tempfile import TemporaryDirectory
        from unittest.mock import patch
        from paired import evaluate

        case = {"id": "case", "condition": "clean", "seed": 1,
                "injected": 0, "attitude": [218, 12, 302]}
        with TemporaryDirectory() as directory:
            path = Path(directory)
            for process, expected in [
                ("failed", "failed"), ("crash", "crash"),
                ("timeout", "timeout"), ("ok", "invalid_output"),
            ]:
                with patch("paired.run", return_value=process):
                    row, inputs = evaluate(path / "lost", path, path,
                                           case, "py", path / "py.dat")
                self.assertEqual(row["status"], expected)
                self.assertEqual(row["outcome"], "")
                self.assertIsNone(inputs)
