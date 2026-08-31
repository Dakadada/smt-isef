import unittest
from pathlib import Path

from rawread_reanalysis.common import read_tsv


ROOT = Path(__file__).resolve().parents[1]
EVENTS = ROOT / "config/event_junctions.GRCh38.tsv"


class EventJunctionDefinitionsTest(unittest.TestCase):
    EXONS = {
        "SMN2_exon7_inclusion": {
            "strand": "+",
            "exons": ((70070641, 70070751), (70076521, 70076574), (70077019, 70077595)),
            "ends": ("TGGCTATCATACTGGCTATTATATG", "GTGCTCACATTCCTTAAATTAAGGA"),
            "starts": ("GGTTTTAGACAAAATCAAAAAGAAG", "GAAATGCTGGCATAGAGCAGCACTA"),
        },
        "DMD_exon53_inclusion": {
            "strand": "-",
            "exons": ((31729631, 31729748), (31679375, 31679586), (31657990, 31658144)),
            "ends": ("GCTAGAACAATCATTACGGATCGAA", "CCAAAAGAAAATCACAGAAACCAAG"),
            "starts": ("TTGAAAGAATTCAGAATCAGTGGGA", "CAGTTGGCCAAAGACCTCCGCCAGT"),
        },
    }

    @staticmethod
    def _bed_boundaries(upstream, downstream, anchor_length=20):
        # quantify.py compares BED12 chromStart/chromEnd, which include the anchors.
        # BED coordinates remain ascending and 0-based even when transcript order is reverse genomic.
        left, right = sorted((upstream, downstream), key=lambda exon: exon[0])
        return left[1] - anchor_length, right[0] - 1 + anchor_length

    def test_markers_span_the_expected_junctions_and_are_mutually_exclusive(self):
        rows = {row["event"]: row for row in read_tsv(EVENTS)}
        self.assertEqual(set(rows), set(self.EXONS))
        for event, definition in self.EXONS.items():
            row = rows[event]
            exon1, exon2, exon3 = definition["exons"]
            end1, end2 = (sequence[-20:] for sequence in definition["ends"])
            start2, start3 = (sequence[:20] for sequence in definition["starts"])
            expected_markers = {
                "inclusion_up": end1 + start2,
                "inclusion_down": end2 + start3,
                "skipping": end1 + start3,
            }
            markers = {kind: row[f"{kind}_marker"] for kind in expected_markers}
            self.assertEqual(markers, expected_markers)
            for marker in markers.values():
                self.assertEqual(len(marker), 40)
            for kind, marker in markers.items():
                self.assertFalse(any(marker in other for other_kind, other in markers.items() if other_kind != kind))

            expected_coordinates = {
                "inclusion_up": self._bed_boundaries(exon1, exon2),
                "inclusion_down": self._bed_boundaries(exon2, exon3),
                "skipping": self._bed_boundaries(exon1, exon3),
            }
            self.assertEqual(row["strand"], definition["strand"])
            for kind, (start, end) in expected_coordinates.items():
                self.assertEqual((row[f"{kind}_start"], row[f"{kind}_end"]), (str(start), str(end)))


if __name__ == "__main__":
    unittest.main()
