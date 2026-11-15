from conference.aso_pipeline.sequences import ASOS, ISS_N1_CONTEXT, derive_target


def test_targets_are_reverse_complementary_unique_context_windows():
    expected = {
        "F18": (10, 27, "CCAGCAUUAUGAAAGUGA"),
        "F14": (10, 23, "CCAGCAUUAUGAAA"),
        "F10": (10, 19, "CCAGCAUUAU"),
        "3UP8": (7, 14, "CUGCCAGC"),
    }
    for aso_id, (start, end, target) in expected.items():
        window = derive_target(ASOS[aso_id].rna_sequence, ISS_N1_CONTEXT)
        assert (window.start_1based, window.end_1based, window.sequence) == (
            start,
            end,
            target,
        )


def test_3up8_target_starts_three_nt_upstream_of_f18():
    upstream = derive_target(ASOS["3UP8"].rna_sequence, ISS_N1_CONTEXT)
    core = derive_target(ASOS["F18"].rna_sequence, ISS_N1_CONTEXT)
    assert core.start_1based - upstream.start_1based == 3
    assert upstream.requires_human_confirmation is True

