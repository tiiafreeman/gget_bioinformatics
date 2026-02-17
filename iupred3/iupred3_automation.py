import os
import iupred3_lib


def run_iupred3_on_sequence(sequence_name, sequence, iupred_type, smoothing, anchor2=False):

    if iupred_type not in {"long", "short", "glob"}:
        raise ValueError("iupred_type must be one of: 'long', 'short', 'glob'.")

    if smoothing not in {"none", "medium", "strong"}:
        raise ValueError("smoothing must be one of: 'none', 'medium', 'strong'.")

    seq = "".join(str(sequence).split()).upper()
    if not seq:
        raise ValueError("Input sequence is empty.")

    scores, glob_text = iupred3_lib.iupred(seq, mode=iupred_type, smoothing=smoothing)

    if len(scores) != len(seq):
        raise ValueError(f"Score length ({len(scores)}) != sequence length ({len(seq)})")

    output = {
        "name": sequence_name,
        "sequence": seq,
        "scores": [float(s) for s in scores],
        "glob_text": glob_text
    }

    if anchor2:
        a2 = iupred3_lib.anchor2(seq)
        if len(a2) != len(seq):
            raise ValueError(f"Anchor2 score length ({len(a2)}) != sequence length ({len(seq)})")
        output["anchor2"] = [float(s) for s in a2]

    return output


result = run_iupred3_on_sequence("TP53", "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP", "long", "medium", "True")
# run_iupred3_on_sequence(sequence_name, sequence, iupred_type, smoothing)

print(result["name"], len(result["scores"]), result["scores"], result["anchor2"] if "anchor2" in result else "No anchor2 data")