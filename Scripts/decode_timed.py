#!/usr/bin/env python3
"""sequential_decode.py – adaptive stream decoder with progress logging

This script **random‑samples** oligos from a RU10 FASTA, feeds them
**one by one** into whatever variant of *RU10Decoder* (or wrapper demo) you have
installed, and records the packet index at which each original chunk becomes
available.

It actively **probes** the decoder object at run‑time, so it works whether your
installation exposes methods like `add_packet`, `add_oligo`, `ingest_read`, … or
something entirely different.

The result is a tidy two‑column CSV:

```
chunk_idx,first_packet
0,8
1,12
...
```

If nothing is decoded the script now prints a diagnostic dump of all methods it
tried—useful for tracking down yet‑another naming change.
"""
from __future__ import annotations

import argparse
import random
import csv
import time
import inspect
from pathlib import Path
from typing import Dict, List, Sequence, Callable, Any

from Bio import SeqIO  # pip install biopython

# NOREC4DNA ---------------------------------------------------------------
from norec4dna.RU10Decoder import RU10Decoder  # direct import still works for the base class

# ──────────────────────────────────────────────────────────────────────────────
# Dynamic method discovery helpers
# ──────────────────────────────────────────────────────────────────────────────

def _find_method(obj: Any, wanted_prefixes: Sequence[str], *, min_arity: int = 1) -> Callable:
    """Return the first callable whose name matches a prefix and arity.

    *prefixes*   – list/tuple of allowed name prefixes (e.g. ("add_", "push_"))
    *min_arity*  – minimum # positional parameters *excluding* ``self``
    """
    for name in dir(obj):
        if any(name.startswith(p) for p in wanted_prefixes):
            fn = getattr(obj, name)
            if callable(fn):
                sig = inspect.signature(fn)
                # Exclude *only* the implicit "self" parameter
                pos_args = [p for p in sig.parameters.values()
                            if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)]
                if len(pos_args) - 1 >= min_arity:
                    return fn
    raise AttributeError("No matching method found.")


# ──────────────────────────────────────────────────────────────────────────────
# Core routine
# ──────────────────────────────────────────────────────────────────────────────

def decode_stream(
    fasta_path: Path,
    config_path: Path,
    *,
    seed: int | None = None,
    sleep: float = 0.0,
) -> Dict[int, int]:
    """Return mapping `chunk_idx → first_packet_seen`."""

    packets = list(SeqIO.parse(str(fasta_path), "fasta"))
    if seed is not None:
        random.seed(seed)
    random.shuffle(packets)

    decoder = RU10Decoder(str(config_path))

    # --- find ingest method ---------------------------------------------------
    try:
        ingest = _find_method(decoder, (
            "add_packet", "add_oligo", "add_read", "input_new_packet", "input_", "push_", "ingest_", "feed_", "receive_", "process_"),
            min_arity=1,
        )
    except AttributeError as e:
        raise RuntimeError(
            "Could not locate a method that accepts incoming oligos. "
            "Update the prefix list in _find_method() to include your decoder's API." ) from e

    # --- find chunk‑query method ---------------------------------------------
    try:
        get_chunks = _find_method(decoder, (
            "get_decoded_chunk_indices", "get_currently_decoded_chunk_indices", "get_decoded_chunks", "decoded_chunks"),
            min_arity=0,
        )
    except AttributeError as e:
        # Fall back to attribute access – some builds expose a *set*
        if hasattr(decoder, "decoded_chunks"):
            get_chunks = lambda: list(getattr(decoder, "decoded_chunks"))  # type: ignore[assignment]
        else:
            raise RuntimeError(
                "Could not locate a method or attribute that returns decoded‑chunk indices." ) from e

    # --- find completion check (optional) ------------------------------------
    try:
        is_done = _find_method(decoder, ("is_decoding_complete", "isComplete", "is_decoded"), min_arity=0)
    except AttributeError:
        is_done = lambda: False  # keep feeding until packets exhausted

    first_seen: Dict[int, int] = {}
    prev_decoded: set[int] = set()

    for pkt_no, rec in enumerate(packets, 1):
        ingest(str(rec.seq))

        current_decoded = set(get_chunks())
        newly = current_decoded - prev_decoded
        for idx in newly:
            first_seen.setdefault(idx, pkt_no)
        prev_decoded = current_decoded

        if sleep > 0:
            time.sleep(sleep)

        if is_done():
            break

    return first_seen

# ──────────────────────────────────────────────────────────────────────────────
# CLI front‑end
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        prog="sequential_decode.py",
        description="Randomly sample RU10 oligos and record chunk recovery times",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    ap.add_argument("--fasta", required=True, type=Path, help="Input FASTA with oligos")
    ap.add_argument("--config", required=True, type=Path, help=".ini produced by encoder")
    ap.add_argument("--out", default="timeline.csv", type=Path, help="CSV output path")
    ap.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    ap.add_argument("--sleep", type=float, default=0.0, help="Seconds to wait between oligos")

    opt = ap.parse_args()

    try:
        mapping = decode_stream(opt.fasta, opt.config, seed=opt.seed, sleep=opt.sleep)
    except RuntimeError as err:
        # Emit the decoder's public API to help the user adapt the prefix lists
        import pprint, textwrap
        decoder = RU10Decoder(str(opt.config))
        methods = [m for m in dir(decoder) if callable(getattr(decoder, m)) and not m.startswith("__")]
        print("\n[!] " + str(err))
        print("\nDecoder exposes the following callables:")
        print(textwrap.fill(", ".join(methods), width=100, subsequent_indent="    "))
        print("\nEdit *sequential_decode.py* and add the correct prefix to _find_method().")
        return

    with opt.out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["chunk_idx", "first_packet"])
        for idx in sorted(mapping):
            writer.writerow([idx, mapping[idx]])

    recovered = len(mapping)
    horizon = max(mapping.values()) if mapping else 0
    print(f"✔ Recovered {recovered} chunks after {horizon} packets → {opt.out}")


if __name__ == "__main__":
    main()

