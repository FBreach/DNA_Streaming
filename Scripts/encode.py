#!/usr/bin/env python3
"""encode_raw.py – **modified version with bias parameters**

This script wraps the RU10Encoder and exposes additional parameters for
temporal biasing.  In addition to prioritising the first N chunks,
the command line now accepts two parameters, ``--a`` and ``--thr_0``,
which control the weighting function used by the biased Raptor
distribution.  These parameters are not interpreted by this script
directly; instead they are propagated to the decoder via the
configuration file produced at the end of the run.  Downstream
decoders can use the encoded values to reproduce the same degree
distribution.

Example
-------
```
python encode_raw.py --path data/bitstream.bin --chunk_size 280 \
                     --insert_header --priority_first 10 \
                     --a 0.5 --thr_0 0.12
```
This pushes chunks 0–9 to the front of the emit queue, encodes using
the default chunk size and bias parameters, and writes
``bitstream_RU10.fasta`` next to the input file.  A companion
``encoding_config.ini`` stores all encoding options, including ``a``
and ``thr_0``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from norec4dna.Encoder import Encoder
from norec4dna.RU10Encoder import RU10Encoder
from norec4dna.ErrorCorrection import nocode, get_error_correction_encode
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.distributions.RaptorDistribution import RaptorDistribution

# ──────────────────────────────────────────────────────────────────────────────
# Constants mirroring RU10Encoder defaults
# ──────────────────────────────────────────────────────────────────────────────
ID_LEN_FORMAT = "H"                 # uint16
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"   # uint32
CRC_LEN_FORMAT = "I"               # uint32
PACKET_LEN_FORMAT = "I"             # uint32
DEFAULT_CHUNK_SIZE = 71             # keeps oligos ≤ 150 nt by default

# ──────────────────────────────────────────────────────────────────────────────


def _encode(
    file_path: Path,
    *,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    error_correction = nocode,
    insert_header: bool = False,
    save_number_of_chunks_in_packet: bool = False,
    mode_1_bmp: bool = False,
    prepend: str = "",
    append: str = "",
    upper_bound: float = 0.5,
    overhead: float = 0.40,
    checksum_len_str: str | None = None,
    xor_by_seed: bool = False,
    mask_id: bool = True,
    id_spacing: int = 0,
    priority_first: int = 0,
    p_thr: float = 0.0,
    a: float | None = None,
    thr_0: float | None = None,
) -> tuple[RU10Encoder, Path]:
    """Run RU10Encoder and return ``(encoder, fasta_path)``.

    Parameters
    ----------
    file_path : Path
        Path to the binary file to encode.
    chunk_size : int, optional
        Payload size in bytes per packet; defaults to ``DEFAULT_CHUNK_SIZE``.
    error_correction : callable, optional
        Error-correction function returned by ``get_error_correction_encode``.
    insert_header : bool, optional
        If true, insert a special header oligo.
    save_number_of_chunks_in_packet : bool, optional
        If true, embed the total chunk count in every oligo.
    upper_bound : float, optional
        Drop threshold for noisy packets.
    overhead : float, optional
        Desired redundancy overhead.
    checksum_len_str : str | None, optional
        Struct format string for header CRC length.
    xor_by_seed : bool, optional
        If true, XOR the payload by the seed before synthesis.
    mask_id : bool, optional
        Apply Gray-code masking to packet identifiers.
    id_spacing : int, optional
        Base spacing for masked identifiers when masking is on.
    priority_first : int, optional
        Give priority to the first N chunks (0 ⇒ disabled).
    p_thr : float, optional
        Probability threshold for the Raptor distribution.
    a : float | None, optional
        Bias parameter controlling the maximum singleton rate.  It is
        recorded in the configuration file so that decoders can use the
        same value.  If ``None``, the default decoder setting is used.
    thr_0 : float | None, optional
        Bias parameter controlling the rate decay across frames.  It is
        recorded in the configuration file so that decoders can use the
        same value.  If ``None``, the default decoder setting is used.
    """

    # determine number of chunks for the given chunk size
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(
        file_path, chunk_size
    )

    # Instantiate a Raptor distribution.  The distribution itself does not
    # interpret the bias parameters directly; they will be passed through
    # to the decoder via the configuration file.
    dist = RaptorDistribution(number_of_chunks, p_thr=p_thr)
    rules = FastDNARules()

    # ⚠  RU10Encoder internally does str-concat on its ``file`` attribute,
    #     so we hand it a plain string instead of a Path.
    enc = RU10Encoder(
        str(file_path),
        number_of_chunks,
        dist,
        insert_header=insert_header,
        pseudo_decoder=None,
        chunk_size=0,  # let encoder derive from number_of_chunks
        rules=rules,
        error_correction=error_correction,
        packet_len_format=PACKET_LEN_FORMAT,
        crc_len_format=CRC_LEN_FORMAT,
        number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT,
        id_len_format=ID_LEN_FORMAT,
        save_number_of_chunks_in_packet=save_number_of_chunks_in_packet,
        mode_1_bmp=mode_1_bmp,
        prepend=prepend,
        append=append,
        drop_upper_bound=upper_bound,
        checksum_len_str=checksum_len_str,
        xor_by_seed=xor_by_seed,
        mask_id=mask_id,
        id_spacing=id_spacing,
    )

    enc.set_overhead_limit(overhead)

    if priority_first > 0:
        enc.set_priority_chunks(priority_first)

    # Encode and write FASTA
    enc.encode_to_packets()
    fasta_path = enc.save_packets_fasta(file_ending="_RU10", seed_is_filename=True)
    return enc, Path(fasta_path)


# ──────────────────────────────────────────────────────────────────────────────
# CLI front-end
# ──────────────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="encode_raw.py",
        description="Encode a raw binary stream into RU10 oligo packets (no JPEG processing)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--path", required=True, help="Binary file to encode")
    parser.add_argument("--chunk_size", type=int, default=DEFAULT_CHUNK_SIZE,
                        help="Chunk size in bytes")
    parser.add_argument("--insert_header", action="store_true",
                        help="Insert special RU10 header oligo")
    parser.add_argument("--save_number_of_chunks", action="store_true",
                        help="Embed total chunk count in every oligo")
    parser.add_argument("--error_correction", default="nocode",
                        help="EC scheme (nocode, rs, ldpc, etc.)")
    parser.add_argument("--repair_symbols", type=int, default=20,
                        help="Repair symbols for chosen EC scheme")
    parser.add_argument("--drop_upper_bound", type=float, default=0.5,
                        help="Dropping threshold for noisy packets")
    parser.add_argument("--overhead", type=float, default=0.4,
                        help="Desired redundancy overhead")
    parser.add_argument("--header_crc_str", default="B",
                        help="struct.format string for header CRC length")
    parser.add_argument("--xor_by_seed", action="store_true",
                        help="XOR payload by seed prior to synthesis")
    parser.add_argument("--no_mask_id", action="store_true",
                        help="Disable Gray-code ID masking")
    parser.add_argument("--id_spacing", type=int, default=0,
                        help="Base spacing for masked IDs when masking is on")
    parser.add_argument("--p_thr", type=float, default=0.0,
                        help="Probability threshold for Raptor distribution")

    # prioritisation – only first-N is actually supported by RU10Encoder
    parser.add_argument("--priority_first", type=int, default=0, metavar="N",
                        help="Give priority to the first N chunks (0 ⇒ disabled)")

    # bias parameters: these are recorded and passed to the decoder via the config
    parser.add_argument("--a", type=float, default=None,
                        help="Bias parameter controlling the maximum singleton rate")
    parser.add_argument("--thr_0", type=float, default=None,
                        help="Bias parameter controlling the rate decay across frames")

    opt = parser.parse_args()

    in_path = Path(opt.path).expanduser().resolve()
    if not in_path.exists():
        raise FileNotFoundError(in_path)

    ec_fn = get_error_correction_encode(opt.error_correction, opt.repair_symbols)

    enc, fasta_out = _encode(
        in_path,
        chunk_size=opt.chunk_size,
        error_correction=ec_fn,
        insert_header=opt.insert_header,
        save_number_of_chunks_in_packet=opt.save_number_of_chunks,
        upper_bound=opt.drop_upper_bound,
        overhead=opt.overhead,
        checksum_len_str=opt.header_crc_str,
        xor_by_seed=opt.xor_by_seed,
        mask_id=not opt.no_mask_id,
        id_spacing=opt.id_spacing,
        p_thr=opt.p_thr,
        priority_first=opt.priority_first,
        a=opt.a,
        thr_0=opt.thr_0,
    )

    # Write a configuration file capturing all parameters used during encoding.
    import configparser
    _cfg = configparser.ConfigParser()
    _cfg['EncodingSettings'] = {
        'input': str(in_path),
        'out_dir': str(fasta_out.parent),
        'chunk_size': str(opt.chunk_size),
        'error_correction': opt.error_correction,
        'repair_symbols': str(opt.repair_symbols),
        'overhead': str(opt.overhead),
        'xor_by_seed': str(opt.xor_by_seed),
        'mask_id': str(not opt.no_mask_id),
        'id_spacing': str(opt.id_spacing),
        'p_thr': str(opt.p_thr),
        'priority_chunks': str(opt.priority_first),
        # record bias parameters even if None so decoders know whether defaults apply
        'a': str(opt.a) if opt.a is not None else '',
        'thr_0': str(opt.thr_0) if opt.thr_0 is not None else '',
    }
    _cfg['Files'] = {
        'fasta': str(fasta_out),
        'config': str(fasta_out.parent / 'encoding_config.ini'),
    }
    cfg_path = fasta_out.parent / 'encoding_config.ini'
    with open(cfg_path, 'w') as _f:
        _cfg.write(_f)

    print(f"✔ Encoding complete. FASTA saved to: {fasta_out}")


if __name__ == "__main__":
    main()