Progressive Audio & Video Streaming from DNA Storage (Biased Raptor Codes)
Tools, configs, and paper artifacts for progressive retrieval of time-dependent media from DNA using a biased Raptor fountain code. We skew degree-1 probability toward early frames (and mildly boost I-frames) so the peeling decoder resolves useful content sooner, enabling playback before full retrieval on a streaming nanopore channel.

Paper PDF: paper/Breach_2025_MRes_Thesis.pdf (see paper/)

Highlights
Frame-aware degree-1 bias for low-latency starts; optional I-frame multiplier for video

Streaming decode with on-the-fly peeling (belief propagation)

Nanopore channel simulation via Badread; configurable packets/s and error profiles

Reproducible Python env (pip/conda), GitHub Actions CI, and a tidy project layout

Repo structure
bash
Copy
Edit
src/raptor_bias/        # Package stubs to extend
  ├─ encoder.py         # biased Raptor encoder (implement here)
  ├─ decoder.py         # peeling decoder (implement here)
  └─ simulate.py        # CLI entry for streaming simulation
scripts/
  └─ run_simulation.py  # Thin wrapper around simulate.py
config/
  └─ experiment.yaml    # Example bias + channel settings
data/{raw,processed}/    # Your inputs/outputs (gitkept)
figures/                 # Plots & JSON artifacts (gitkept)
paper/                   # Thesis / paper PDF
tests/                   # Minimal sanity tests
Installation
Conda (recommended)

bash
Copy
Edit
mamba env create -f environment.yml
mamba activate dna-streaming
pip install -e .
Pip (virtualenv)

bash
Copy
Edit
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .
Quick start
Run a tiny end-to-end simulation stub and write a demo artifact:

bash
Copy
Edit
python scripts/run_simulation.py --config config/experiment.yaml --out figures/demo.json
Smoke tests:

bash
Copy
Edit
pytest -q
python -m raptor_bias.simulate --help
Configuration
Edit config/experiment.yaml:

yaml
Copy
Edit
packets_per_second: 200
bias:
  theta0: 0.12        # max degree-1 rate at frame 0
  alpha: 0.5          # decay across frames
  iframe_multiplier: 1.2
  theta_max: 0.5
Swap in your own bias schedule, error model, or frame weighting. The stub prints a JSON summary; extend it to emit metrics/plots used in the paper.

Usage patterns
Audio: encode chunks by time; bias early frames to minimize startup.

Video: add a modest I-frame boost so dependent frames reconstruct smoothly.

Channel: simulate nanopore reads with Badread (installed via requirements) to approximate ONT-like indel errors and streaming arrival.

Reproducing figures
Place your plotting code under scripts/ or extend simulate.py to emit per-frame recovery curves, startup latency, and overhead. Save outputs in figures/ and commit small thumbnails if desired.

Contributing
PRs welcome—especially for:

A reference biased-Raptor encoder/decoder

Plotting scripts to reproduce paper figures

Adapters for real ONT runs (Fast5/FastQ → packet stream)

Citation
Please cite the paper and this repository (see CITATION.cff).

License
MIT (see LICENSE).

Related & credits
RaptorPJPEG — progressive JPEG over DNA with on-the-fly error handling; this README’s structure borrows from that project’s clear, pipeline-oriented docs. 
GitHub

