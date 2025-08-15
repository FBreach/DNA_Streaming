# Progressive Audio & Video Streaming from DNA Storage (Biased Raptor Codes)

Tools, configs, and paper artifacts for progressive retrieval of time‑dependent media from DNA using a biased Raptor fountain code. We skew degree‑1 probability toward early frames (and mildly boost I‑frames) so the peeling decoder resolves useful content sooner, enabling playback before full retrieval on a streaming nanopore channel.

** Video and Audio Demos will be uploaded shortly (issues with file sizes and type permissions)**

## Highlights
- **Frame‑aware degree‑1 bias** for low‑latency starts; optional **I‑frame multiplier** for video.
- **Streaming decode** with on‑the‑fly peeling (belief propagation).
- **Nanopore channel simulation** via Badread; configurable packets/s and error profiles.
- **Reproducible setup** (pip/conda friendly), CI via GitHub Actions, and a tidy project layout.

## Dependencies
- CMake
- Badread

### Install Badread
```bash
cd DNA_Streaming
git clone 'https://github.com/rrwick/Badread'
```

## Usage

### 1) Encode
```bash
cd NOREC4DNA
python encode.py data/Audio.bmp   --chunk_size 72   --error_correction reedsolomon   --overhead 0.5   --thr0 0.4   --a 0.01
```

### 2) Sequencing Simulation
```bash
badread simulate   --reference data/Audio_FASTA.fasta   --quantity 60x   --identity 97,99,1.0 | gzip > reads.fastq.gz
```

### 3) Decode
```bash
python decode_timed.py data/Audio_FASTA.ini --badread data/Audio_SIM.fasta
```

---

> **Note:** Paths and filenames above are illustrative and may need to be adapted to your local directory structure and outputs produced by the encode step.

