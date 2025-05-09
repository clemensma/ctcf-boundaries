# ctcf-boundaries

**CTCF-boundaries** is a command-line tool that identifies CTCF-associated boundaries from Hi-C `.cool` files and CTCF ChIP-seq peak files in BED format. It uses [cooltools](https://github.com/open2c/cooltools) and [bioframe](https://github.com/open2c/bioframe) under the hood.

---

## 🚀 Installation

Clone the repo and install it:

```bash
git clone https://github.com/YOUR_USERNAME/ctcf-boundaries.git
cd ctcf-boundaries
pip install .
```

This installs the `ctcf-boundaries` command.

---

## 🔧 Usage

```bash
ctcf-boundaries COOL_FILE PEAKS_FILE OUTPUT_BED [--resolution RES] [--nproc N] [--genome GENOME]
```

### Example:

```bash
ctcf-boundaries data/sample.cool data/ctcf_peaks.bed ctcf_boundaries.bed --resolution 10000 --nproc 4
```

### Arguments:

| Argument       | Description                            |
|----------------|----------------------------------------|
| `COOL_FILE`    | Path to Hi-C `.cool` file              |
| `PEAKS_FILE`   | BED file with header row (CTCF peaks)  |
| `OUTPUT_BED`   | Output file for boundaries             |
| `--resolution` | Hi-C resolution (default: 10000)       |
| `--nproc`      | Number of threads (default: 4)         |
| `--genome`     | Genome assembly (default: hg38)        |

---

## 🧬 Requirements

- Python ≥ 3.8
- cooler
- cooltools
- bioframe
- pandas

Install dependencies using pip:

```bash
pip install cooler cooltools bioframe pandas
```

---

## 🛠 Features

- Calls insulation-based boundaries
- Filters to boundaries overlapping CTCF peaks
- BED-compatible output

---

## 📜 License

MIT License
