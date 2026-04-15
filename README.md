# RNA Probe Designer

A genetic-algorithm-based tool for designing DNA probe sequences optimised to adopt a user-specified secondary structure at a target melting temperature. Candidates are evaluated using [UNAFold](https://unafold.rna.albany.edu/) thermodynamic folding and scored against a configurable structural target. A Tkinter GUI lets you interactively build structure constraints, run the GA, and review ranked results.

---

## Features

- Interactive dot-bracket structure editor with per-base sequence constraints and mismatch weights
- Section-level sequence alternatives (e.g. full probe vs. –1 nt variant) evaluated in a single run
- Dynamic penalty spans that bias the GA toward keeping stems or terminal ends structurally stable
- Two-level warm-start pre-optimisation (fine then coarse segment passes) before the main GA
- IDT OligoAnalyzer API integration for accurate nearest-neighbour melting-temperature scoring
- Placeholder strand design for iMS/miRNA probe coupling
- Results exported to CSV; structures rendered with ViennaRNA (optional) or circular fallback layout

---

## Repository structure

```
gui.py                # Tkinter UI — program entry point
genetic_algorithm.py  # GA core: UNAFold I/O, scoring, evolution loop, GA() driver
idt_api.py            # IDT API auth, Tm lookups, placeholder-strand design
ensemble_helper.py    # Dot-bracket segmentation helpers used by the warm-start pass
dna_utils.py          # Shared nucleotide utilities (rnatoDNA, revcomp, DNAtorna)
requirements.txt      # Python package dependencies
.env.example          # Template for IDT API credentials
.gitignore
```

---

## Dependencies

### System tools (not installed by pip)

These must be installed separately before running the program.

| Tool | Version tested | Purpose |
|---|---|---|
| [UNAFold](https://unafold.rna.albany.edu/) | 4.0.1 | RNA/DNA folding — provides `hybrid-ss-min` |
| [Clustal Omega](http://www.clustal.org/omega/) | 1.2.2 | Multiple sequence alignment — provides `clustalo` |
| [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) *(optional)* | any | Improved structure layout in the GUI |

### Python packages

| Package | Purpose |
|---|---|
| `biopython` | Sequence alignment, file I/O |
| `matplotlib` | Segment distribution plotting |
| `numpy` | Numerical operations |
| `requests` | IDT API HTTP calls |

---

## Installation

### 1. Python

Python **3.11** (64-bit) is recommended. Download from [python.org](https://www.python.org/downloads/).

> **Windows:** on the installer's first screen, check **"Add Python to PATH"** before clicking Install Now.

### 2. UNAFold

1. Run the UNAFold `.msi` installer and accept the default install location (`C:\UNAfold`).
2. Confirm `C:\UNAfold\hybrid-ss-min.exe` exists after installation.

### 3. Clustal Omega

1. Download and extract the Clustal Omega win64 zip.
2. **Add the extracted folder to your system PATH** (recommended), or copy `clustalo.exe` to any folder already on your PATH (e.g. `C:\`).

To verify: open a terminal and run `clustalo --version`. You should see the version number, not an error.

### 4. Python packages

```bash
pip install -r requirements.txt
```

### 5. IDT API credentials

Melting temperatures are computed via the [IDT OligoAnalyzer REST API](https://www.idtdna.com/pages/tools/oligoanalyzer). A free IDT account with OAuth2 client credentials is required.

1. Copy `.env.example` to `.env` in the same folder.
2. Fill in your credentials:

```
IDT_USERNAME=your_username
IDT_PASSWORD=your_password
IDT_CLIENT_ID=your_client_id
IDT_CLIENT_SECRET=your_client_secret
```

The program reads these at runtime via environment variables. `.env` is already listed in `.gitignore` — never commit it.

> **Note:** an internet connection is required at runtime because every candidate's melting temperature is fetched from the IDT API during each GA generation.

---

## Running the program

```bash
python gui.py
```

Always launch from a terminal rather than double-clicking, so that any error messages remain visible.

---

## Usage

### 1. Enter and render a structure

Paste a valid dot-bracket string into the **Desired structure** box and click **Render**.

- Valid characters: `.` (unpaired), `(` and `)` (base-paired, must be balanced)
- The canvas draws the folded structure: gray backbone, blue base-pairs, numbered nodes
- If you see a validation error, your parentheses are unbalanced — fix and re-render

### 2. Per-base constraints *(optional)*

Click any node on the canvas to open its editor.

**Allowed bases** — controls which nucleotides the GA may place at this position:

| Setting | Meaning |
|---|---|
| Single base (e.g. `A`) | Hard-lock: only this base allowed |
| Subset (e.g. `A` + `G`) | Soft constraint: GA chooses among these |
| All four (`ACGT`) | Unconstrained (shown as `N`) |

**Positional mismatch weight (0–9)** — penalty applied when the GA produces the wrong structure at this position:

| Range | Meaning |
|---|---|
| `0` | Don't care |
| `1–3` | Gentle preference |
| `4–6` | Important region |
| `7–9` | Critical — do not deviate |

Click **Apply** to save. The node label updates to show the allowed set and the weight digit appears below the node.

### 3. Section alts *(optional)*

Use **Section alts** to define spans where the GA can evaluate whole-sequence alternatives in a single run — for example, a full probe alongside a –1 nt variant.

1. Right panel → **Section alts** → **Add section alt…**
2. Enter start and end indices (1-based), then pipe-separated alternatives, e.g.:
   - `TATTAACTG|TATTAA` (original + one alt)
   - `AAA|AA|A` (length sweep)
3. The span is exported as `[ORIGINAL|ALT1|ALT2]`, where `ORIGINAL` is whatever the per-base constraints currently specify for that region.
4. An orange badge marks the span on the canvas (drag it if it overlaps other labels — cosmetic only).

### 4. Gap groups / dynamic penalty spans *(optional)*

Gap groups add targeted structural penalties that encourage the GA to keep specific regions stable — useful for stems, terminal ends, or any region where secondary structure is critical.

1. **Shift+click** a start point:
   - A base node, or the `+` ghost node (5′ terminal — prevents extra bases being added at the start)
2. **Click** an end point:
   - A base node, or the `−` ghost node (3′ terminal — prevents extra bases being added at the end)
3. Enter a penalty (0–9) when prompted.
4. A coloured band appears above the span:
   - **Green** = includes a terminal (`+` or `−`)
   - **Peach** = internal span

Overlapping spans are merged automatically; the exporter keeps the maximum penalty and preserves terminal flags.

**Examples:**

| Goal | Span | Penalty |
|---|---|---|
| Stabilise 5′ end | `+` → base 5 | 3 |
| Protect 3′ end | last base → `−` | 9 |
| Emphasise a stem | base 20 → base 36 | 6 |

### 5. Export

Click **Export strings**. Four constraint strings are generated:

| Field | Description | Example |
|---|---|---|
| `desired_structure` | The dot-bracket string | `....((((...))))` |
| `seq_constraint` | Per-base sets + section-alt tokens | `NNN[ORIG\|ALT1]N[A\|G]NN` |
| `struct_constraints_base` | Mismatch weights, one digit per base | `000112222333300` |
| `dynamic_penalty` | Weights with penalty-span brackets overlaid | `[+000](3)112...9` |

A CSV snapshot of the exported design is also saved automatically to your Downloads folder.

> Re-export any time you change constraints before running the GA to keep these fields in sync. **Run Genetic Algorithm…** is disabled until you export at least once.

<p align="center">
  <img width="728" height="600" alt="image" src="https://github.com/user-attachments/assets/017aa076-694a-4a42-8b3e-a77b3e0065c6" />
</p>
### 6. Run the GA

Click **Run Genetic Algorithm…** and set parameters:

| Parameter | Description |
|---|---|
| Use melting temperature constraint | Disable to skip Tm scoring (faster; useful for structure-only optimisation) |
| Target melting temp (°C) | Desired probe Tm |
| Penalize off-target structures | Penalise candidates that fold into multiple suboptimal structures |
| Population size | Candidates per generation |
| Generations | Evolution steps per run |
| Runs | Independent GA runs; results are merged and ranked by fitness |
| Runs at a time | Parallel batch size (increase carefully — each run spawns a UNAFold subprocess) |
| [Na⁺] / [Mg²⁺] (M) | Salt concentrations passed to UNAFold |
| Generate placeholder strand | iMS/miRNA mode — requires a miRNA set via the iMS shortcut first |
| CSV output | File name and folder for saving results |

Click **Run**. A progress dialog opens; click **Hide** to dismiss it without stopping the computation, or **Cancel** to request an early stop (the current GA step finishes before cancellation takes effect).

### 7. Review results

A **Genetic Algorithm Results** window opens on completion:

- Navigate candidates with **◀ Prev** / **Next ▶** (sorted by fitness, best first)
- Top canvas: structure drawing (ViennaRNA layout if installed, circular fallback otherwise)
  - Blue lines = base pairs; gray lines = backbone; nodes show the nucleotide identity
- Details panel:
  - **Melting Temp (°C)** — IDT-computed Tm of the candidate
  - **Fitness** — GA score (lower = better)
  - **Length** — sequence length in nt
  - **Alt structures** — number of suboptimal structures predicted (excluding MFE)
  - **Placeholder / Placeholder Tm / miRNA Tm** — shown only in iMS/miRNA mode

All results are also written to the CSV file you specified before running.

---

## Configuration

The following environment variables can be set in your `.env` file or shell:

| Variable | Default | Description |
|---|---|---|
| `UNAFOLD_WORK_DIR` | `<system-temp>/unafold_runs` | Directory for UNAFold temporary run folders |
| `IDT_USERNAME` | — | IDT account username |
| `IDT_PASSWORD` | — | IDT account password |
| `IDT_CLIENT_ID` | — | IDT OAuth2 client ID |
| `IDT_CLIENT_SECRET` | — | IDT OAuth2 client secret |

`UNAFOLD_WORK_DIR` is created automatically if it does not exist. Override it if you want UNAFold's working files in a specific location (e.g. a fast local drive).

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `python` not recognised | Python not on PATH | Re-run the Python installer and check "Add Python to PATH" |
| `clustalo` not found | Clustal Omega not on PATH | Add its folder to PATH, or copy `clustalo.exe` to a folder that is already on PATH |
| `hybrid-ss-min` not found | UNAFold not installed | Run the UNAFold `.msi` installer |
| IDT authentication error | Missing or wrong credentials | Check your `.env` file against `.env.example` |
| GUI closes immediately | Unhandled runtime error | Launch from a terminal — the error message will be printed there |
| Structure renders as a circle | ViennaRNA not installed | Non-fatal; install ViennaRNA Python bindings for a better layout |
| `pip install` fails for ViennaRNA | Binary not available for your platform | ViennaRNA is optional — remove it from `requirements.txt` and re-run |
