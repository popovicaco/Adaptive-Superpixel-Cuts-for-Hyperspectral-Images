# Adaptive Superpixel Cuts for Hyperspectral Images

An unsupervised, graph-based segmentation algorithm for hyperspectral images (HSI). It combines SLIC superpixels, normalized cuts, and graph-regularized spectral unmixing so that materials with similar spectra but different spatial locations can still be separated, at a fraction of the cost of pixel-level clustering.

This repository holds the Python implementation (`SuperpixelCutsPy`), the experiment notebooks, and the LaTeX source of the accompanying thesis by Aleksandar Popovic, supervised by Dr. You Liang and Dr. Na Yu.

![Algorithm overview](thesis%20writing/algorithm_view.png)

## How it works

1. **Preprocessing** – optional SVD/PCA denoising, then per-band normalization of the cube.
2. **Superpixel creation** – a hyperspectral SLIC variant groups pixels into `n_s` compact, spectrally homogeneous superpixels; each is represented by its mean spectrum.
3. **Initial segmentation** – a weighted graph is built over superpixels, with edge weights combining spectral similarity (spectral angle, scaled by `σ`) and a spatial cutoff (`κ`). Normalized cuts splits it into `n_e` segments.
4. **Graph-regularized unmixing** – the mean spectrum of each initial segment is treated as an endmember, and abundances are estimated for every superpixel with an ADMM solver (fully constrained least squares plus a graph penalty `β` that keeps nearby superpixels' abundances similar).
5. **Final segmentation** – the abundances are appended to each superpixel's spectrum to form a new feature vector, the graph is rebuilt, and normalized cuts produces the final label map.

## Repository layout

```
SuperpixelCutsPy/                 Core library
  Preprocesser.py                 Normalization, Gaussian blur, SVD and PCA denoising
  superpixel.py                   Hyperspectral SLIC (assignments + superpixel library)
  normalized_cuts.py              Affinity matrices, normalized cuts, full adaptive pipeline
  unmixing.py                     NNLS / FCLS and (graph-regularized) ADMM unmixing
  segmentation_evaluation.py      Label matching, confusion matrices, noise utilities
  utility.py                      Cube <-> matrix reshaping, image export helpers

graph_adaptivencuts_*.ipynb       Experiments: Salinas A/B, Samson A/B/C, BHSI and others
result_compilation*.ipynb         Metrics and figures for the thesis results
explore_*.ipynb                   Exploration of the full Salinas and Samson scenes
adaptive_ncuts/                   Earlier iteration of the method (notebooks + jupytext .md exports)
thesis writing/                   LaTeX source, figures and compiled PDF of the thesis
```

## Installation

The code was developed with Python 3.11. There is no `setup.py`; clone the repo and work from its root so `SuperpixelCutsPy` is importable.

```bash
git clone https://github.com/popovicaco/Adaptive-Superpixel-Cuts-for-Hyperspectral-Images.git
cd Adaptive-Superpixel-Cuts-for-Hyperspectral-Images
pip install numpy scipy scikit-learn matplotlib h5py
```

## Quick start

```python
import numpy as np
import scipy as sp
from SuperpixelCutsPy import Preprocesser, superpixel, normalized_cuts

# Load a hyperspectral cube of shape (nx, ny, n_bands)
cube = sp.io.loadmat("data/Salinas_corrected.mat")["salinas_corrected"][35:90, 25:175, :]

# 1. Preprocess
pp = Preprocesser.Preprocesser(data=cube)
pp.singular_value_decomposition(n_svd=5)
pp.layer_normalization()
X = pp.data.copy()

# 2. Superpixels
assignments, centers = superpixel.generate_SLIC_assignments(data=X, n_superpixels=250, slic_m_param=2)
superpixeled_cube, superpixel_library = superpixel.generate_SLIC_superpixels(data=X, assignments=assignments)

# 3–5. Adaptive superpixel cuts
labels, segment_spectra, intermediate = normalized_cuts.graph_regularized_ncuts_admm(
    data=X,
    superpixel_library=superpixel_library,
    superpixel_centers=centers,
    superpixel_assignments=assignments,
    n_endmembers=4,              # n_e: number of segments
    spectral_sigma2_param=0.05,  # σ: spectral similarity scale
    spatial_kappa_param=10,      # κ: spatial neighbourhood cutoff
    spatial_dmax_param=10,       # distance cutoff for the unmixing graph
    spatial_beta_param=0.025,    # β: graph regularization strength
    n_unmixing_iters=200,
    spectral_metric="SAM",
)

# labels: (nx, ny) segmentation map
# intermediate: initial_labels, initial_spectra, abundance_results, unmixing_history
```

The `graph_adaptivencuts_test_*.ipynb` notebooks walk through this pipeline end to end for each dataset, including plots and evaluation against ground truth.

## Choosing parameters

| Parameter | Meaning | Guidance |
|---|---|---|
| `n_superpixels` (`n_s`) | Number of superpixels | Controls runtime. About 16 pixels per superpixel (`n_pixels / n_s ≈ 16`) works well without much tuning of `m`. |
| `slic_m_param` (`m`) | SLIC compactness | 2–3 was used in all experiments. |
| `n_endmembers` (`n_e`) | Number of segments | Set from prior knowledge of the scene. |
| `spectral_sigma2_param` (`σ`) | Spectral similarity scale | Main tuning knob, searched over 0.001–0.1. |
| `spatial_kappa_param` (`κ`) | Spatial cutoff between superpixels | Main tuning knob, searched over 10–40. |
| `unmixing_mu_param` (`μ`) | ADMM penalty | Results are largely insensitive to it; default 1. |
| `spatial_beta_param` (`β`) | Graph regularization of abundances | Small values (≈0.0025–0.025) help separate spectrally similar but spatially distinct regions. |

## Datasets

Data files are not included (`data/` is git-ignored). The experiments use two standard remote-sensing scenes, saved as `.mat` files under `data/`:

- **Salinas** – AVIRIS, 512 × 217 pixels, 204 bands after removing noisy channels, 16 ground-truth classes. Subsets used: Salinas A (86 × 83, 6 classes) and Salinas B (55 × 150 from pixel (35, 25), 5 classes). Available from the [Hyperspectral Remote Sensing Scenes](https://www.ehu.eus/ccwintco/index.php/Hyperspectral_Remote_Sensing_Scenes) collection.
- **Samson** – 952 × 952 pixels, 156 bands (400–900 nm), no ground truth. Three 95 × 95 shoreline subsets (A, B, C) containing water, soil and trees.

The `*_bhsi*` notebooks apply the method to biological hyperspectral images, the application that motivated this work in collaboration with Dr. Yeni Yucel's Ophthalmic Pathology lab at St. Michael's Hospital.

## Results

| Dataset | `n_s` | `σ` | `κ` | Overall accuracy |
|---|---|---|---|---|
| Salinas A | 420 | 0.0025 | 20 | 0.993 |
| Salinas B | 330 | 0.0025 | 10 | 0.998 |
| Samson A/B/C | 961 | 0.015 | 30 | No ground truth; assessed visually |

Per-class accuracy and IoU tables, along with segmentation maps, are in the thesis PDF.

![Samson results](thesis%20writing/samsonabc-results.png)

## Thesis

The full write-up, covering background on hyperspectral unmixing, SLIC, spectral clustering and ADMM, the algorithm derivation, and the experimental results, is in [`thesis writing/Adaptive Superpixel Cuts.pdf`](thesis%20writing/Adaptive%20Superpixel%20Cuts.pdf).

## Acknowledgements

Thanks to Dr. You Liang and Dr. Na Yu for supervising this work, and to Dr. Yeni Yucel and the Ophthalmic Pathology lab at St. Michael's Hospital for the motivating application.
