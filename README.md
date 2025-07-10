# Proving Stability of the Active Flux Method in the Framework of Summation-by-Parts Operators

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/TODO.svg)](https://zenodo.org/doi/TODO)

This repository contains information and code to reproduce the results presented
in the article
```bibtex
@online{barsukow2025stability,
  title={Stability of the Active Flux Method in the
         Framework of Summation-by-Parts Operators},
  author={Barsukow, Wasilij and Klingenberg, Christian and
          Lechner, Lisa and Nordstr{\"o}m, Jan and Ortleb, Sigrun
          and Ranocha, Hendrik},
  year={2025},
  month={07},
  eprint={TODO},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{barsukow2025provingRepro,
  title={Reproducibility repository for
         "{S}tability of the Active Flux Method in the
         Framework of Summation-by-Parts Operators"},
  author={Barsukow, Wasilij and Klingenberg, Christian and
          Lechner, Lisa and Nordstr{\"o}m, Jan and Ortleb, Sigrun
          and Ranocha, Hendrik},
  year={2025},
  howpublished={\url{https://github.com/ranocha/2025_active_flux_sbp}},
  doi={TODO}
}
```

## Abstract

The Active Flux method is a numerical method for conservation laws using a combination of cell averages and point values, based on ideas from finite volumes and finite differences. This seemingly unusual mix has been shown to work well in many situations. We expand the theoretical justifications of the Active Flux method by analyzing it from the point of view of summation-by-parts (SBP) operators, which have been used routinely to analyze finite difference, finite-volume and discontinuous Galerkin schemes. We demonstrate that the Active Flux method can be formulated using degenerate SBP operators, yielding a first and novel approach to show energy stability of the Active Flux method.


## Numerical experiments

To reproduce the numerical experiments presented in this article, you need
to install [Julia](https://julialang.org/).
The numerical experiments presented in this article were performed using
Julia v1.10.9.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface. Then, you need to start
Julia in the `code` directory of this repository and follow the instructions
described in the `README.md` file therein.


## Authors

- Wasilij Barsukow
- Christian Klingenberg
- Lisa Lechner
- Jan Nordström
- Sigrun Ortleb
- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
