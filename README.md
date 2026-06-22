# pulsar-lensing
Calculate gravitational lensing time delay on pulsar signals

## Abstract

Dark matter may form low-mass gravitationally bound halos within the Galaxy, but detecting such structures is difficult because dark matter interacts primarily through gravity. This thesis develops a method for detecting and characterising dark matter halos through their gravitational lensing signatures in millisecond pulsar observations.

A numerical framework is constructed to calculate lensing-induced magnification and relativistic time delays for both single and multiple halo systems. The method is applied to a range of halo density profiles, including the Navarro-Frenk-White model, and incorporates a novel use of Hankel transforms to improve the efficiency of time-delay calculations. Simulations are used to determine the observational signatures expected from lensing events and to assess their detectability.

The results show that low-mass dark matter halos can produce measurable timing and flux variations in millisecond pulsars under realistic conditions. Because pulsar timing observations already achieve exceptional precision, the proposed technique can be applied to existing and future pulsar surveys without substantial changes to data-processing pipelines. This provides a practical observational method for mapping dark matter substructure in the Galaxy and, potentially, in nearby globular clusters and galaxies.

## Compilation

With a valid TeX installation:

```sh
    latexmk -pdf main.tex
```

You will need the following packages:

- amssymb
- booktabs
- cleveref
- graphicx
- hyperref
- lscape
- multirow
- subcaption
- threeparttable
- xparse

The thesis useds a custom class file based on `amsbook`: if you have difficulties, start with `amsbook` instead.

## Licence

The thesis text, figures, tables and accompanying documentation are licensed under CC BY 4.0.  For the full text (which is reporduced in `./LICENCE-THESIS`) see https://creativecommons.org/licenses/by/4.0/

The source code is licenced under the MIT Licence, as in `./LICENCE-CODE`.

## Citation

Please cite the following thesis.  If you are using `biblatex` use the `thesis` class:

```bibtex
@thesis{VonBraunBates2010,
  author      = {{von Braun-Bates}, Francesca},
  title       = {Gravitational lensing of pulsars as a probe of dark matter halos},
  type        = {Masters Thesis},
  institution = {University of Sydney},
  year        = {2014},
  month       = apr,
}
```
Otherwise use `misc`:

```bibtex
@misc{VonBraunBates2010,
  author = {{von Braun-Bates}, Francesca},
  title  = {Gravitational lensing of pulsars as a probe of dark matter halos,
  year   = {2014},
  note   = {Masters Thesis, University of Sydney}
}
```
