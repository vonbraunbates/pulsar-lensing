# pulsar-lensing
Calculate gravitational lensing time delay on pulsar signals

## Abstract

A key question in cosmology is the properties of dark matter. A particular open problem is
whether dark matter on small scales is clumpy, forming gravitationally-bound halos distributed
within the Galaxy. The practical difficulties inherent in testing this hypothesis stem from
the fact that, on astrophysical scales, dark matter is solely observable via its gravitational
interaction with other objects.

This thesis presents a gravitational-lensing-based solution for the mapping and characterisation of low-mass dark matter halos via their signature in millisecond pulsar observations.
This involves numerical calculations in three stages: first, determining the time delay and
magnification surfaces generated in the frame of reference of the halo; second, obtaining the
corresponding pulsar signature in the reference frame of the observer; and last, generalising
the method to multiple halos at varying distances. In both the single-lens and multiple-lens
cases, we discuss whether the delay is observationally detectable.

Dark matter halos act as gravitational lenses which produce a variable flux and induce ad-
ditional time delays in (tangent) bundles of photons passing near or through the halo. The key
dependency of the mass estimate is the density profile adopted for the halo. I utilise a variety
of proposed halo mass profiles — namely the elliptical model of Kochanek &al., the axially
symmetric Schwarzschild and homogenous disc lenses and the Navarro–Frenk–White
density profile — which are applicable over a broad range of halo masses. The pulsar simulations use the most realistic and sophisticated of these, the empirically-derived profile of Navarro, Frenk & White. 

I justify the adoption of a radially-symmetric density profile by showing that this greatly simplifies the calculation of the lens convergence. Moreover, I demonstrate that the use
of Hankel transforms is a novel way to increase the efficiency of the relativistic time delay.
The observational signatures of such halos are best identified using millisecond pulsars. This remarkable subset of the pulsar population has both the highest rotational frequencies and the most
period stability of all known pulsars. Furthermore, the potential for gravitational wave detection
using millisecond pulsars will result in an abundance of new data from pulsar surveys. I propose
that observational techniques do not require major adjustments when searching for signs of gravitational lensing, thus it is unnecessary to implement specialist data reduction pipelines, which enable the data from existing and future surveys to be examined for lensing with relative ease.

This thesis provides a practical method to search for dark matter halos within our Galaxy
and is readily extensible to nearby globular clusters and galaxies, pending the discovery of
millisecond pulsars in these more distant systems.

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
