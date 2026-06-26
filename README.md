# pulsar-lensing
Calculate gravitational lensing time delay on pulsar signals

## Abstract

Dark matter may form low-mass gravitationally bound halos within the Galaxy, but detecting such structures is difficult because dark matter interacts primarily through gravity. This thesis develops a method for detecting and characterising dark matter halos through their gravitational lensing signatures in millisecond pulsar observations.

A numerical framework is constructed to calculate lensing-induced magnification and relativistic time delays for both single and multiple halo systems. The method is applied to a range of halo density profiles, including the Navarro-Frenk-White model, and incorporates a novel use of Hankel transforms to improve the efficiency of time-delay calculations. Simulations are used to determine the observational signatures expected from lensing events and to assess their detectability.

The results show that low-mass dark matter halos can produce measurable timing and flux variations in millisecond pulsars under realistic conditions. Because pulsar timing observations already achieve exceptional precision, the proposed technique can be applied to existing and future pulsar surveys without substantial changes to data-processing pipelines. This provides a practical observational method for mapping dark matter substructure in the Galaxy and, potentially, in nearby globular clusters and galaxies.

## Compilation

Load the files into MATLAB Online.  No Makefile nor compilation is needed.

## Licence

The thesis text, figures, tables and accompanying documentation are licensed under CC BY 4.0.  For the full text (which is reporduced in `./LICENCE-THESIS`) see https://creativecommons.org/licenses/by/4.0/

The source code is licenced under the MIT Licence, as in `./LICENCE-CODE`.

## Citation

Please cite the following thesis.  If you are using `biblatex` use the `thesis` class:

```bibtex
@thesis{VonBraunBates2014,
  author       = {Francesca von Braun-Bates},
  title        = {Gravitational Lensing of Pulsars as a Probe of Dark Matter Halos},
  type         = {Master's thesis},
  institution  = {University of Sydney},
  date         = {2014-04},
  eprint       = {2210.06151},
  eprinttype   = {arXiv},
  eprintclass  = {astro-ph.CO},
  url          = {https://arxiv.org/abs/2210.06151},
}```
```
