# pulsar-lensing
Calculate gravitational lensing time delay on pulsar signals

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
