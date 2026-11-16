# pyEnvipat

Python port of the R package [enviPat](https://github.com/blosloos/enviPat)
(Loos, Gerber et al., Eawag): isotope fine structure patterns, theoretical mass
spectrometry envelopes, centroid detection and instrument resolution
interpolation.

The distribution is called `pyEnvipat`; the import name keeps the original
package name:

```python
import envipat as ep
import envipat.data as data
```

The algorithms run in a compiled C extension built from the original,
unmodified enviPat C sources. Results are checked against reference outputs
frozen from enviPat v2.8.

```{toctree}
:maxdepth: 2

installation
quickstart
tutorial
api/index
data
```

## Citation

If you use pyEnvipat in your work, please cite the original algorithm authors:

Loos, M.; Gerber, C.; Corona, F.; Hollender, J.; Singer, H. (2015).
Accelerated isotope fine structure calculation using pruned transition trees,
*Analytical Chemistry*, 87(11), 5738-5744.
<https://doi.org/10.1021/acs.analchem.5b00941>

## License

GPL-2.0-only, like the original enviPat package. See the `LICENSE`, `AUTHORS`
and `NOTICE` files in the [repository](https://gitlab.eawag.ch/uchem/msf-server/pyenvipat).
