Reference datasets (``envipat.data``)
=====================================

.. py:module:: envipat.data

The package bundles the R package's reference datasets as plain CSV/TXT files.
Each dataset is a module attribute of :mod:`envipat.data` that is read from disk
on first access and cached afterwards, so importing the module is cheap:

.. code-block:: python

   import envipat.data as data

   data.isotopes.shape                      # (308, 5)
   data.resolution_list["Elite_R240000@400"]

The files were exported once from enviPat's R datasets, with floating point
columns written at full precision, so every value round-trips bit-exactly.
The datasets are not re-exported at the top level (``envipat.isotopes`` does not
exist), because that would load them at import time.

.. py:data:: isotopes
   :type: pandas.DataFrame

   Isotope masses and natural abundances, 308 rows, one per isotope; R's
   ``data(isotopes)``. Columns:

   ``element``
      Element symbol, e.g. ``"C"``. Labelled isotopes have their own element
      entries in brackets, e.g. ``"[13]C"``.
   ``isotope``
      Isotope name, e.g. ``"13C"``.
   ``mass``
      Exact isotope mass (Da).
   ``abundance``
      Natural abundance as a fraction (``1.0`` for the labelled entries).
   ``ratioC``
      Element-to-carbon ratio column of the R dataset; carried over unchanged
      and not used by the pipeline functions.

   This table is the ``isotopes`` argument of
   :func:`~envipat.formula.check_chemform`, :func:`~envipat.isopattern.isopattern`
   and :func:`~envipat.isowrap.isowrap`. Pass a modified copy to change
   abundances, e.g. for labelling experiments.

.. py:data:: chemforms
   :type: pandas.Series

   107 example chemical formulas (strings, Series name ``"chemform"``); R's
   ``data(chemforms)``. Useful as test input for the pipeline functions.

.. py:data:: adducts
   :type: pandas.DataFrame

   51 common ESI adducts, positive and negative mode; R's ``data(adducts)``.
   Columns:

   ``Name``
      Adduct name, e.g. ``"M+H"``, ``"2M-H"``.
   ``calc``
      The mass calculation as text, e.g. ``"M+1.007276"``.
   ``Charge``
      Ion charge (signed integer).
   ``Mult``
      Molecule multiplier (``2`` for dimers).
   ``Mass``
      Mass added to (or, if negative, removed from) ``Mult`` times the
      neutral mass (Da).
   ``Ion_mode``
      ``"positive"`` or ``"negative"``.
   ``Formula_add``, ``Formula_ded``
      Formula added to / deducted from the molecule, or the string
      ``"FALSE"`` if nothing is added / deducted. Use them with
      :func:`~envipat.formula.mergeform`, :func:`~envipat.formula.check_ded` and
      :func:`~envipat.formula.subform`.
   ``Multi``
      Multiplier applied to the molecular formula, for use with
      :func:`~envipat.formula.multiform`.

.. py:data:: resolution_list
   :type: dict[str, pandas.DataFrame]

   31 instrument resolution calibration curves (Orbitrap, Q Exactive,
   TOF, ...); R's ``data(resolution_list)``. Keyed by curve name, e.g.
   ``"Elite_R240000@400"``, in the R order. Each value is a two-column
   DataFrame: ``mz`` (the R column ``m/z``, renamed) and ``R``, the resolving
   power at that m/z. Pass one curve as ``resmass`` to
   :func:`~envipat.resolution.get_resolution`.

   The available curves:

   .. code-block:: python

      list(data.resolution_list)
