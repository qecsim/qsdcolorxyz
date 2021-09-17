qsdcolorxyz
===========

.. image:: https://github.com/qecsim/qsdcolorxyz/workflows/CI/badge.svg?branch=main
    :target: https://github.com/qecsim/qsdcolorxyz/actions?workflow=CI
    :alt: CI Status

.. image:: https://codecov.io/gh/qecsim/qsdcolorxyz/branch/main/graph/badge.svg?token=ZF3QNFIN9J
    :target: https://codecov.io/gh/qecsim/qsdcolorxyz
    :alt: Coverage

**qsdcolorxyz** is a Python 3 package that extends `qecsim`_ with additional
components for a variant of the color code.

.. _qecsim: https://github.com/qecsim/qecsim

The qecsim package is a quantum error correction simulator, which can be
extended with new codes, error models and decoders that integrate into its
command-line interface. This qsdcolorxyz package includes implementations of a
variant of the color code and a corresponding tensor-network decoder.


Installation
------------

TODO


Usage
-----

TODO


Developer Notes
_______________

Tasks for running tests with coverage, checking style, generating documentation
and building source and binary distributables can be executed using tox_. See
`<./tox.ini>`__ for more details.

.. _tox: https://tox.readthedocs.io/

For example, distributables can be built as follows:

.. code-block:: text

    $ source venv/bin/activate              # activate venv (Windows: venv\Scripts\activate)
    (venv) $ tox -ebuild                    # build qsdcolorxyz distributables
    ...
    (venv) $ ls ./dist/                     # list qsdcolorxyz distributables
    qsdcolorxyz-0.1b1-py3-none-any.whl    qsdcolorxyz-0.1b1.tar.gz


License / Citing
----------------

qsdcolorxyz is released under the BSD 3-Clause license; see `<LICENSE>`__.

If you use qecsim in your research, please see the `qecsim documentation`_ for
citing details.

.. _qecsim documentation: https://qecsim.github.io/


Links
-----

* Source code: https://github.com/qecsim/qsdcolorxyz
* qecsim source code: https://github.com/qecsim/qecsim
* qecsim documentation: https://qecsim.github.io/
* Contact: qecsim@gmail.com

----

Copyright 2021, David K. Tuckett.
