qsdcolorxyz
===========

.. image:: https://img.shields.io/badge/docs-dev-blue.svg
    :target: https://qecsim.github.io/qsdcolorxyz
    :alt: Documentation

.. image:: https://github.com/qecsim/qsdcolorxyz/workflows/CI/badge.svg?branch=main
    :target: https://github.com/qecsim/qsdcolorxyz/actions?workflow=CI
    :alt: CI Status

.. image:: https://codecov.io/gh/qecsim/qsdcolorxyz/branch/main/graph/badge.svg?token=ZF3QNFIN9J
    :target: https://codecov.io/gh/qecsim/qsdcolorxyz
    :alt: Coverage

**qsdcolorxyz** is a Python 3 package that extends `qecsim`_ with additional
components for the XYZ variant of the color code.

.. _qecsim: https://github.com/qecsim/qecsim

The qecsim package is a quantum error correction simulator, which can be
extended with new codes, error models and decoders that integrate into its
command-line interface. This qsdcolorxyz package includes implementations
relevant to the XYZ variant of the color code.


Installation
------------

* Single-line install and upgrade using `pip`_:

.. code-block:: text

    $ pip install -U git+https://github.com/qecsim/qsdcolorxyz.git#egg=qsdcolorxyz

.. _pip: https://pip.pypa.io/en/stable/quickstart/

* Recommended install into a virtual environment:

.. code-block:: text

    $ python3 --version                     # qecsimext requires Python 3.5+
    Python 3.9.5
    $ python3 -m venv venv                  # create virtual environment
    $ source venv/bin/activate              # activate venv (Windows: venv\Scripts\activate)
    (venv) $ pip install -U setuptools pip  # install / upgrade setuptools and pip
    ...
    Successfully installed pip-21.3 setuptools-58.2.0
    (venv) $ # install qsdcolorxyz from github repository
    (venv) $ pip install git+https://github.com/qecsim/qsdcolorxyz.git#egg=qsdcolorxyz
    ...
    Successfully installed ... qsdcolorxyz-0.1b1 ...
    (venv) $ deactivate                     # deactivate venv
    $

* Logging configuration (optional): See `qecsim logging configuration`_ for details.

.. _qecsim logging configuration: https://qecsim.github.io/installation.html#logging-configuration-optional


Usage
-----

* Display run help via console script:

.. code-block:: text

    $ source venv/bin/activate              # activate venv (Windows: venv\Scripts\activate)
    (venv) $ qecsim run --help
    Usage: qecsim run [OPTIONS] CODE ERROR_MODEL DECODER ERROR_PROBABILITY...

      Simulate quantum error correction.

      Arguments:

       CODE                  Stabilizer code in format name(<args>)
        ...
        color666xyz           Color 6.6.6 XYZ (size INT odd >=3)
        ...

       ERROR_MODEL           Error model in format name(<args>)
        ...
        generic.depolarizing  Pr I,X,Y,Z is 1-p,p/3,p/3,p/3
        ...

       DECODER               Decoder in format name(<args>)
        ...
        color666xyz.mps       MPS ([chi] INT, ...)
        ...

       ERROR_PROBABILITY...  One or more probabilities as FLOAT in [0.0, 1.0]
    ...

* Run simulation example via console script:

.. code-block:: text

    $ source venv/bin/activate              # activate venv (Windows: venv\Scripts\activate)
    (venv) $ qecsim run -r100 "color666xyz(3)" "generic.depolarizing" "color666xyz.mps(36)" 0.1
    ...
    [{"code": "Color 6.6.6 XYZ 3", "custom_totals": null, "decoder": "Color 6.6.6 XYZ MPS (chi=36)", "error_model": "Depolarizing", "error_probability": 0.1, "error_weight_pvar": 0.4356, "error_weight_total": 62, "logical_failure_rate": 0.08, "measurement_error_probability": 0.0, "n_fail": 8, "n_k_d": [7, 1, 3], "n_logical_commutations": [7, 4], "n_run": 100, "n_success": 92, "physical_error_rate": 0.08857142857142858, "time_steps": 1, "wall_time": 0.38195787700000006}]

* Run simulation example via module script with Python optimize flag:

.. code-block:: text

    $ source venv/bin/activate              # activate venv (Windows: venv\Scripts\activate)
    (venv) $ python3 -O -m qecsim run -r100 "color666xyz(3)" "generic.depolarizing" "color666xyz.mps(36)" 0.1
    ...
    [{"code": "Color 6.6.6 XYZ 3", "custom_totals": null, "decoder": "Color 6.6.6 XYZ MPS (chi=36)", "error_model": "Depolarizing", "error_probability": 0.1, "error_weight_pvar": 0.5416, "error_weight_total": 72, "logical_failure_rate": 0.09, "measurement_error_probability": 0.0, "n_fail": 9, "n_k_d": [7, 1, 3], "n_logical_commutations": [9, 4], "n_run": 100, "n_success": 91, "physical_error_rate": 0.10285714285714287, "time_steps": 1, "wall_time": 0.345211924}]


Developer Notes
---------------

If you want to make code changes to qsdcolorxyz, you can `clone`_ or `fork`_
the repository.

.. _clone: https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository
.. _fork: https://docs.github.com/en/github/getting-started-with-github/fork-a-repo

Developer tasks for running tests with coverage, checking style, generating
documentation and building source and binary distributables can be executed
using tox_. See `<./tox.ini>`__ for more details.

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
* Documentation: https://qecsim.github.io/qsdcolorxyz
* qecsim source code: https://github.com/qecsim/qecsim
* qecsim documentation: https://qecsim.github.io/
* Contact: qecsim@gmail.com

----

Copyright 2021, David K. Tuckett.
