*********************
pyHiM - Documentation
*********************

**Date**: |today| **Version**: |release|

`pyHiM` is a software package developed by the `Nollmann Lab <http://www.nollmannlab.org>`_ at the `Center of Structural Biology <http://www.cbs.cnrs.fr>`_, a department of the `CNRS <http://www.cnrs.fr>`_ and the `INSERM <http://www.inserm.fr>`_. 
`pyHiM` implements the analysis of multiplexed DNA-FISH data, as described in our `protocols paper <https://github.com/NollmannLab/HiM_protocol>`_.

pyHiM was entirely written in python and makes extensive use of the following packages:

- `Astropy <https://www.astropy.org/>`_
- `scikit-image <https://scikit-image.org/>`_
- `starDist <https://github.com/stardist/stardist>`_
- `Big-FISH <https://github.com/fish-quant/big-fish>`_ to perform Gaussian 3D fits

TEMPO: `panda` template:

.. panels::
    :card: + intro-card text-center
    :column: col-lg-6 col-md-6 col-sm-6 col-xs-12 p-2

    ---
    :img-top: _static/index_getting_started.svg

    Getting started
    ^^^^^^^^^^^^^^^

    New to *pyHiM*? Check out the getting started guides. They contain an
    introduction to *pyHiM'* main concepts and links to additional tutorials.

    +++

    .. link-button:: getting_started/welcome
            :type: ref
            :text: To the getting started guides
            :classes: btn-block btn-secondary stretched-link

    ---
    :img-top: _static/index_user_guide.svg

    User guide
    ^^^^^^^^^^

    The user guide provides in-depth information on the
    key concepts of pyHiM with useful background information and explanation.

    +++

    .. link-button:: user_guide/User_Guide_pyHiM
            :type: ref
            :text: To the user guide
            :classes: btn-block btn-secondary stretched-link

    ---
    :img-top: _static/index_api.svg

    API reference
    ^^^^^^^^^^^^^

    The reference guide contains a detailed description of
    the pyHiM API. The reference describes how the methods work and which parameters can
    be used. It assumes that you have an understanding of the key concepts.

    +++

    .. link-button:: reference/pyhim_class_diagram
            :type: ref
            :text: To the reference guide
            :classes: btn-block btn-secondary stretched-link

    ---
    :img-top: _static/index_contribute.svg

    Developer guide
    ^^^^^^^^^^^^^^^

    Saw a typo in the documentation? Want to improve
    existing functionalities? The contributing guidelines will guide
    you through the process of improving pyHiM.

    +++

    .. link-button:: contributor/todolist_doc_pyhim
            :type: ref
            :text: To the development guide
            :classes: btn-block btn-secondary stretched-link



.. toctree::
   :caption: Getting Started
   :hidden:

   getting_started/welcome

.. toctree::
   :caption: User Guide
   :hidden:

   user_guide/pyhim_presentation
   user_guide/pyhim_installation
   user_guide/fundamental
   user_guide/script
   user_guide/tutorial
   user_guide/reporting_bug
   user_guide/license
   user_guide/User_Guide_pyHiM

.. toctree::
   :caption: Reference
   :hidden:

   reference/pyhim_class_diagram
   reference/infoList_comprehension

.. toctree::
   :caption: Contributor's Guide
   :hidden:

   contributor/dev_process
   contributor/todolist_doc_pyhim
   contributor/recovery_bigfish_update

