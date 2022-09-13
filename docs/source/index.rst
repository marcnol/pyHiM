*********************
pyHiM - Documentation
*********************

**Date**: |today| **Version**: |release|

`pyHiM` implements the analysis of multiplexed DNA-FISH data, as described in our `protocols paper <https://github.com/NollmannLab/HiM_protocol>`_.

=============
Main features
=============

* Register nuclei and barcode images to remove or minimize drift
* Segment nuclei
* Localize barcodes
* Construct chromatin traces
* Build single-trace and ensemble pair-wise distance maps
* Work with 2D or 3D images

==============
How to proceed
==============


.. panels::
    :card: + intro-card text-center
    :column: col-lg-6 col-md-6 col-sm-6 col-xs-12 p-2

    ---
    .. image:: _static/index_getting_started.svg
       :width: 50

    Getting started
    ^^^^^^^^^^^^^^^

    New to *pyHiM*? Check out the getting started guides. They contain an
    introduction to *pyHiM'* main concepts and links to additional tutorials.

    +++

    .. link-button:: getting_started/quick_install
            :type: ref
            :text: To the getting started guides
            :classes: btn-block btn-info stretched-link text-white

    ---
    .. image:: _static/index_user_guide.svg
       :width: 60

    User guide
    ^^^^^^^^^^

    The user guide provides in-depth information on the
    key concepts of *pyHiM* with background information and explanation.

    +++

    .. link-button:: user_guide/pyhim_presentation
            :type: ref
            :text: To the user guide
            :classes: btn-block btn-info stretched-link text-white

    ---
    .. image:: _static/index_api.svg
       :width: 56

    API reference
    ^^^^^^^^^^^^^

    The reference guide contains a detailed description of
    the *pyHiM* API. The reference describes how the methods work and which parameters can
    be used.

    +++

    .. link-button:: reference/pyhim_class_diagram
            :type: ref
            :text: To the reference guide
            :classes: btn-block btn-info stretched-link text-white

    ---
    .. image:: _static/index_contribute.svg
       :width: 50

    Developer guide
    ^^^^^^^^^^^^^^^

    Saw a typo in the documentation? Want to improve
    existing functionalities? The contributing guidelines will guide
    you through the process of improving *pyHiM*.

    +++

    .. link-button:: contributor/dev_env
            :type: ref
            :text: To the development guide
            :classes: btn-block btn-info stretched-link text-white


`pyHiM` is a software package developed by the `Nollmann Lab <http://www.nollmannlab.org>`_ at the `Center of Structural Biology <http://www.cbs.cnrs.fr>`_, a department of the `CNRS <http://www.cnrs.fr>`_ and the `INSERM <http://www.inserm.fr>`_. 


.. toctree::
   :caption: Getting Started
   :hidden:

   getting_started/quick_install
   getting_started/typical_run
   getting_started/configuration_file
   getting_started/personalised_run
   getting_started/full_pyHiM_run_2D

.. toctree::
   :caption: User Guide
   :hidden:

   user_guide/pyhim_presentation
   user_guide/pipeline_overview
   user_guide/fundamental
   user_guide/script
   user_guide/license

.. toctree::
   :caption: Reference
   :hidden:

   reference/pyhim_class_diagram
   reference/infoList_comprehension
   reference/apidoc/modules

.. toctree::
   :caption: Contributor's Guide
   :hidden:

   contributor/dev_env
   contributor/dev_installation
   contributor/coding_style
   contributor/dev_process
   contributor/good_commit
   contributor/build_doc
   contributor/issues
   contributor/work_with_apifish
   contributor/ressources

