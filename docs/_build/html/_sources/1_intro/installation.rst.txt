Installation
============

The simplest way to download the most recent stable version of GalPRIME is to use pip. 
You can install it by running the following command:

.. code-block:: bash

    pip install galprime

You can also install the current development tab through github:

.. code-block:: bash

    git clone https://github.com/HSouch/GalPRIME.git /path/to/dir
    cd /path/to/dir
    pip install .

.. note::

    You can include the `-e` tag to make the distribution editable.

To test that the code is working, you can run the following command in python:

.. code-block:: python

   import galprime as gp
   print(gp.__version__)

