.. highlight:: shell

.. role:: bash(code)
   :language: bash

Installation
------------




System requirements
>>>>>>>>>>>>>>>>>>>

* Linux/Unix
* Python >= 3.8


We recommend to create an independent conda environment for STRIDE. If users do not have conda, please install Miniconda first:
::
   
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh


Install the stable version
>>>>>>>>>>>>>>>>>>>>>>>>>>

Step 1 Prepare conda environment for STRIDE.
::::::::::::::::::::::::::::::::::::::::::::
:: 

   conda create -n stride python=3.8
   conda activate stride

Step 2 Install STRIDE package from :bash:`pypi`.
::::::::::::::::::::::::::::::::::::::::::::::::
::

   pip install -U stride


Install the developing version
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Step 1 Prepare conda environment for STRIDE.
::::::::::::::::::::::::::::::::::::::::::::
:: 

   conda create -n stride python=3.8
   conda activate stride

Step 2 Download STRIDE package from github.
:::::::::::::::::::::::::::::::::::::::::::
::

   git clone https://github.com/DongqingSun96/STRIDE.git

Step 3 Install dependencies of STRIDE.
::::::::::::::::::::::::::::::::::::::
::

   cd STRIDE
   pip install -r requirements.txt

Step 4 Install STRIDE.
::::::::::::::::::::::
::
  
   python setup.py install