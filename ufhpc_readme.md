**Setup instructions on UF HPG2**

* Check that you have python3/3.5.0 loaded in your modules


* In your directory of choice under your home, install the virtualenv with
  virtualenv -p python3 <venv_name>

  Lets call it mpint_venv

* Activate it with
  source mpint_venv/bin/activate

* pip install --ignore pymatgen==3.7.1

Finally,

* python hpg2_setup.py develop


If this does not work please contact me at joshgabriel92@gmail.com