**Setup instructions on UF HPG2**

* Check that you have python2.7.8 modules loaded on ufhpc

* In your directory of choice under your home, install the virtualenv with
  python virtualenv.xxx/virtualenv.py mpint_venv as mentioned in the setup section 
  of README.md

  Lets call it mpint_venv

* Activate it with
  source mpint_venv/bin/activate

* pip install --ignore pymatgen==3.7.1

* git clone https://github.com/joshgabriel/MPInterfaces.git

Finally,

* cd MPInterfaces

* git checkout ufhpc_py27_compat

* python hpg2_setup.py develop

* cp config.yaml mpinterfaces/config_mine.yaml

* update the config_mine.yaml 


-----


If you want to use the latest features (mostly post-processing) 
of pymatgen (>v4.0) like pymatgen4.6, you will need to create
a separate virtual environment for python3 on ufhpc by doing 

* module load python3/3.5.0 

* virtualenv -p python3 <venv_name>

* source mpint_venv/bin/activate

NOTE: spglib still gives library conflicts, so any post-processing
on ufhpc with symmetry operations may still give errors. 


-----

pymatgen now requires python3.5 and above and so development on 
MPInterfaces proceeds also with this requirement with limited 
backward compatibility support confined only to this branch. 

On ufhpc the temporary solution is: 
* submit jobs using this branch 
* do higher order post-processing jobs in your own python3.5 environment
  preferrably on your own computer 

A BYOE permanent solution is in the works. 

If this does not work please contact me at joshgabriel92@gmail.com
