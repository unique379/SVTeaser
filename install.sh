rm -rf dist/
python3 setup.py sdist bdist_wheel
python3 -m pip install dist/SVTeaser*.tar.gz

#git clone https://github.com/fritzsedlazeck/SURVIVOR.git
#cd SURVIVOR/Debug
#make
#export PATH=$PATH:`pwd`/
#
#echo "Place $(pwd) in your PATH"


