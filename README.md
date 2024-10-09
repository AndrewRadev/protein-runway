## Setup

```.sh-session
# Set up and activate a virtualenv
python -mvenv .venv
source .venv/bin/activate

# Install requirements
pip install -r requirements.txt
```

## Run

For now, just a simple test script taken from <https://github.com/MDAnalysis/mdanalysis>

```.sh-session
# Move "crambin_md.pdb" and "crambin_md.xtc" into "data/"
python md_test.py
```
