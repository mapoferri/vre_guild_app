# VRE Tool Executor for GUILD

[![Documentation Status](https://readthedocs.org/projects/vre-template-tool/badge/?version=latest)](https://vre-template-tool.readthedocs.io/en/latest/?badge=latest)

Tool ready to run in the VRE matching the code in the [documentation](https://vre-template-tool.readthedocs.io/en/latest/?badge=latest).

Starting from the template for a VRE Tool, this is the adaptation for the GUILD - DisGenet draft workflow - retrieval for Gene Disease associations and generation of a PPI prioritazed network based on the results.

## Requirements

- Python 3.6 or later
- Conda env
- [Git](https://git-scm.com/downloads)

## Installation

Directly from GitHub:

```bash
cd $HOME
git clone https://github.com/inab/mapoferri/guild_app_vre.git
cd vre_template_tool
```

Create the Conda environment (attention to changing the paths accordingly):

```bash
source /anaconda3/etc/profile.d/conda.sh
conda env create -f environment.yml 
```

It may be possible to download/update the OpenVRE API directly in the Conda ENV:

```bash
pip install git+https://github.com/inab/openvre-tool-api.git
```

Or: 

```bash
pip install -r requirments.txt
```



## Run the Wrapper

```bash
./VRE_RUNNER --config tests/basic/config.json --in_metadata tests/basic/in_metadata.json --out_metadata out_metadata.json --log_file VRE_RUNNER.log
```

Look for the results in `tests/basic/run000/`.

To test the applications, all the attemps were made from the test/basics dir. To replicate them, access directly the directory. 

```bash
./test_VRE_RUNNER.sh 
```


## License
* © 2020-2022 Barcelona Supercomputing Center (BSC), ES

Licensed under the Apache License [Version 2.0](https://www.apache.org/licenses/LICENSE-2.0), see the file `LICENSE` for details.
