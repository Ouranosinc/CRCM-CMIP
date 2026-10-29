## Access the CRCM5-CMIP6 data

### `intake-esgf`
The [intake-esgf](https://intake-esgf.readthedocs.io/en/latest/) python package allows an easy search to the ESGF indexes and will download the files in the background upon request. For example, to search, download and open as an xarray dataset the `tasmax` dataset of the evaluation run, one can do:
```python
from intake_esgf import ESGFCatalog

cat = ESGFCatalog()
res = cat.search(project='CORDEX-CMIP6', variable_id='tasmax', driving_experiment_id='evaluation', source_id='CRCM5-SN', version_realization='v1-r1', frequency='day')
dsd = cat.to_dataset_dict()
```
Please see the package's documentation for more.

### Generate a custom wget script
The CRCM5-CMIP6 data set can be easily downloaded using the following python script:
```python
import threddsclient
from pathlib import Path
outfolder = Path('/desired/output/path/') # save to current directory
url = "https://esgf.ouranos.ca/thredds/catalog/mrcc5-cmip6/CORDEX-CMIP6/DD/NAM-12/OURANOS/ERA5/evaluation/r1i1p1f1/CRCM5-SN/v1-r1/day/tasmax/v20260318/catalog.html"
with open(r"/path/to/bash/script/test.sh", 'w') as f:
    for ds in threddsclient.crawl(url, depth=10):
        outfile = outfolder.joinpath('CORDEX-CMIP6', ds.download_url().split('/CORDEX-CMIP6/')[-1])
        cmd = f"wget --limit-rate=20m {ds.download_url()} -P {outfile.parent.as_posix()}"
        print(cmd)
        f.write(f"{cmd}\n")
```
The python script generates a `test.sh` bash script containing the wget commands needed to download the variable listed by the url. The structure of the url is the following:

https://esgf.ouranos.ca/thredds/catalog/mrcc5-cmip6/CORDEX-CMIP6/DD/NAM-12/OURANOS/{PILOT}/{SCENARIO}/{DRIVING_MEMBER}/CRCM5-SN/{VERSION_RELIZATION}/{FREQUENCY}/{VARIABLE}/{FILE_VERSION}/catalog.html"

To naviguate the different combinations of pilot, scenarios, frequency and variables, visit the root catalog : https://esgf.ouranos.ca/thredds/catalog/mrcc5-cmip6/CORDEX-CMIP6/catalog.html

Once the desired url is selected, simply run the python script above and then run the generated bash script. Please leave the `--limit-rate=20m` option as is to ensure that you are not overwelming our servers.
