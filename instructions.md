## Access the CRCM5-CMIP6 data
The CRCM5-CMIP6 data set can be easily downloaded using the following python script:
```python
import threddsclient
from pathlib import Path
outfolder = Path('/desired/output/path/') # save to current directory
url = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/catalog/birdhouse/disk2/ouranos/CORDEX/CMIP6/DD/NAM-12/OURANOS/ERA5/evaluation/r1i1p1/CRCM5/v1-r1/day/tasmax/catalog.html"
with open(r"/path/to/bash/script/test.sh", 'w') as f:
    for ds in threddsclient.crawl(url, depth=10):
        outfile = outfolder.joinpath('CORDEX', ds.download_url().split('/CORDEX/')[-1])
        cmd = f"wget --limit-rate=20m {ds.download_url()} -P {outfile.parent.as_posix()}"
        print(cmd)
        f.write(f"{cmd}\n")
```
The python script generates a `test.sh` bash script containing the wget commands needed to download the variable listed by the url. The structure of the url is the following:

https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/catalog/birdhouse/disk2/ouranos/CORDEX/CMIP6/DD/NAM-12/OURANOS/{PILOT}/{SCENARIO}/r1i1p1/CRCM5/v1-r1/{FREQUENCY}/{VARIABLE}/catalog.html

To naviguate the different combinations of pilot, scenarios, frequency and variables, visit this url https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/catalog/birdhouse/disk2/ouranos/CORDEX/CMIP6/DD/NAM-12/OURANOS/catalog.html

Once the desired url is selected, simply run the python script above and then run the generated bash script. Please leave the `--limit-rate=20m` option as is to ensure that you are not overwelming our servers.