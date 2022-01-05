# Data Ingestion 

## Dependencies for download script
#### OSM extraction setup (requires homebrew at https://brew.sh/)
```
brew install wget
brew install osmium-tool
```
#### Python environment setup (requires miniconda at https://docs.conda.io/en/latest/miniconda.html)
```
conda create --name geo_download_env python=3.9.7 --yes
source activate geo_download_env
conda install -c conda-forge pygeos=0.12.0 --yes
conda install -c conda-forge geopandas=0.10.2 --yes
conda install -c conda-forge urlpath --yes
```

## Steps to download OSM linestrings and GADMs boundaries and write to local directory

#### Activate environment, clone repo, cd into repo, make output directory
```
source activate geo_download_env
cd /users/repos/
git clone git@github.com:mansueto-institute/kblock.git
cd /users/repos/kblock/ingestion
mkdir -p /users/downloads/folder
```
#### Deploy download job (pass in output directory with -o argument)
> **Caution:** Running this script will download LOTS of data for all countries in the `gadm_list=()` in the `geofabrik_list=()` in `download_job.sh`.
```
bash download_job.sh -o /users/downloads/folder
```

## Country array arguments in `download_job.sh`

##### Potential country codes for GADM array
```
gadm_list=(DZA AGO BEN BWA BFA BDI CPV CMR CAF TCD COM COG CIV COD DJI EGY GNQ ERI SWZ ETH GAB GMB GHA GIN GNB KEN LSO LBR LBY MDG MWI MLI MRT MUS MAR ESH MOZ NAM NER NGA RWA SHN STP SEN SYC SLE SOM ZAF SSD SDN TZA TGO TUN UGA ZMB ZWE AFG ARM AZE BHR BGD BTN BRN MYS SGP KHM CHN IND IDN IRQ IRN PSE ISR JPN JOR KAZ KGZ LAO LBN MDV MNG MMR NPL PRK PAK PHL KOR LKA SYR TWN TJK THA TKM UZB VNM YEM ALB AND AUT BLR BEL BIH BGR HRV CYP CZE DNK EST FRO FIN FRA GEO DEU GRC HUN ISL IRL IMN ITA KOS LVA LIE LTU LUX MLT MDA MCO MNE NLD MKD NOR POL PRT ROU RUS SRB SVK SVN ESP SWE CHE TUR UKR GBR CAN GRL MEX USA AUS COK FJI KIR MHL FSM NRU NCL NZL NIU PLW PNG WSM SLB TON TUV VUT ARG BOL BRA CHL COL ECU PRY PER SUR URY VEN)
```

##### Potential country names for Geofabrik array 
```
geofabrik_list=(algeria angola benin botswana burkina-faso burundi cape-verde cameroon central-african-republic chad comores congo-brazzaville ivory-coast congo-democratic-republic djibouti egypt equatorial-guinea eritrea swaziland ethiopia gabon senegal-and-gambia ghana guinea guinea-bissau kenya lesotho liberia libya madagascar malawi mali mauritania mauritius morocco morocco mozambique namibia niger nigeria rwanda saint-helena-ascension-and-tristan-da-cunha sao-tome-and-principe senegal-and-gambia seychelles sierra-leone somalia south-africa south-sudan sudan tanzania togo tunisia uganda zambia zimbabwe afghanistan armenia azerbaijan gcc-states bangladesh bhutan malaysia-singapore-brunei malaysia-singapore-brunei malaysia-singapore-brunei cambodia china india indonesia iraq iran israel-and-palestine israel-and-palestine japan jordan kazakhstan kyrgyzstan laos lebanon maldives mongolia myanmar nepal north-korea pakistan philippines south-korea sri-lanka syria taiwan tajikistan thailand turkmenistan uzbekistan vietnam yemen albania andorra austria belarus belgium bosnia-herzegovina bulgaria croatia cyprus czech-republic denmark estonia faroe-islands finland france georgia germany greece hungary iceland ireland-and-northern-ireland isle-of-man italy kosovo latvia liechtenstein lithuania luxembourg malta moldova monaco montenegro netherlands macedonia norway poland portugal romania russia serbia slovakia slovenia spain sweden switzerland turkey ukraine great-britain canada greenland mexico usa australia cook-islands fiji kiribati marshall-islands micronesia nauru new-caledonia new-zealand niue palau papua-new-guinea samoa solomon-islands tonga tuvalu vanuatu argentina bolivia brazil chile colombia ecuador paraguay peru suriname uruguay venezuela)
```
