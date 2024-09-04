# obstools_ipgp

Tools for evalauating and manipulating obs data.

Data must be in SDS format and metadata in StationXML format for all of these tools


## Command-line programs

Type ``{command} -h`` to get a list of parameters and options

### Existing

| Program     | description                                             |
| ----------- | ------------------------------------------------------- |
| plotPSDs    | plot power spectral densities of all stations           |
| obstest     | plot different representations of tests (noise, taps... |

### Future

| Program      | description                                                   |
| ------------ | ------------------------------------------------------------- |
| data_extent  | plot/list the extent of data for all channels????             |
| drift_correl | Use inter-station cross-correlations to calculate clock drift |
| to_SEGY      | Transform data + shotfiles into SEGY                          |


For details, see the [documentation](obstools_ipgp.readthedocs.io)
