# twixreader: python reader for Siemens MRI raw data

Loosely based on mapVBVD (twix reader for matlab).

Warning: There are still a long list of missing features (e.g. protocol header information is currently not parsed), and the existing features are not tested very well.


## Installation
Enter directory and type:
```bash
python setup.py install
```
You may require super-user privileges. If this is a problem, I suggest you to install anaconda in user-space.

## Supported Files
Twix raw data files ('.dat') from Siemens' VB and VD/VE software lines are supported (including multi-raid files).

## How to use?
### Step 1: Parse file
Examples:
* use filename (relative or full path to file)
```python
from twixreader import Twix
twix = Twix(filename)
```


* use measurement ID (if unique and file in directory)
```python
twix = Twix(10)
```

Twix() returns a list with one entry per measurement of a multi-raid file (a single entry if not multi-raid). 
Each entry in the 'twix' list consists of a python dictionary with entries for each raw data type that was found in the measurement (e.g. 'ima', 'noise', 'phasecor', ...):


```python
print(twix)
Out[]: [{'noise': <twixreader.twixreader.generic_twix_object at 0x7f56def3b518>},
 {'noise': <twixreader.twixreader.generic_twix_object at 0x7f56def3b630>,
  'refscan': <twixreader.twixreader.generic_twix_object at 0x7f56def3b358>,
  'ima': <twixreader.twixreader.generic_twix_object at 0x7f56def3b3c8>}]
```


### Step 2: Read data

Examples:

* All noise data from first scan 
```python
noise = twix[0]['noise'][:]
```


* All image data from second scan (multi-raid file):
```python
imdata = twix[1]['ima'][:]
```


* All data for second phase-encoding line:
```python
imdata = twix[1]['ima'][...,1,:,:]
```



### Additional Features

* Control oversampling removal
```python
twix[0]['ima'].args.removeOS = True
```

* Print dimension order of virtual array
```python
print(twix[0]['noise'].dataDims)
Out: ['Ide', 'Idd', 'Idc', 'Idb', 'Ida', 'Seg', 'Set', 'Rep', 'Eco', 'Phs', 'Ave', 'Sli', 'Par', 'Lin', 'Cha', 'Col']
```


* Change order of virtual array (colfirst or collast)
```python
print(twix[0]['noise'].shape)
Out: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 128, 52, 256]
twix[0]['noise'].args.dataOrder = 'colfirst'
Out: [256, 52, 128, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
```


* Obtain all raw data in acquisition order:
```python
twix[0]['ima'].raw()
```


* Obtain slice of raw data in acquisition order:
```python
twix[0]['ima'].raw(slice(10, 50))
```


* Obtain information about line/partition acquisition order:
```python
lin_order = twix[0]['ima'].lin
par_order = twix[0]['ima'].par
```
The rest of the measurement data header (mdh) information is stored as well in other member variables, such as 'ave', 'rep', 'sliceData', 'icePara', 'freePara',...


## Authors
*Philipp Ehses* (philipp.ehses@dzne.de)
