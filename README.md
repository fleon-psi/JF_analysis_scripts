# JF_analysis_scripts
Scripts to analyse data obtained with JUNGFRAU detector and processed with XDS

# System requirements

* Tested on Linux CentOS 6, 7 and MacOS X
* Python 2.7 with standard scientific packages (matplotlib, numpy, math, sys).
* If python is installed, should also operate on Windows machine.

# Installation guide

```
git clone https://github.com/fleon-psi/JF_analysis_scripts
```

# Running

## Rpixel.py

```
python Rpixel.py XDS_ASCII.HKL 50 2.8 20
```

```
python Rpixel.py XDS_ASCII.HKL 50 2.8 20
```


## pixelmap.py

```
python pixelmap.py XDS_ASCII.HKL 50 2.0
```

## CorrDataset.py

```
python CorrDataset.py dataset1/INTEGRATE.HKL dataset2/INTEGRATE.HKL
```


## Rmeas.py

```
python Rmeas.py XDS_ASCII.HKL 50 2.8 20
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

# Authors
Filip Leonarski (Paul Scherrer Institut)
Meitian Wang (Paul Scherrer Institut)

## Acknowledgments

* Kay Diederichs (Uni Konstanz)
