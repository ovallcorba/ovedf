# ovedf: Convert EDF image to 1D-XRD data pattern

`ovedf.py` is a small python script to convert EDF images to 1D XRD patterns (i.e. powder pattern). It allows several options such as external calibration files, consider sever azimuthal bins, use of masks, etc...

The script is used in the ALBA-MSPD beamline with the Rayonix SX165 detector that uses Lima to generates the EDF data images. Although the defaults are set for this beamline conditions, it can be adapted elsewhere. It can use instrumental calibration from D2Dplot or Fit2D.

*Note:* Matplotlib is needed if you want the script to also plot the data)

## Author(s)

  - **Oriol Vallcorba** (ovallcorba@cells.es)

## Disclaimer

This software is distributed WITHOUT ANY WARRANTY. The authors (or their institutions) have no liabilities in respect of errors in the software, in the documentation and in any consequence of erroneous results or damages arising out of the use or inability to use this software. Use it at your own risk.

## License

This project is licensed under the [GPL-3.0 license](LICENSE.txt)
