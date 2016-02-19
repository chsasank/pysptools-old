@echo on

cls

rem win32
rem set path=c:\tools\Python27-32;%path%
rem C:\tools\Python27-32\python.exe setup_sam.py build_ext --inplace

rem amd64
python setup_nfindr.py build_ext --inplace
rem python setup_atgp.py build_ext --inplace
