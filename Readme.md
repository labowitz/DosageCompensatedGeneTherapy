# Synthetic dosage-compensating miRNA circuits for quantitative gene therapy

This git repository contains the code related to our paper "Synthetic
dosage-compensating miRNA circuits for quantitative gene therapy".

This code package contains C++ code that is compiled together using
`cmake`. The code is dependent on Windows-specific concurrency
libraries, and so must be compiled on a Windows machine. There are 3
executables produced: `analysisU2OSHCR2023-07-11_2.exe`, the code for
analyzing the in-vitro smFISH data contained in the
`U2OSFishvsHCR2023-06-29` folder used for constructing figure 2;
`analysisBrain2023-10-30.exe`, the code for analyzing ectopic vs
endogenous MeCP2 expression in the folder `BrainHCR2023-10-26` which
is featured in figure 3; and `analysisBrain2023-11-04.exe` which
analyzes the data in `BrainHCR2023-11-04` which is used in Figure 4.

The code depends on `libTIFF` and `openCV4` that must be installed using
`vcpckg`. To install vcpckg:

```
git clone https://github.com/microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.bat
```

Then to install the dependencies:

```
.\vcpkg\vcpkg install --triplet=x64-windows opencv4 
```

The `opencv4` library requires `tiff` which will be installed alongside
it. 

To check if the installation has succeeded, the directory

`./vcpkg/packages/opencv4_x64-windows/share/opencv/`

should be created and populated with `.cmake` files. Additionally, the
directories

```
./vcpkg/packages/tiff_x64-windows/include
./vcpkg/packages/tiff_x64-windows/lib
```

should be populated with header and lib files, respectively.

To build the system, run:

`cmake -DCMAKE_TOOLCHAIN_FILE=.\vcpkg\scripts\buildsystems\vcpkg.cmake .`

This should produce a file called `Analysis.sln` which can be opened
with Visual Studio. In Visual studio, change the build type from
'Debug' to 'Release' and then use `Build > Build Solution`. 

This should then create the files 

```
./Release/analysisBrain2023-10-30.exe
./Release/analysisBrain2023-11-04.exe
./Release/analysisU2OSHCR2023-07-11_2.exe 
```

which can be called directly to run the analysis code on the images,
provided that this repository is a subdirectory of the Dosage
Compensated Gene Therapy magic folder. 
