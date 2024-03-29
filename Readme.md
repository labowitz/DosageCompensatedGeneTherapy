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
.\vcpkg\vcpkg install --triplet=x64-windows-static opencv4 
```

The `opencv4` library requires `tiff` which will be installed alongside
it. 

To check if the installation has succeeded, the directory

`./vcpkg/packages/opencv4_x64-windows-static/share/opencv/`

should be created and populated with `.cmake` files. Additionally, the
directories

```
./vcpkg/packages/tiff_x64-windows-static/include
./vcpkg/packages/tiff_x64-windows-static/lib
```

should be populated with header and lib files, respectively.

To build the system, run:

`cmake -DCMAKE_TOOLCHAIN_FILE=.\vcpkg\scripts\buildsystems\vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows-static -S . -B ./build -G "Visual Studio 17 2022"`

`cmake --build build --config Release`

This will then create the files 

```
./build/Release/analysisBrain2023-10-30.exe
./build/Release/analysisBrain2023-11-04.exe
./build/Release/analysisU2OSHCR2023-07-11_2.exe 
```

which can be called directly to run the analysis code on the
images. The files are dependent on the Nd2SDK `.dll` files in the root
directory, so these programs should be called from the root.


## File descriptions

`ImageAnalysis.hpp`: an image analysis library 

`analysisU2OSHCR2023-07-11_2.cpp`: analysis to generate data from
images for figure 2 

`analysisBrain2023-10-30.cpp`: analysis to generate data from images
for Figure 3

`analysisBrain2023-11-04.cpp`: analysis to generate data from images for
Figure 4

All other files are support files.
