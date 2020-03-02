# DimmiLitho

Include pixel-based mask synthesis, imaging model for optical lithography in Python. This is a simplified package, maybe OUT-OF-DATA.

# USAGE

1. These files base on PIL, Numpy, Sci, pyFFTW packages. GDSii package is already included
2. It can be run directly or included as a lib
3. Just run this test.py for inverse mask synthesis application
3. The NanGateLibGDS was accessed in (http://projects.si2.org/openeda.si2.org/projects/nangatelib)

# NOTE:

1. This package only contains files applied for scalar optics, thin mask model assumption
2. Image Equations refer to A.K.Wong's book
3. Process Equations refer to C.Mack's book
4. Opimization Equations can be found in the paper: 
   [1] Wen Lv, Shiyuan Liu, Xiaofei Wu, and Edmund Y. Lam, “Illumination source optimization in optical lithography via derivative-free optimization,” JOSAA 31(12): B19-B26 (2014). 
   [2] Wen Lv, Edmund Y. Lam, Haiqing Wei, and Shiyuan Liu, “Cascadic multigrid algorithm for robust inverse mask synthesis in optical lithography,” JM3 13(2): 023003 (2014).
   [3] Wen Lv, Qi Xia, and Shiyuan Liu, “Mask-filtering-based inverse lithography,” JM3 12(4): 043003 (2013).
   [4] Wen Lv, Shiyuan Liu, Qi Xia, Xiaofei Wu, Yijiang Shen, and Edmund Y. Lam, “Level-set-based inverse lithography for mask synthesis using the conjugate gradient and an optimal time step,” JVSTB 31(4): 041605(2013).


# Models:

## Mask

1. You can just run the mask.py file
2. Basd on the GDSii Python Packages (BY Eugeniy Meshcheryakov), slightly modified. It can convert to a pixel image from GDSII
3. The method openGDS file works well with the test GDS file in folder NanGateLibGDS, but it does not be verified for others

## Lens

1. You can just run the lens.py file


## Source


## Image 


## ILT 


# Notes

GDS needs to be flat

