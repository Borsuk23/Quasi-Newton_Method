========================================================================
    CONSOLE APPLICATION : CLAPACK-EXAMPLE Project Overview
========================================================================

CLAPACK-EXAMPLE.vcproj
    This is the main project file for VC projects

CLAPACK-EXAMPLE.c
    This is the main application source file.
    It coontains a very basic dgesv call using the CLAPACK library

========================================================================
Please copy the library into the lib directory 
or 
change the properties "Linker > General > Additional Library Directory" to where the libraries are

Libraries needed: a BLAS library, libf2c and lapack.lib (clapack) 

Do not forget to add any dll that is needed to run the program.
/////////////////////////////////////////////////////////////////////////////
Feedback and question welcome on the LAPACK forum
http://icl.cs.utk.edu/lapack-forum/index.php

/////////////////////////////////////////////////////////////////////////////
