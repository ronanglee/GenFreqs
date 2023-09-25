# GenFreqs
Calculate vibrational frequencies from a user-supplied Hessian and geometry files

Note  
- Only works for non-linear molecules
- Should not be used in production calculations
- No format checking is done on input files

  ## Installation
  Firstly, clone into the repository  

  `# git clone https://github.com/ronanglee/GenFreqs.git`

  Make the build folder  

  `# mkdir build`

  Enter build folder  

  `# cd build`

  Generate make files

  `# cmake ..`

  Compile source files and link to the final executable

  `# make`

  ## Usage

  The `GenFreq` executable needs to be run with two other arguments, `example.dat` `example.hess` containing
  the input geometry and Hessian respectively  

  `# GenFreq example.dat example.hess`
  

  
