NOT PART OF CURRENT CHANNELFLOW

 What I did: 
  The header flowfield.h now contains several references to CUDA as well.
  Every flowfield knows also arrays on the GPU, if cmake defined 
  HAVE_CUDA, i.e. the CUDA libs have been found. 
  All in all things are rather similar to the use of cfmpi(), there is 
  a call usecuda() which is defined in any case, but should only be true
  if CUDA is working. Otherwise there should be an error. 

  For the integration of CUDA I touched the following files:
  - flowfield.h
  - flowfield.cpp
  - benchmark.cpp
  - dns.cpp
  - diffops.cpp
  - chebychev.cpp, maybe
  - and cuda_flowfield.cu, obviously

  First, I just did the ffts via cufft, but as every makeSpectral and the 
  like contained two memcopy executions, this was not faster than a 
  single thread.
  Second, I took navierstokesNL and ported it entirely to CUDA, including
  the routines curl and cross. This resulted in a speedup of about 2fold
  compared to a single thread. Not embarassing, but a bit.

  ToDo list
  - have usecuda controlled from the calling function and not the source 
  code of the flowfield constructor. Maybe adding a constructor parameter.
  - get the results right. There are still huge bugs in the calculation.
  - port the linear part of the time step to CUDA
  - extend to other algorithms than the defaults of benchmark
