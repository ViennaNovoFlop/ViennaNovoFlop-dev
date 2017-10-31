Two-dimensional periodic boundary conditions (2D PBC) code for OOMMF.

Version: April-2011-a, dependencies: OOMMF release 1.2.0.4, must older than snapshot 2009.12.16. 

Version: April-2011-b, dependencies: OOMMF release 1.2.0.4, snapshot 2009.12.16. 


The basic assumption of 2D PBC is that the "primary" sample simulation window 
was repeated infnitely many times in one plane (xy plane was assumed here), each 
repeated "image" have the same magnetization with the primary window.

The mainly changed parts for 2D PBC are demagnetization interaction and exchange 
interaction, which according to PBC_Demag_2D and PBC_Exchange_2D in this code, respectively.


For any suggestion please contact me, email:

	wangww2011@gmail.com, or wangww08@lzu.edu.cn

thank you!

Install:

  1) Please remove the previous 2D PBC codes in oommf/app/oxs/local (if exists).

  2) Copy/Move 2D PBC code into oommf/app/oxs/local.

  3) Run commands

         tclsh oommf.tcl pimake distclean

     and 

         tclsh oommf.tcl pimake 

     from the oommf root directory oommf/ .
 

Usages:

# Demag
Specify PBC_Demag_2D {
  tensor_file_name ..
  tensor_error ..
  sample_repeat_number ..
  asymptotic_radius ..
  dipolar_radius ..
}

where

   tensor_file_name gives the name of expected file to be read/written the demagnetization tensors.

   tensor_error gives the maximum error allowed at calculating the demagnetization tensor, default value is 1e-10. The code will estimate the sample_repeat_number automatically according to this option.

   sample_repeat_number gives the number of the sample expected to be repeated while calculating the demagnetization tensors. If this option is given then the option tensor_error is invalid.

	 asymptotic_radius gives the critical radius within which the demagnetization tensors are calculated using Newell's formula. Default value of asymptotic_radius is 32.0 .

   dipolar_radius gives the critical radius within which the demagnetization tensors are calculated using the asymptotic methods. Default value is very large, 10000.0 current now, just to disable the usage of direct dipolar method.


For instance,

Specify PBC_Demag_2D {
  tensor_file_name just_test
  tensor_error 1e-11
}

# Exchange
Specify PBC_Exchange_2D {
  A   13e-12
}

Jsut uniform six neighbor exchange energy on rectangular mesh for 2D PBC is implemented, now.
