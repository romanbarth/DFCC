 This section includes third-party license information for certain 
 third-party functions included with the 

 Dense Flow reConstruction and Correlation (DFCC)
 ----------------------------------------------------------------------- 

 Reference to the publication:
   Haitham A Shaban, Roman Barth, Kerstin Bystricky; Formation of correlated chromatin domains at nanoscale dynamic resolution during transcription, Nucleic Acids Research, , gky269, https://doi.org/10.1093/nar/gky269

 developed at:  
       Laboratoire de Biologie Moléculaire Eucaryote (LBME), 
       Centre de Biologie Intégrative (CBI), CNRS; 
       University of Toulouse, UPS; 31062 
       Toulouse; France


for the following files:

Optical Flow
- estimate_flow_interface.m
- load_of_method.m
- scale_image.m
- compute_flow.m
- structure_texture_decomposition_rof.m
- compute_image_pyramid.m
- resample_flow.m
- compute_flow_base.m
- robust_function.m
- generalized_charbonnier.m
- quadratic.m
- partial_deriv.m
- interp2_bicubic.m
- flow_operator.m
- make_convn_mat.m
- convmtxn.m
- deriv_over_x.m
- detect_occlusion.m
- denoise_color_weighted_medfilt2.m
- weighted_median.m

Others
- bfilter2.m
- inpaint_nans.m

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
- bfilter2.m
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Copyright (c) 2006, Douglas Lanman 
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

* Redistributions of source code must retain the above copyright 
notice, this list of conditions and the following disclaimer. 
* Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in 
the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
- inpaint_nans.m
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Copyright (c) 2009, John D'Errico 
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

* Redistributions of source code must retain the above copyright 
notice, this list of conditions and the following disclaimer. 
* Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in 
the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Modified or unmodified files within the Optical Flow determination are the following:
- estimate_flow_interface.m
- load_of_method.m
- scale_image.m
- compute_flow.m
- structure_texture_decomposition_rof.m
- compute_image_pyramid.m
- resample_flow.m
- compute_flow_base.m
- robust_function.m
- generalized_charbonnier.m
- quadratic.m
- partial_deriv.m
- interp2_bicubic.m
- flow_operator.m
- make_convn_mat.m
- convmtxn.m
- deriv_over_x.m
- detect_occlusion.m
- denoise_color_weighted_medfilt2.m
- weighted_median.m
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Authors: Deqing Sun, Department of Computer Science, Brown University
Contact: dqsun@cs.brown.edu
$Date: $
$Revision: $

Copyright 2007-2010, Brown University, Providence, RI. USA

                         All Rights Reserved

All commercial use of this software, whether direct or indirect, is
strictly prohibited including, without limitation, incorporation into in
a commercial product, use in a commercial service, or production of other
artifacts for commercial purposes.     

Permission to use, copy, modify, and distribute this software and its
documentation for research purposes is hereby granted without fee,
provided that the above copyright notice appears in all copies and that
both that copyright notice and this permission notice appear in
supporting documentation, and that the name of the author and Brown
University not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior permission.        

For commercial uses contact the Technology Venture Office of Brown University

THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
