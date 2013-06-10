/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of MBA.
 *
 * MBA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * MBA is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with MBA. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using MBA.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the MBA library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

//! \mainpage Multilevel B-splines Reference Manual
//!
//! \image html spockAndTerrain.jpg "From scattered data to B-spline surfaces at increasing level of detail"
//!
//! <p>
//! <a name="intro"><h2>Short introduction</h2></a>
//! This documentation contains a brief reference manual for the
//! <b>SINTEF</b> <b>Multilevel</b> <b>B-spline</b> <b>Library</b>
//! developed at <a
//! href="http://www.sintef.no/content/page3____342.aspx">SINTEF
//! Applied Mathematics</a>. The main interface is through the class
//! <a class="el" href="classMBA.html">MBA</a> which also contains a
//! small example in its documentation. Note that only the
//! non-adaptive/non-parametric version of the algorithms are
//! documented here. Parametric and adaptive versions, as used to
//! generate the head-model in the figure above, are implemented in
//! the classes MBApar (parametric), MBAadaptive (adaptive) and
//! MBAadaptivePar (adaptive and parametric).
//! <p>
//! For a thorough description of the basic schemes, the papers by the
//! originators of Multilevel B-splines, S. Lee, G. Wolberg and
//! S. Y. Shin, should be consulted.
//! <p>
//! A detailed description of the mathematical content of the
//! algorithms, with extensions made in the SINTEF library and with
//! numerical examples, can be found in the report:<br>
//! <em>Ø. Hjelle. Approximation of Scattered Data with Multilevel
//! B-splines. SINTEF report, 2001.</em><br>
//! <a href="http://www.sintef.no/upload/IKT/9011/geometri/MBA/MBA.pdf">Full report (780 K, pdf)</a></p>
//! <p>
//! <a name="gettingstarted"><h2>Getting started</h2></a>
//! It's very easy - just look at the small <a class="el"
//! href="classMBA.html#mba_example1">example</a> in class <a
//! class="el" href="classMBA.html">MBA</a> and read the documentation
//! for the functions that are called there. Then, copy the complete
//! <a href="mainsimplest.html#mainsimplest">main program</a> to a
//! file - compile it and run it! The program samples the resulting
//! spline surface and writes the result to a VRML-file. You can
//! download a free VRML-viewer from <a
//! href="http://www.km.kongsberg.com/sim">Kongsberg SIM</a>.
//! <p>
//! If you have got the Visual C++ workspace with all the source code,
//! then you can just build the application and run it with the
//! current main program.
//! <p>
//! <a name="examples"><h2>Numerical examples</h2></a>
//! I will make numerical examples and discuss pros and cons
//! thoroughly later. Consult the SINTEF report for mathematical
//! details.
//! <p>
//! <a name="furtherWork"><h2>Further work</h2></a>
//! Note that the basic algorithms in this library does not solve "the
//! scattered data approximation problem" in general. The algorithms
//! will produce anomalies near the data if the underlying grid of
//! B-spline coefficients is dense and the data points are unevenly
//! distributed in the domain. The remedy, as far as I have concluded
//! through numerical experiments, is to combine the basic schemes
//! with smoothing operators. This is now being implemented and will
//! be available in the next version.
//! <p>
//! <a name="download"><h2>Download</h2></a> 
//! A GPL-version of the library (for Linux/Unix) can be downloaded from <a
//! href="http://www.sintef.no/math_software">http://www.sintef.no/math_software</a>.
//! <p>
//! Please report any problems or comments to <a
//! href="mailto:jan.b.thomassen@sintef.no">jan.b.thomassen@sintef.no</a>.
//! <p>
//! Øyvind Hjelle, June 2001
//! <p>
//! <em>Last modified 28.11.2007 by <a
//! href="mailto:jan.b.thomassen@sintef.no">Jan Thomassen</a>
//!
//! \image html SINTEFlogo.gif
//!
//! \page lsmglib Least Squares Approximation of Scattered Data with B-splines
//! <a href="http://www.sintef.no/upload/IKT/9011/geometri/LSMG/lsmg_doc/index.html">Least Squares Approximation of Scattered Data with B-splines</a>
