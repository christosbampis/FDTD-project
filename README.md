FDTD-project
============

Project in Electronic Compatibility course, using Matlab and openEMS software

  
 openEMS software can be found at http://openems.de/start/index.php

 Run movie.m using predefined parameters, postprocessing.m can also be used for displaying more results.
 
 Some results are accseible ffrom the link file where pictures for testcases are uploaded. See the video
 for a complete view of the phenomenon for different z-slices.

 The basic concepts are to estimate how the radiation patterns of an antenna (modelled as a dipole)
 are subject to change in presence of a human head. The shape of the human head is considered a sperical one,
 and its size can be modified whether we are studying an adult or a child case. The parameters under examination
 are SAR, input impedance Z_in and the reflection coefficient S11. Related ideas can be studied in the following
 links:
 
 http://en.wikipedia.org/wiki/Reflection_coefficient
 
 http://en.wikipedia.org/wiki/Input_impedance
 
 http://en.wikipedia.org/wiki/Impedance_matching
 
 Numerous parameters can be adjusted, like the dielectrical properties of the human model, the distance between
 the antenna and the human model, the layers in the human model and a small inhomegenous area can also be placed 
 inside the human head to represent a possible tumorous tissue.
 
 The coordinate system is Cartesian and the antenna's axis can be either x, y or z.
 
