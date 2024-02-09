// This work is licensed under the Creative Commons Attribution 3.0 Unported License.
// To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
// or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
// California, 94041, USA.

// Persistence Of Vision raytracer sample file.
// File by Dieter Bayer.
// This scene shows fog with transmittance used.
//
// -w320 -h240
// -w800 -h600 +a0.3

#version 3.6;

global_settings {
  assumed_gamma 1
}

#include "colors.inc"

camera {
   location  <0, 20, -100>
   direction <0,  0,    1>
   up        <0,  1,    0>
   right   x*image_width/image_height
}

background { colour SkyBlue }

fog{
    color rgbt <.7,.7,.7,.25>
    fog_type 2
    fog_alt 0.5
    fog_offset 0
    distance 1.5
    turbulence <.15, .15, .15>
    omega 0.35
    lambda 1.25
    octaves 5
}

// Put down the beloved famous raytrace green/yellow checkered floor
plane { y, -10
   pigment {
      checker colour Yellow colour Green
      scale 20
   }
   finish {
      ambient 0.2
      diffuse 0.8
   }
}

sphere { <0, 25, 0>, 40
   pigment {Red}
   finish {
      ambient 0.2
      diffuse 0.6
      phong 1.0
      phong_size 20
   }
}

sphere { <-100, 150, 200>,  20
   pigment {Green}
   finish {
      ambient 0.2
      diffuse 0.6
      phong 1.0
      phong_size 20
   }
}

sphere { <100, 25, 100>, 30
   pigment {Blue}
   finish {
      ambient 0.2
      diffuse 0.6
      phong 1.0
      phong_size 20
   }
}

light_source {<100, 120, 40> colour White}
