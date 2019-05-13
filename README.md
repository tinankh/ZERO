ZERO
===

JPEG Grid Detection based on the Number of DCT Zeros
and its Application to Automatic and Localized Forgery Detection

================================================

Version 1 - May 13th, 2019

by Tina Nikoukhah <tina.nikoukhah@gmail.com>

and Rafael Grompone von Gioi <grompone@gmail.com>

joint work with Jérémy Anger, Thibaud Ehret, Miguel Colom and Jean-Michel Morel


Introduction
------------

ZERO is an implementation of the JPEG grid detector applied to forgery
detection in images described in the paper:

     "JPEG Grid Detection based on the Number of DCT Zeros and its
     Application to Automatic and Localized Forgery Detection" by Tina
     Nikoukhah, Jérémy Anger, Thibaud Ehret, Miguel Colom, Jean-Michel
     Morel and Rafael Grompone von Gioi.


Files
-----

- README.txt: this file.

- LICENSE: GNU AFFERO GENERAL PUBLIC LICENSE Version 3.

- Makefile: Compilation instructions.

- iio.{c,h}: [iio](https://github.com/mnhrdt/iio) code and header.

- zero.c: Main code.

- *.{ppm,pgm}: Test images.


Compiling
---------
The compiling instruction is just
```bash
    make
```
or if you want the code to be parallel
```bash
   make openmp
```
from the directory where the source codes and the Makefile are located.


Running ZERO Command
-------------------
The command execution is just
```bash
    ./zero
```

Copyright and License
---------------------

Copyright (c) 2018-2019 Rafael Grompone von Gioi <grompone@gmail.com>

ZERO is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

ZERO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.


Thanks
------

I would be grateful to receive any comment, especially about errors,
bugs, or strange results.
