ZERO
===

JPEG Grid Detection based on the Number of DCT Zeros
and its Application to Automatic and Localized Forgery Detection

================================================

Version 1 - May 13th, 2019

by Tina Nikoukhah <tina.nikoukhah@cmla.ens-cachan.fr>

and Rafael Grompone von Gioi <grompone@cmla.ens-cachan.fr>

joint work with Jérémy Anger, Thibaud Ehret, Miguel Colom and Jean-Michel Morel


Introduction
------------

ZERO is an implementation of the JPEG grid detector applied to forgery
detection in images described in the paper:

     "JPEG Grid Detection based on the Number of DCT Zeros and its
     Application to Automatic and Localized Forgery Detection" by Tina
     Nikoukhah, Jérémy Anger, Thibaud Ehret, Miguel Colom, Jean-Michel
     Morel and Rafael Grompone von Gioi.
[PDF](http://openaccess.thecvf.com/content_CVPRW_2019/papers/Media%20Forensics/Nikoukhah_JPEG_Grid_Detection_based_on_the_Number_of_DCT_Zeros_CVPRW_2019_paper.pdf)     


Online Demo
------------

[IPOL](https://ipolcore.ipol.im/demo/clientApp/demo.html?id=77777000073)

Files
-----

- zero.c: Main code.

- README.txt: this file.

- LICENSE: GNU AFFERO GENERAL PUBLIC LICENSE Version 3.

- Makefile: Compilation instructions.

- iio.{c,h}: [iio](https://github.com/mnhrdt/iio) code and header.

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


To verify a correct compilation you can apply the algorithm to the
test images. This can be done by executing:
```bash
    make test
```

This should print the following message:
```
test on roma.pgm
----------------
./zero roma.pgm
no overall JPEG grid found

test on pelican.ppm
-------------------
./zero pelican.ppm
main grid: #6 [6 0] log(nfa) = -6821.13

test on tampered1.pgm
---------------------
./zero tampered1.pgm
no overall JPEG grid found
forgery found: 104 94 - 153 159 [50x66] grid: #0 [0 0] n 68 k 30 log(nfa) = -25.2505

test on tampered2.ppm
---------------------
./zero tampered2.ppm
main grid: #6 [6 0] log(nfa) = -6618.52
forgery found: 330 68 - 401 104 [72x37] grid: #34 [2 4] n 81 k 36 log(nfa) = -32.2729
```


Copyright and License
---------------------

Copyright (c) 2018-2019 Rafael Grompone von Gioi <grompone@gmail.com>
Copyright (c) 2018-2019 Tina Nikoukhah <nikoukhah@cmla.ens-cachan.fr>
Copyright (c) 2018-2019 Jérémy Anger <anger@cmla.ens-cachan.fr>
Copyright (c) 2018-2019 Thibaud Ehret <ehret@cmla.ens-cachan.fr>

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
