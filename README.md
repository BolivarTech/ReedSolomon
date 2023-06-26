# Reed–Solomon Error-Correcting Codes

# Project Details

![OS](https://img.shields.io/badge/OS-JRE-darkgreen)
![licence](https://img.shields.io/badge/language-Java-brightgreen.svg?style=flat-square)
![license](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)

Implements the Reed–Solomon error correction protocol on Java

This project is provided "AS IS", no warranty, no responsibilities, no more documentation that the one included in this
repository.

## What is Reed–Solomon error correction

Reed–Solomon codes are a group of error-correcting codes that were introduced by Irving S. Reed and Gustave Solomon in 1960.

Reed–Solomon codes operate on a block of data treated as a set of finite-field elements called symbols. Reed–Solomon
codes are able to detect and correct multiple symbol errors. By adding t = n − k check symbols to the data, a Reed–Solomon
code can detect (but not correct) any combination of up to t erroneous symbols, or locate and correct up to ⌊t/2⌋ erroneous
symbols at unknown locations. As an erasure code, it can correct up to t erasures at locations that are known and provided
to the algorithm, or it can detect and correct combinations of errors and erasures. Reed–Solomon codes are also suitable
as multiple-burst bit-error correcting codes, since a sequence of b + 1 consecutive bit errors can affect at most two
symbols of size b. The choice of t is up to the designer of the code and may be selected within wide limits.

## Dependencies

This library depends on [BTUtils](https://github.com/BolivarTech/BTUtils) library.

## Credits

- [Julian Bolivar](https://www.linkedin.com/in/jbolivarg/)

## License

Copyright © [BolivarTech](https://www.bolivartech.com) 2022 All Rights Reserved.

ISO WIM Costumizer © 2022 by [Julian Bolivar](https://www.bolivartech.com) is licensed under [Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/legalcode)

Please see [License File](LICENSE.md) for more information.