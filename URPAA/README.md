![Image](img/logo_urb_long.png)

URPAA: Urbauramon Psycho-acoustic Annoyance Analyzer
=============================================

### Table of Contents

**[Requirements](#requirements)**  
**[What can I do with URPAAnalyzer?](#what-can-i-do-with-urpaanalyzer)**   
**[How to use URPAAnalyzer?](#how-to-use-urpaanalyzer)**   
**[Credits and License](#credits-and-license)**    
**[Authors](#authors)**

Requirements
------------

* **Matlab R2018b** (or later version)
* **Signal Processing Toolbox** (Mathworks)
* **DSP System Toolbox** (Mathworks)


What can I do with URPAA?
----------------------------------

URPAAnalyzer allows to analyze the psycho-acoustic annoyance of an audio file supplied as input parameter. 
The nuisance model is based on the Zwicker model, described in:
Zwicker E., Fastl H. ‘Psychoacoustics: Facts and Models’(1990).

It should be noted that only single-channel with 1 second duration audio is accepted and that  
it is resampled at 16Khz before analysis.
But it is very simple to select channels using basic Matlab functions and create *.wav files 
with the necessary characteristics.


How to use URPAAnalyzer?
----------------------------------

Together with the URPAAnalizer function we have supplied 2 files with which you can perform example tests, one with a high annoyance value and one with a low value:
highPAsound.wav
lowPAsound.wav

If from the Matlab command line you are in the URPAA folder, you can call the URPAAnalizer function to analyze an audio file as follows:

	PA = URPAAnalyzer('lowPAsound.wav')

And you'll get:

	PA = 6.7710

In the same way you can analyze any audio file of 1 second duration.


Credits and License
----------------------------------

This is the source distribution of **URPAA: Urbauramon Psycho-acoustic Annoyance Analyzer** licensed
under the GPLv3+. Please consult the file COPYING for more information about
this license.

Website: https://github.com/jausegar/urbauramon/URPAA

If you have questions, bug reports or feature requests, please use the [Issue
Section on the website](https://github.com/jausegar/urbauramon/URPAA/issues) to report them. 

If you use our Audio Analyzer for your publications please cite our SEA papers:

"Zwicker’s Annoyance model implementation in a WASN node"
A. Pastor-Aparicio, J. Lopez-Ballester, J. Segura-Garcia, S. Felici-Castell, M. Cobos-Serrano, R. Fayos-Jordán, J.J. Pérez-Solano.
INTERNOISE 2019, 48th International Congress and Exhibition on Noise Control Engineering

"Visualization of nuisance information in acoustic environments using an IoT system"
J. Segura-Garcia, J. Lopez-Ballester, A. Pastor-Aparicio, S. Felici-Castell, M. Cobos-Serrano, J.J. Pérez-Solano, A. Soriano-Asensi, M. Garcia-Pineda.
INTERNOISE 2019, 48th International Congress and Exhibition on Noise Control Engineering


Copyright (c) 2018-2019   
Department of Computer Science    
Universitat de València  
Av. de la Universitat s/n, 46100, Valencia, Spain.  

Authors
------------
 
* Jaume Segura-Garcia
* Jesus Lopez-Ballester
* Adolfo Pastor-Aparicio