%	Program to calculate loudness based on DIN 45631
%	Based on BASIC Program Published in J. Acoust. Soc. Jpn (E) 12, 1 (1991)
%	by E. Zwicker, H. Fastl, U. Widmann, K. Kurakata, S. Kuwano, and S. Namba
%
%  This file is part of DIN45631.m
%  Copyright 2003 Aaron Hastings, Ray W. Herrick Laboratories and Purdue University
%  This program is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
%
%
%	"Re-Author":	Aaron Hastings, Herrick Labs, Purdue University
%	Date Started: 29 October 00
%	Last Modified: 29 November 01
%	Status: Program Correctly Calculates Loudness for a 70 dB 1000 Hz sine 
%			  filtered using 1/3 octave band filters
%
%	Syntax:
%	[N, Ns, err]=DIN45631(LT, MS)
%
%	This is a loudness function which:
%		Calculates loudness based on DIN 45631 / ISO 532 B (Zwicker)
%		Accepts 1/3 octave band levels (SPL Linear Weighting)
%		* This data must be calibrated using a separate calibration function
%
%	Input Variables
%	LT(28)			Field of 28 elements which represent the 1/3 OB levels in dB with 
%						fc = 25 Hz to 12.5 kHz
%	MS					String variable to distinguish the type of sound field ( free / diffues )
%	
%	Output Variables
%	N					Loudness in sone G
%	NS					Specific Loudness
%	err				Error Code

%	Working Variables
%	Ltbew 			Adjusted 1/3 octave band intensity for fc <= 315
%	FR(28)			Center frequencies of 1/3 OB
%	RAP(8)			Ranges of 1/3 OB levels for correction at low frequencies according 
%						to equal loudness contours
%	DLL(11,8)		Reduction of 1/3 OB levels at low frequencies according to equal 
%						loudness contours within the 8 ranges defined by RAP
%	LTQ(20)			Critical Band Rate level at absolute threshold without taking into 
%						account the transmission characterisics of the ear
%	AO(20)			Correction of levels according to the transmission characteristics 
%						of the ear
%	DDF(20)			Level difference between free and diffuse sound fields
%	DCB(20)			Adaptation of 1/3 OB levels to the corresponding critical band level
%	ZUP(21)			Upper limits of approximated critical bands in terms of critical 
%						band rate
%	RNS(18)			Range of specific loudness for the determination of the steepness of 
%						the upper slopes in the specific loudness -critical band rate pattern
%	USL(18,8)		Steepness of the upper slopes in the specific loudness - critical 
%						band rate pattern for the ranges RNS as a function of the number of 
%						the critical band
%	
%	Working Variables (Uncertain of Definitions)
%	XP					Equal Loudness Contours
%	TI					Intensity of LT
%	LCB				Lower Critical Band
%	LE					Level Excitation 
%	NM					Critical Band Level 
%	KORRY			Correction Factor
%	N					Loudness (in sones)
%	DZ					Speparation in CBR 
%	N2					Main Loudness 
%	Z1					Critical Band Rate for Lower Limit 
%	N1					Loudness of previous band 
%	IZ					Center "Frequency" Counter, used with NS
%	Z					Critical band rate 
%	J,IG				Counters used with USL 
%
%   err                 -1   ->    The level of the sound exceeded the range for which the metric is defined
%
%********************************************************************************
% Copyright (c) 2018-2019 Department of Computer Science                        *
%                         ETSE, Universitat de València                         *
%                         46100, Burjassot, Valencia, Spain                     *
%                                                                               *
% This file is part of the URPAA: URBAURAMON Psycho-acoustic Annoyance Analizer *
%                                                                               *
% URPAA is free software:  you can redistribute it and/or modify it  under      *
% the terms of the  GNU  General  Public  License  as published by the  Free    *
% Software Foundation, either version 3 of the License,  or (at your option)    *
% any later version.                                                            *
%                                                                               *
% URPAA  is distributed in the hope that it will be useful, but WITHOUT ANY     *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS    *
% FOR A PARTICULAR PURPOSE.                                                     *
% See the GNU General Public License for more details.                          *
%                                                                               *
% You should  have received a copy  of the GNU General Public License  along    *
% with this program.  If not, see <http://www.gnu.org/licenses/>.               *
%                                                                               *
% URPAAnalyzer is a function to perform psycho-acoustic annoyance analisys of   *
% an audio input file in *.wav format.                                          *
% The psycho-acoustic model used to perform the analysis is based on Zwicker's  *
% psycho-acoustic annoyance model.												*
%                                                                               *
% In case of use of the implemented code, please cite us.                       *
%                                                                               *
% https://github.com/jausegar/urbauramon/tree/master/URPAA  jaume.segura@uv.es  *
%********************************************************************************

function [N, NS, err]=DIN45631(LT, MS)

%%	Begin initializing the working variables

FR=[25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1.0 1.25 1.6 ...
      2.0 2.5 3.15 4.0 5.0 6.3 8.0 10.0 12.5];

RAP=[45 55 65 71 80 90 100 120];

DLL=[-32 -24 -16 -10 -5 0 -7 -3 0 -2 0
   -29 -22 -15 -10 -4 0 -7 -2 0 -2 0
   -27 -19 -14 -9  -4 0 -6 -2 0 -2 0
   -25 -17 -12 -9  -3 0 -5 -2 0 -2 0
   -23 -16 -11 -7  -3 0 -4 -1 0 -1 0
   -20 -14 -10 -6  -3 0 -4 -1 0 -1 0
   -18 -12 -9  -6  -2 0 -3 -1 0 -1 0
   -15 -10 -8  -4  -2 0 -3 -1 0 -1 0]';	%%	BASIC code does this oddly, hence the transpose

LTQ=[30 18 12 8 7 6 5 4 3 3 3 3 3 3 3 3 3 3 3 3];   

AO=[0 0 0 0 0 0 0 0 0 0 -0.5 -1.6 -3.2 -5.4 -5.6 -4.0 -1.5 2.0 5.0 12.0];

DDF=[0 0 0.5 0.9 1.2 1.6 2.3 2.8 3.0 2.0 0.0 -1.4 -2.0 -1.9 -1.0 0.5 ... 
      3.0 4.0 4.3 4.0];

DCB=[-0.25 -0.6 -0.8 -0.8 -0.5 0.0 0.5 1.1 1.5 1.7 1.8 1.8 1.7 1.6 1.4 ...
      1.2 0.8 0.5 0.0 -0.5];

ZUP=[0.9 1.8 2.8 3.5 4.4 5.4 6.6 7.9 9.2 10.6 12.3 13.8 15.2 16.7 18.1 ... 
      19.3 20.6 21.8 22.7 23.6 24.0];

RNS=[21.5 18.0 15.1 11.5 9.0 6.1 4.4 3.1 2.13 1.36 0.82 0.42 0.30 0.22 ... 
      0.15 0.10 0.035 0.0];

USL=[13.00 8.20 6.30 5.50 5.50 5.50 5.50 5.50
   9.00 7.50 6.00 5.10 4.50 4.50 4.50 4.50
   7.80 6.70 5.60 4.90 4.40 3.90 3.90 3.90
   6.20 5.40 4.60 4.00 3.50 3.20 3.20 3.20
   4.50 3.80 3.60 3.20 2.90 2.70 2.70 2.70
   3.70 3.00 2.80 2.35 2.20 2.20 2.20 2.20
   2.90 2.30 2.10 1.90 1.80 1.70 1.70 1.70
   2.40 1.70 1.50 1.35 1.30 1.30 1.30 1.30
   1.95 1.45 1.30 1.15 1.10 1.10 1.10 1.10
   1.50 1.20 0.94 0.86 0.82 0.82 0.82 0.82
   0.72 0.67 0.64 0.63 0.62 0.62 0.62 0.62
   0.59 0.53 0.51 0.50 0.42 0.42 0.42 0.42
   0.40 0.33 0.26 0.24 0.22 0.22 0.22 0.22
   0.27 0.21 0.20 0.18 0.17 0.17 0.17 0.17
   0.16 0.15 0.14 0.12 0.11 0.11 0.11 0.11
   0.12 0.11 0.10 0.08 0.08 0.08 0.08 0.08
   0.09 0.08 0.07 0.06 0.06 0.06 0.06 0.05
   0.06 0.05 0.03 0.02 0.02 0.02 0.02 0.02];

%%	Begin Loudness Calcultation

%%	Correction of 1/3 OB levels according to equal loudness contours (XP) and 
%%	calculation of the intensities for 1/3 OB's up to 315 Hz

for I=1:11
   J=1;
   while J<8  %%  ALH I think that this should be changed from J<8 to J<=8  (it never gets to RAP(8))
      if LT(I)<=RAP(J)-DLL(I,J)
         XP=LT(I)+DLL(I,J);
         TI(I)=10^(0.1*XP);
         J=9;	%%	To exit from the while loop
      else
         J=J+1;
      end
   end
end


%%	Determination of Levels LCB(1), LCB(2), and LCB(3) within the first three
%%	critical bands

if length(TI)<11
    disp(' ')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    disp('WARNING!!!! ERROR!!!!! One of the 1/3 OBs level exceded the metrics range.')
    disp('If you are using this code to normalize loudness, reduce the initial level and try again.')
    disp('Otherwise any results if obtained will be incorrect.')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    disp(' ')
    N=-1;
    NS=-1*ones(1,240);
    err=-1;
    return
end  %%  ALH 7 July 03

GI(1)=TI(1)+TI(2)+TI(3)+TI(4)+TI(5)+TI(6);
GI(2)=TI(7)+TI(8)+TI(9);
GI(3)=TI(10)+TI(11);

for I=1:3
   if GI(I)>0
      LCB(I)=10*log10(GI(I));
      %%	Note: The BASIC code uses a divide by "log(10)" to gauruntee that
      %%	the log is base 10
   end
end

%%	Calculation of Main Loudness

for I=1:20
   LE(I)=LT(I+8);
   if I<=3  %%  Corrected 26 Nov 01 from "+" to "=" 
      LE(I)=LCB(I);
   end
   LE(I)=LE(I)-AO(I);
   NM(I)=0;
   if MS=='d' | MS=='D'
      LE(I)=LE(I)+DDF(I);
   end
   if LE(I)>LTQ(I)
      LE(I)=LE(I)-DCB(I);
      S=0.25;
      MP1=0.0635*10^(0.025*LTQ(I));
      MP2=(1-S+S*10^(0.1*(LE(I)-LTQ(I))))^0.25-1;
      NM(I)=MP1*MP2;
      if NM(I)<=0
         NM(I)=0;
      end
   end
end

NM(21)=0;

%%	Correction of specific loudness in the lowest critical band taking 
%%	into account the dependence of absolute threshold within this critical band      

KORRY=0.4+0.32*NM(1)^0.2;
if KORRY>1
   KORRY=1;
end
NM(1)=NM(1)*KORRY;

%%	Start Values

N=0;
Z1=0;
N1=0;
IZ=1;
Z=0.1;
short=0;

%%	Step to first and subsequent critical bands

for I=1:21
   ZUP(I)=ZUP(I)+0.0001;
   IG=I-1;
   if IG>8
      IG=8;
   end
   while Z1<ZUP(I)	%%	Note, Z1 will always be < ZUP(I) when line is first reached for each I
      if N1>NM(I)	
         %%	Decide whether the critical band in question is completely or
         %%	partly masked by accessory loudness
         N2=RNS(J);
         if N2<NM(I)
            N2=NM(I);
         end
         DZ=(N1-N2)/USL(J,IG);
         Z2=Z1+DZ;
         if Z2>ZUP(I)
            Z2=ZUP(I);
            DZ=Z2-Z1;
            N2=N1-DZ*USL(J,IG);
         end
         %%	Contribtion of accessory loudness to total loudness
         N=N+DZ*(N1+N2)/2;
         while Z<Z2
            NS(IZ)=N1-(Z-Z1)*USL(J,IG);
            IZ=IZ+1;
            Z=Z+0.1;
         end
      elseif N1==NM(I)
         %%	Contribution of umasked main loudness to total loudness and calculation
         %%	of values NS(IZ) with a spacing of Z=IZ*0.1 Bark
         Z2=ZUP(I);
         N2=NM(I);
         N=N+N2*(Z2-Z1);
         while Z<Z2
            NS(IZ)=N2;
            IZ=IZ+1;
            Z=Z+0.1;
         end
      else
         %%	Determination of the number J corresponding to the range of specific
         %%	loudness
         for J=1:18
            if RNS(J)<NM(I)
               break
            end
         end
         %%	Contribution of umasked main loudness to total loudness and calculation
         %%	of values NS(IZ) with a spacing of Z=IZ*0.1 Bark
         Z2=ZUP(I);
         N2=NM(I);
         N=N+N2*(Z2-Z1);
         while Z<Z2
            NS(IZ)=N2;
            IZ=IZ+1;
            Z=Z+0.1;
         end
      end
      %%	Step to next segment
      while J<18
         if N2<=RNS(J)
            J=J+1;
         else
            break
         end
      end
      if N2<=RNS(J) & J>=18
         J=18;
      end
      Z1=Z2;
      N1=N2;
   end
end

%%	Now apply rounding types of corrections

if N<0
   N=0;
elseif N<=16
   N=floor(N*1000+0.5)/1000;
else
   N=floor(N*100+0.5)/100;
end

err=0;





