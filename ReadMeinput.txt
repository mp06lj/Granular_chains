**********************************************************
		The Taper Chain Program  

		Version 1-31-14 

		Version 3-08-14 : corrected file read issues for cygwin

**********************************************************
	Welcome to the Readme File
	I will first spray everything I read in.
	Then I will decode this and tell you what I found.
	Then I will begin the simulations and keep some output here.
**********************************************************
Getting all my setup from file : input.txt

Read line 1: c: comment line, this is ignored
Read line 2: c: not a bad idea to have some info here for your own notes
Read line 3: c: we are outputting all files, extension.dat
Read line 4: outputfiles: fxbvkp
Read line 5: output file extension: .dat
Read line 6: outputinterval: 5.0E-6
Read line 7: outputprecision: 6
Read line 8: output grains:1-100
Read line 9: fileName: monoDisperseTest
Read line 10: DEFAULTGRAIN: sphere#0.006|rho#7800.0|poisson#0.28|youngs#200.0E+9|
Read line 11: ChainPattern:[a#100]
Read line 12: subunit:a
Read line 13: loadstep: number#1|solver:velocityverlet|start#0|step#10.0E-12|stop#0.01|
Read line 14: rightwall: loadstep#1|subunit:a|radius:inf
Read line 15: leftwall: loadstep#1|subunit:a|radius:inf
Read line 16: initialcondition: loadstep#1|variable:velocity|addvalue#0.099188|grains:1
Read line 17: end
Read line 18: 
Read line 19: 
*******************************************
So I have read your input file and now find the following:
*******************************************
You want outputs to have the base name : monoDisperseTest
You want these output files: 
	Force output :  yes 
	Relative Location output :  yes 
	Absolute Location output :  yes 
	Velocity output :  yes 
	Potential Energy output :  yes 
	Kinetic Energy output :  yes 
	Overlaps output :  no 
	Interface forces output :  no 
And you want these outputs for the following grains (if they exist): 
	1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100

These were all the subunits you defined: 
	Unit Name: DEFAULTGRAIN

Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


	Unit Name: a

Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2



I read in a pattern which I decoded to be: 
	a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|a|

So that all in all these are the grains I have added to the chain: 
From left to right (going in the positive direction)

Grain 1
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 2
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 3
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 4
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 5
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 6
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 7
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 8
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 9
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 10
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 11
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 12
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 13
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 14
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 15
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 16
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 17
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 18
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 19
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 20
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 21
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 22
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 23
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 24
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 25
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 26
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 27
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 28
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 29
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 30
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 31
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 32
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 33
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 34
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 35
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 36
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 37
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 38
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 39
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 40
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 41
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 42
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 43
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 44
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 45
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 46
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 47
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 48
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 49
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 50
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 51
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 52
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 53
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 54
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 55
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 56
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 57
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 58
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 59
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 60
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 61
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 62
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 63
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 64
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 65
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 66
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 67
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 68
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 69
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 70
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 71
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 72
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 73
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 74
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 75
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 76
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 77
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 78
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 79
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 80
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 81
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 82
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 83
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 84
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 85
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 86
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 87
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 88
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 89
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 90
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 91
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 92
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 93
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 94
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 95
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 96
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 97
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 98
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 99
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2


Grain 100
	
Grain Youngs Modulus : 2e+11
Grain Poisson Ratio : 0.28
Grain Density : 7800
Grain Volume : 9.04779e-07
Grain mass : 0.00705727
Grain Shape : Shape Parameter
	 Spherical Grain R : 0.006
	 m : 1
	 alpha : 2





Now let's look at the interfaces between grains.
	Interface: (1,2)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (2,3)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (3,4)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (4,5)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (5,6)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (6,7)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (7,8)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (8,9)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (9,10)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (10,11)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (11,12)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (12,13)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (13,14)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (14,15)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (15,16)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (16,17)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (17,18)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (18,19)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (19,20)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (20,21)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (21,22)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (22,23)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (23,24)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (24,25)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (25,26)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (26,27)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (27,28)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (28,29)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (29,30)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (30,31)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (31,32)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (32,33)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (33,34)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (34,35)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (35,36)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (36,37)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (37,38)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (38,39)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (39,40)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (40,41)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (41,42)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (42,43)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (43,44)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (44,45)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (45,46)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (46,47)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (47,48)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (48,49)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (49,50)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (50,51)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (51,52)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (52,53)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (53,54)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (54,55)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (55,56)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (56,57)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (57,58)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (58,59)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (59,60)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (60,61)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (61,62)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (62,63)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (63,64)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (64,65)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (65,66)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (66,67)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (67,68)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (68,69)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (69,70)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (70,71)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (71,72)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (72,73)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (73,74)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (74,75)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (75,76)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (76,77)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (77,78)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (78,79)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (79,80)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (80,81)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (81,82)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (82,83)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (83,84)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (84,85)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (85,86)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (86,87)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (87,88)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (88,89)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (89,90)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (90,91)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (91,92)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (92,93)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (93,94)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (94,95)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (95,96)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (96,97)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (97,98)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (98,99)	 has a force prefactor :1.24473e+10 and force power law : 1.5
	Interface: (99,100)	 has a force prefactor :1.24473e+10 and force power law : 1.5
Wall interfaces (if used) will be treated as a special case
The undeformed length of the chain before any load steps is : 1.2
Details about this LoadStep  : 
	 start time : 0
	 stop time : 0.01
	 time step: 1e-11
	 body acceleration : 0
	 velocity based damping  : 0
	 force restitution (unload force is reduced by this factor)  : 1
Wall info: 
	 LEFT Wall : on
	 LEFT Wall Offset: 0
	 LEFT Wall Force Prefactor: 1.76032e+10
	 LEFT Wall Force Power: 1.5
	 RIGHT Wall : on
	 RIGHT Wall Offset: 0
	 RIGHT Wall Force Prefactor: 1.76032e+10
	 RIGHT Wall Force Power: 1.5
The following grains are movable ? 
Grain : 1 is movable : yes
Grain : 2 is movable : yes
Grain : 3 is movable : yes
Grain : 4 is movable : yes
Grain : 5 is movable : yes
Grain : 6 is movable : yes
Grain : 7 is movable : yes
Grain : 8 is movable : yes
Grain : 9 is movable : yes
Grain : 10 is movable : yes
Grain : 11 is movable : yes
Grain : 12 is movable : yes
Grain : 13 is movable : yes
Grain : 14 is movable : yes
Grain : 15 is movable : yes
Grain : 16 is movable : yes
Grain : 17 is movable : yes
Grain : 18 is movable : yes
Grain : 19 is movable : yes
Grain : 20 is movable : yes
Grain : 21 is movable : yes
Grain : 22 is movable : yes
Grain : 23 is movable : yes
Grain : 24 is movable : yes
Grain : 25 is movable : yes
Grain : 26 is movable : yes
Grain : 27 is movable : yes
Grain : 28 is movable : yes
Grain : 29 is movable : yes
Grain : 30 is movable : yes
Grain : 31 is movable : yes
Grain : 32 is movable : yes
Grain : 33 is movable : yes
Grain : 34 is movable : yes
Grain : 35 is movable : yes
Grain : 36 is movable : yes
Grain : 37 is movable : yes
Grain : 38 is movable : yes
Grain : 39 is movable : yes
Grain : 40 is movable : yes
Grain : 41 is movable : yes
Grain : 42 is movable : yes
Grain : 43 is movable : yes
Grain : 44 is movable : yes
Grain : 45 is movable : yes
Grain : 46 is movable : yes
Grain : 47 is movable : yes
Grain : 48 is movable : yes
Grain : 49 is movable : yes
Grain : 50 is movable : yes
Grain : 51 is movable : yes
Grain : 52 is movable : yes
Grain : 53 is movable : yes
Grain : 54 is movable : yes
Grain : 55 is movable : yes
Grain : 56 is movable : yes
Grain : 57 is movable : yes
Grain : 58 is movable : yes
Grain : 59 is movable : yes
Grain : 60 is movable : yes
Grain : 61 is movable : yes
Grain : 62 is movable : yes
Grain : 63 is movable : yes
Grain : 64 is movable : yes
Grain : 65 is movable : yes
Grain : 66 is movable : yes
Grain : 67 is movable : yes
Grain : 68 is movable : yes
Grain : 69 is movable : yes
Grain : 70 is movable : yes
Grain : 71 is movable : yes
Grain : 72 is movable : yes
Grain : 73 is movable : yes
Grain : 74 is movable : yes
Grain : 75 is movable : yes
Grain : 76 is movable : yes
Grain : 77 is movable : yes
Grain : 78 is movable : yes
Grain : 79 is movable : yes
Grain : 80 is movable : yes
Grain : 81 is movable : yes
Grain : 82 is movable : yes
Grain : 83 is movable : yes
Grain : 84 is movable : yes
Grain : 85 is movable : yes
Grain : 86 is movable : yes
Grain : 87 is movable : yes
Grain : 88 is movable : yes
Grain : 89 is movable : yes
Grain : 90 is movable : yes
Grain : 91 is movable : yes
Grain : 92 is movable : yes
Grain : 93 is movable : yes
Grain : 94 is movable : yes
Grain : 95 is movable : yes
Grain : 96 is movable : yes
Grain : 97 is movable : yes
Grain : 98 is movable : yes
Grain : 99 is movable : yes
Grain : 100 is movable : yes
At the beginning of the load step I will add this vector to the velocity:
0.099188	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
At the beginning of the load step I will add this vector to the displacement:
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Details regarding the Monte Carlo Params (if used): 
acceptance ratio : 0.4
max iterations : 100
max cutbacks for scale : 10
MC avgs : 1500
tolerance  : 0.001
random step scale : 0.1
random step shift : 0.25




