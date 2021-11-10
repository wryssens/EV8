#-------------------------------------------------------------------------------
# Example script for EV8
#  1) Make compilation parameters for Zr84 in a modest box.
#  2) Compile Nil8
#  3) Run Nil8 to get a fort.13 with Nilson wavefunctions.
#  4) Move the output of Nilson to the input of EV8
#  5) Compile EV8 in the same box.
#  6) Run unconstrained EV8 to get the global minimum of the energy surface of
#     Zr84.
#  7) Use the output fort.13 as an input fort.12 for a constrained calculation
#     of Zr84.
#-------------------------------------------------------------------------------

#Compiler to use
compiler=gfortran-4.4

cp ../Codes/*.f ../Examples

echo '-----------------------------------------------'
echo ' Creating the compilation parameters.'

cat <<eof  >param8.h
      parameter (mx=16,my=16,mz=16,mc=10,mv=mx*my*mz,mq=4*mv,mw=130)
      parameter (meven=5,modd=5)
eof
echo ' Done!'
echo '-----------------------------------------------'
echo 'Compiling nil8!'

$compiler -O3 -o xnil8 nil8.f &> nil8.diag

echo 'Done!'
echo '-----------------------------------------------'
echo 'Running nil8!'
touch fort.13 zr84.nil8.out
\rm fort.13 zr84.nil8.out

 ./xnil8 <<EOF  > zr84.nil8.out
Zr84  Sly4  Nil8   
0060 0040 0044 0040
0000 0000 0000
Sly4
0016 0016 0016
0.80000000E+00
0000 00000
0.23300000E+00 1.050000000E+00
0.10000000E+03
EOF
echo 'Done!'
\rm xnil8
echo '-----------------------------------------------'
echo 'Compiling ev8!'

$compiler -O3 -o xev8 ev8.f &> ev8.diag
echo 'Done!'

echo '-----------------------------------------------'
echo 'Running unconstrained ev8.'

mv fort.13 fort.12

time ./xev8 <<EOF  > zr84.out
 Zr84 Sly4 Ev8									  
 0005 0005 0005 0002								  
 0.01000000E+00								        
 0500 0025 0000                                                    
 0100 0000                                                         
 0044 0040                                                          
Sly4                                                                                 
 0004 0001 0000                                                     
 1.25000000E+03 0.05000000E+02 0.00000000E+00 1.00000000E+00        
 1.25000000E+03 0.05000000E+02 0.00000000E+00    
 5.00000000E-08 3.00000000E-03 0.01000000E+00 1.00000000E-05        
 0001 0000 0000                                                     
 0.10000000E+00 0.02000000E+00 4.00000000E+00                        
 0.00000000E+00 0.00000000E+00                                      
 0.00000000E-04 1.00000000E+00                                      
 0.00000000E+02 0.00000000E+00                                      
EOF

echo 'Done!'
echo '-----------------------------------------------'
echo 'Running constrained ev8.'


mv fort.13 zr84.ev8
cp zr84.ev8 fort.12


time ./xev8 << EOF  > zr84.ev8.q200.out
 Zr84 Sly4 Ev8
 0000 0000 0000 0002
 0.01000000E+00
 0500 0025 0000
 0100 0000 
 0044 0040
Sly4
 0004 0001 0000 0000 0000
 1.25000000E+03 0.05000000E+02 0.00000000E+00 1.00000000E+00
 1.25000000E+03 0.05000000E+02 0.00000000E+00 
 5.00000000E-08 2.00000000E-03 0.01000000E+00 0.01000000E-03
 0001 0000 0000
 0.10000000E+00 0.02000000E+00 4.00000000E+00
 0.00000000E+00 0.00000000E+00
 2.00000000E-03 1.00000000E+02
 2.00000000E+00 0.00000000E+00
EOF

\mv fort.13 zr84.ev8.q200
rm fort.12
echo 'Done!'
echo '-----------------------------------------------'

\rm *.f
\rm xev8
