#-------------------------------------------------------------------------------
# Example script for EV8.
#
# This script executes the following steps:
#
# 1) Compile nil8 & ev8 for computations for Pb208.
# 2) Run nil8 to create an initial guess on a file fort.12
# 3) Run ev8 with fort.12 to get a (bad) result for Pb208.
# 4) Recompile ev8.
# 5) Rerun the calculation in a bigger box.
#-------------------------------------------------------------------------------

#Compiler to use
compiler=gfortran

cp ../Codes/*.f ../Examples

echo '-----------------------------------------------'
echo ' Creating the compilation parameters.'
cat << eof > param8.h
      parameter (mx=10,my=10,mz=10,mc=10,mv=mx*my*mz,mq=4*mv,mw=160)
      parameter (meven=5,modd=5)
eof


echo ' Done!'
echo '-----------------------------------------------'
echo 'Compiling nil8!'
#Compiling NIL8
$compiler -o nil8.exe nil8.f &> nil8.diag
echo 'Compiled!'
echo 'Done!'
echo '-----------------------------------------------'
# Making a runtime data file for nil8
echo 'Creating runtime NIL8 parameters.'
cat << eof > nil8.data
Pb208 Sly4 Nil8   
0080 0080 0126 0082
0000 0000 0000
Sly4
 0010 0010 0010
 1.0       E+00
 0000 0000
 0.16      E+00 1.05      E+00
 0.100     E+03
eof
echo 'Done!'
echo '-----------------------------------------------'
# Running nil8
echo 'Running NIL8' 
./nil8.exe < nil8.data > Pb208.nil8.out
echo 'Done!'

# Cleaning up
rm nil8.data
rm nil8.exe

# Preparing fort.12 for EV8 run
mv fort.13 Pb208.nil8.wf
cp Pb208.nil8.wf fort.12

# Compiling EV8
echo '-----------------------------------------------'
echo 'Compiling EV8'
$compiler -O3 -o ev8.exe ev8.f &> ev8.diag
echo 'Done!'
echo '-----------------------------------------------'
# Making a runtime data file for ev8
cat << eof > ev8.data
 Pb208 Sly4 Ev8
 0010 0010 0010 0000
 0.015     E+00
 0500 0000 0000
 0100 0000
 0126 0082
Sly4
 0005 0001
 1.25000000E+03 5.00000000E+00 0.00000000E+00 1.00000000E+00
 1.25000000E+03 5.00000000E+00 0.00000000E+00 
 1.00000000E-09 1.00000000E-07 1.00000000E-07 1.00000000E-07
 0000 0000 0000 0000
 0.00010000E+00 0.00000000E+00 4.00000000E+00
 0.00000000E+00 0.00000000E+00
 0.00000000E+00 1.00000000E+00
 0.00000000E+00 0.00000000E+00
eof

# Running EV8
echo 'Running EV8; in the small box.'
./ev8.exe < ev8.data > Pb208.ev8.Sly4.SmallBox.out
echo 'Done!'
echo '-----------------------------------------------'
# Cleaning up
mv fort.13 Pb208.ev8.Sly4.Smallbox.wf
rm fort.12

#Making a new param8.h file for the bigger box
echo 'Making a param8.h for a bigger box!'
cat << eof > param8.h
      parameter (mx=15,my=15,mz=15,mc=15,mv=mx*my*mz,mq=4*mv,mw=160)
      parameter (meven=5,modd=5)
eof
echo 'Done!'
echo '-----------------------------------------------'
#Compiling size8
echo 'Compiling ev8 in a bigger box.'
$compiler -o ev8.exe ev8.f &> ev8.diag
echo 'Done'
echo '-----------------------------------------------'

# Making a runtime data file for ev8
cat << eof > ev8.data
 Pb208 Sly4 Ev8
 0015 0015 0015 0000
 0.015     E+00
 0500 0000 0000
 0100 0000
 0126 0082
Sly4
 0005 0001
 1.25000000E+03 5.00000000E+00 0.00000000E+00 1.00000000E+00
 1.25000000E+03 5.00000000E+00 0.00000000E+00 
 1.00000000E-09 1.00000000E-07 1.00000000E-07 1.00000000E-07
 0000 0000 0000 0000
 0.00010000E+00 0.00000000E+00 4.00000000E+00
 0.00000000E+00 0.00000000E+00
 0.00000000E+00 1.00000000E+00
 0.00000000E+00 0.00000000E+00
eof

cp Pb208.ev8.Sly4.Smallbox.wf fort.12

# Running EV8
echo 'Running EV8; in the big box.'
./ev8.exe < ev8.data > Pb208.ev8.Sly4.BigBox.out
# Cleaning up
mv fort.13 Pb208.ev8.Sly4.Bigbox.wf
rm fort.12
echo 'Done!'
echo '-----------------------------------------------'
rm ev8.exe
rm *.f
