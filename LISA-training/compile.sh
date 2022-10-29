echo '-----------------------------------------------'
echo ' Creating the compilation parameters.'

cat <<eof  >param8.h
      parameter (mx=10,my=10,mz=14,mc=10,mv=mx*my*mz,mq=4*mv,mw=220)
      parameter (meven=5,modd=5)
eof
echo ' Done!'
echo '-----------------------------------------------'
echo 'Compiling nil8!'
echo '-----------------------------------------------'
gfortran -O3 -o nil8.exe nil8.f &> nil8.diag
echo '-----------------------------------------------'
echo 'Compiling ev8!'
echo '-----------------------------------------------'
gfortran -O3 -o ev8.exe ev8.f&> ev8.diag
echo 'Done!'

