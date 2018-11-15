module load intel impi fftw3 

cd ..
make clean
./configure --enable-fft --with-particle_self_gravity=fft_disk --enable-mpi --enable-shearing-box --enable-fargo --with-particles=feedback --with-gas=hydro --with-eos=isothermal --with-problem=par_strat3d_turb --with-order=3p 
make all MACHINE=stampede2
cd bin

