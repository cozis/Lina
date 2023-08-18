gcc bench_dot.c ../src/lina.c -O3 -march=native -ffast-math -funroll-loops -o bench_dot
./bench_dot
python3 py_dot.py