# linpack_bench_s.ll
#
#  Modified:
#
#    24 October 2003
#
#  Author:
#
#    John Burkardt
#
#@ job_name = linpack_bench_s
#@ class = short
#@ output = linpack_bench_s.output
#@ job_type = serial
#@ queue
#
#  Compile and load the program.
#
xlf90 linpack_bench_s.f
mv a.out linpack_bench_s
#
#  Run the compiled program.
#
linpack_bench_s
rm linpack_bench_s
