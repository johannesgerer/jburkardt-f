! 
!  for bounded flows
!
   integer ( i2b ), parameter :: c_b = 0
   integer ( i2b ), parameter :: b_n = 1
   integer ( i2b ), parameter :: b_s = 2
   integer ( i2b ), parameter :: b_w = 4
   integer ( i2b ), parameter :: b_e = 8
   integer ( i2b ), parameter :: b_nw = b_n + b_w
   integer ( i2b ), parameter :: b_sw = b_s + b_w
   integer ( i2b ), parameter :: b_ne = b_n + b_e
   integer ( i2b ), parameter :: b_se = b_s + b_e

   integer ( i2b ), parameter :: c_f = 16
   integer ( i2b ), parameter :: c_x = c_f - 1
   integer ( i2b ), parameter :: c_a = c_f + b_n + b_s + b_e + b_w
! 
!  for free-boundary flows
!
   integer ( i2b ), parameter :: c_o = 256
   integer ( i2b ), parameter :: c_w = 512
   integer ( i2b ), parameter :: c_s = 1024
   integer ( i2b ), parameter :: c_n = 2048

   integer ( i2b ), parameter :: c_e = 4096

   integer ( i2b ), parameter :: c_wo = c_w + c_o
   integer ( i2b ), parameter :: c_ns = c_n + c_s
   integer ( i2b ), parameter :: c_sw = c_s + c_w
   integer ( i2b ), parameter :: c_nw = c_n + c_w
   integer ( i2b ), parameter :: c_no = c_n + c_o
   integer ( i2b ), parameter :: c_so = c_s + c_o

   integer ( i2b ), parameter :: c_swo = c_s + c_w + c_o
   integer ( i2b ), parameter :: c_nsw = c_n + c_s + c_w
   integer ( i2b ), parameter :: c_nwo = c_n + c_w + c_o
   integer ( i2b ), parameter :: c_nso = c_n + c_s + c_o

   integer ( i2b ), parameter :: c_nswo = c_n + c_s + c_w + c_o
