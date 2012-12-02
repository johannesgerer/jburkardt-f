!     Ercolessi-Adams glue potential for Al.
!     Ref.: F. Ercolessi and J. B. Adams, Europhys. Lett. 26, 583 (1994).
!

subroutine v2 ( arg, func, dfunc, d2func )

!*******************************************************************************
!
!! V2 stores data for the Aluminum pair potential.
!
!        Aluminum  : pair potential   and its first two derivatives.
!        Generated automatically by PoCo, version 04-may-93           
!        Hamiltonian type #  2, run on 93/06/09 at 15.04.43
!        Uses subroutine seval from netlib@ornl.gov [to get it,
!        use 'send seval from sfmm'], trivially modified to
!        compute also dfunc and d2func and use double precision.
!
!  Modified:
!
!    28 October 2005
!
  implicit real ( kind = 8 ) (a-h,o-z)

  parameter (nv2= 17)
  parameter (argmax=   0.555805441821810D+01)
  real ( kind = 8 ) xv2(nv2),yv2(nv2),bv2(nv2),cv2(nv2),dv2(nv2)
  save xv2,yv2,bv2,cv2,dv2
  data xv2(  1) /   0.202111069753385D+01 /
  data xv2(  2) /   0.227374953472558D+01 /
  data xv2(  3) /   0.252638837191732D+01 /
  data xv2(  4) /   0.277902720910905D+01 /
  data xv2(  5) /   0.303166604630078D+01 /
  data xv2(  6) /   0.328430488349251D+01 /
  data xv2(  7) /   0.353694372068424D+01 /
  data xv2(  8) /   0.378958255787597D+01 /
  data xv2(  9) /   0.404222139506771D+01 /
  data xv2( 10) /   0.429486023225944D+01 /
  data xv2( 11) /   0.454749906945117D+01 /
  data xv2( 12) /   0.480013790664290D+01 /
  data xv2( 13) /   0.505277674383463D+01 /
  data xv2( 14) /   0.530541558102636D+01 /
  data xv2( 15) /   0.555805441821810D+01 /
  data xv2( 16) /   0.555807968210182D+01 /
  data xv2( 17) /   0.555810494598553D+01 /
  data yv2(  1) /   0.196016472197158D+01 /
  data yv2(  2) /   0.682724240745344D+00 /
  data yv2(  3) /   0.147370824539188D+00 /
  data yv2(  4) /  -0.188188235860390D-01 /
  data yv2(  5) /  -0.576011902692490D-01 /
  data yv2(  6) /  -0.519846499644276D-01 /
  data yv2(  7) /  -0.376352484845919D-01 /
  data yv2(  8) /  -0.373737879689433D-01 /
  data yv2(  9) /  -0.531351030124350D-01 /
  data yv2( 10) /  -0.632864983555742D-01 /
  data yv2( 11) /  -0.548103623840369D-01 /
  data yv2( 12) /  -0.372889232343935D-01 /
  data yv2( 13) /  -0.188876517630154D-01 /
  data yv2( 14) /  -0.585239362533525D-02 /
  data yv2( 15) /   0.000000000000000D+00 /
  data yv2( 16) /   0.000000000000000D+00 /
  data yv2( 17) /   0.000000000000000D+00 /
  data bv2(  1) /  -0.702739315585347D+01 /
  data bv2(  2) /  -0.333140549270729D+01 /
  data bv2(  3) /  -0.117329394261502D+01 /
  data bv2(  4) /  -0.306003283486901D+00 /
  data bv2(  5) /  -0.366656699104026D-01 /
  data bv2(  6) /   0.588330899204400D-01 /
  data bv2(  7) /   0.384220572312032D-01 /
  data bv2(  8) /  -0.390223173707191D-01 /
  data bv2(  9) /  -0.663882722510521D-01 /
  data bv2( 10) /  -0.312918894386669D-02 /
  data bv2( 11) /   0.590118945294245D-01 /
  data bv2( 12) /   0.757939459148246D-01 /
  data bv2( 13) /   0.643822548468606D-01 /
  data bv2( 14) /   0.399750987463792D-01 /
  data bv2( 15) /   0.177103852679117D-05 /
  data bv2( 16) /  -0.590423369301474D-06 /
  data bv2( 17) /   0.590654950414731D-06 /
  data cv2(  1) /   0.877545959718548D+01 /
  data cv2(  2) /   0.585407125495837D+01 /
  data cv2(  3) /   0.268820820643116D+01 /
  data cv2(  4) /   0.744718689404422D+00 /
  data cv2(  5) /   0.321378734769888D+00 /
  data cv2(  6) /   0.566263292669091D-01 /
  data cv2(  7) /  -0.137417679148505D+00 /
  data cv2(  8) /  -0.169124163201523D+00 /
  data cv2(  9) /   0.608037039066423D-01 /
  data cv2( 10) /   0.189589640245655D+00 /
  data cv2( 11) /   0.563784150384640D-01 /
  data cv2( 12) /   0.100486298765028D-01 /
  data cv2( 13) /  -0.552186092621482D-01 /
  data cv2( 14) /  -0.413902746758285D-01 /
  data cv2( 15) /  -0.116832934994489D+00 /
  data cv2( 16) /   0.233610871054729D-01 /
  data cv2( 17) /   0.233885865725971D-01 /
  data dv2(  1) /  -0.385449887634130D+01 /
  data dv2(  2) /  -0.417706040200591D+01 /
  data dv2(  3) /  -0.256425277368288D+01 /
  data dv2(  4) /  -0.558557503589276D+00 /
  data dv2(  5) /  -0.349316054551627D+00 /
  data dv2(  6) /  -0.256022933201611D+00 /
  data dv2(  7) /  -0.418337423301704D-01 /
  data dv2(  8) /   0.303368330939646D+00 /
  data dv2(  9) /   0.169921006301015D+00 /
  data dv2( 10) /  -0.175759761362548D+00 /
  data dv2( 11) /  -0.611278214082881D-01 /
  data dv2( 12) /  -0.861140219824535D-01 /
  data dv2( 13) /   0.182451950513387D-01 /
  data dv2( 14) /  -0.995395392057973D-01 /
  data dv2( 15) /   0.184972909229936D+04 /
  data dv2( 16) /   0.362829766922787D+00 /
  data dv2( 17) /   0.362829766922787D+00 /

  if ( argmax <= arg ) then
    func =   0.000000000000000D+00
    dfunc  = 0.0D+00
    d2func = 0.0D+00
    return
  end if

  call seval ( nv2, arg, xv2, yv2, bv2, cv2, dv2, func, dfunc, d2func )

  return
end
subroutine rh ( arg, func, dfunc, d2func )

!*******************************************************************************
!
!! RH stores data for the Aluminum atomic density.
!
!        Aluminum  : atomic density   and its first two derivatives.
!        Generated automatically by PoCo, version 04-may-93           
!        Hamiltonian type #  2, run on 93/06/09 at 15.04.43
!        Uses subroutine seval from netlib@ornl.gov [to get it,
!        use 'send seval from sfmm'], trivially modified to
!        compute also dfunc and d2func and use double precision.
!
!  Modified:
!
!    28 October 2005
!
  implicit real ( kind = 8 ) (a-h,o-z)

  parameter (nrh= 17)
  parameter (argmax=   0.555805441821810D+01)
  real ( kind = 8 ) xrh(nrh),yrh(nrh),brh(nrh),crh(nrh),drh(nrh)
  save xrh,yrh,brh,crh,drh
  data xrh(  1) /   0.202111069753385D+01 /
  data xrh(  2) /   0.227374953472558D+01 /
  data xrh(  3) /   0.252638837191732D+01 /
  data xrh(  4) /   0.277902720910905D+01 /
  data xrh(  5) /   0.303166604630078D+01 /
  data xrh(  6) /   0.328430488349251D+01 /
  data xrh(  7) /   0.353694372068424D+01 /
  data xrh(  8) /   0.378958255787597D+01 /
  data xrh(  9) /   0.404222139506771D+01 /
  data xrh( 10) /   0.429486023225944D+01 /
  data xrh( 11) /   0.454749906945117D+01 /
  data xrh( 12) /   0.480013790664290D+01 /
  data xrh( 13) /   0.505277674383463D+01 /
  data xrh( 14) /   0.530541558102636D+01 /
  data xrh( 15) /   0.555805441821810D+01 /
  data xrh( 16) /   0.555807968210182D+01 /
  data xrh( 17) /   0.555810494598553D+01 /
  data yrh(  1) /   0.865674623712589D-01 /
  data yrh(  2) /   0.925214702944478D-01 /
  data yrh(  3) /   0.862003123832002D-01 /
  data yrh(  4) /   0.762736292751052D-01 /
  data yrh(  5) /   0.606481841271735D-01 /
  data yrh(  6) /   0.466030959588197D-01 /
  data yrh(  7) /   0.338740138848363D-01 /
  data yrh(  8) /   0.232572661705343D-01 /
  data yrh(  9) /   0.109046405489829D-01 /
  data yrh( 10) /   0.524910605677597D-02 /
  data yrh( 11) /   0.391702419142291D-02 /
  data yrh( 12) /   0.308277776293383D-02 /
  data yrh( 13) /   0.250214745349505D-02 /
  data yrh( 14) /   0.147220513798186D-02 /
  data yrh( 15) /   0.000000000000000D+00 /
  data yrh( 16) /   0.000000000000000D+00 /
  data yrh( 17) /   0.000000000000000D+00 /
  data brh(  1) /   0.608555214104682D-01 /
  data brh(  2) /  -0.800158928716306D-02 /
  data brh(  3) /  -0.332089451111092D-01 /
  data brh(  4) /  -0.521001991705069D-01 /
  data brh(  5) /  -0.618130637429111D-01 /
  data brh(  6) /  -0.529750064268036D-01 /
  data brh(  7) /  -0.442210477548108D-01 /
  data brh(  8) /  -0.473645664984640D-01 /
  data brh(  9) /  -0.390741582571631D-01 /
  data brh( 10) /  -0.101795580610560D-01 /
  data brh( 11) /  -0.318316981110289D-02 /
  data brh( 12) /  -0.281217210746153D-02 /
  data brh( 13) /  -0.236932031483360D-02 /
  data brh( 14) /  -0.683554708271547D-02 /
  data brh( 15) /  -0.638718204858808D-06 /
  data brh( 16) /   0.212925486831149D-06 /
  data brh( 17) /  -0.212983742465787D-06 /
  data crh(  1) /  -0.170233687052940D+00 /
  data crh(  2) /  -0.102317878901959D+00 /
  data crh(  3) /   0.254162872544396D-02 /
  data crh(  4) /  -0.773173610292656D-01 /
  data crh(  5) /   0.388717099948882D-01 /
  data crh(  6) /  -0.388873819867093D-02 /
  data crh(  7) /   0.385388290924526D-01 /
  data crh(  8) /  -0.509815666327127D-01 /
  data crh(  9) /   0.837968231208082D-01 /
  data crh( 10) /   0.305743500420042D-01 /
  data crh( 11) /  -0.288110886134041D-02 /
  data crh( 12) /   0.434959924771674D-02 /
  data crh( 13) /  -0.259669459714693D-02 /
  data crh( 14) /  -0.150816117849093D-01 /
  data crh( 15) /   0.421356801161513D-01 /
  data crh( 16) /  -0.842575249165724D-02 /
  data crh( 17) /  -0.843267014952237D-02 /
  data drh(  1) /   0.896085612514625D-01 /
  data drh(  2) /   0.138352319847830D+00 /
  data drh(  3) /  -0.105366473134009D+00 /
  data drh(  4) /   0.153300619856764D+00 /
  data drh(  5) /  -0.564184148788224D-01 /
  data drh(  6) /   0.559792096400504D-01 /
  data drh(  7) /  -0.118113795329664D+00 /
  data drh(  8) /   0.177827488509794D+00 /
  data drh(  9) /  -0.702220789044304D-01 /
  data drh( 10) /  -0.441413511810337D-01 /
  data drh( 11) /   0.954024354744484D-02 /
  data drh( 12) /  -0.916498550800407D-02 /
  data drh( 13) /  -0.164726813535368D-01 /
  data drh( 14) /   0.754928689733184D-01 /
  data drh( 15) /  -0.667110847110954D+03 /
  data drh( 16) /  -0.912720300911022D-01 /
  data drh( 17) /  -0.912720300911022D-01 /

  if ( argmax <= arg ) then
    func =   0.000000000000000D+00
    dfunc  = 0.0D+00
    d2func = 0.0D+00
    return
  end if

  call seval ( nrh, arg, xrh, yrh, brh, crh, drh, func, dfunc, d2func )

  return
end
subroutine uu ( arg, func, dfunc, d2func )

!*******************************************************************************
!
!! UU stores data for the Aluminum glue function.
!
!        Aluminum  : glue function    and its first two derivatives.
!        Generated automatically by PoCo, version 04-may-93           
!        Hamiltonian type #  2, run on 93/06/09 at 15.04.43
!        Uses subroutine seval from netlib@ornl.gov [to get it,
!        use 'send seval from sfmm'], trivially modified to
!        compute also dfunc and d2func and use double precision.
!
!  Modified:
!
!    28 October 2005
!
  implicit real ( kind = 8 ) (a-h,o-z)

  integer, parameter :: nuu = 13

  real ( kind = 8 ) :: argmin = 0.000000000000000D+00
  real ( kind = 8 ), dimension ( nuu ), save :: buu
  real ( kind = 8 ), dimension ( nuu ), save :: cuu
  real ( kind = 8 ), dimension ( nuu ), save :: duu
  real ( kind = 8 ), dimension ( nuu ), save :: xuu
  real ( kind = 8 ), dimension ( nuu ), save :: yuu

  data xuu(  1) /   0.000000000000000D+00 /
  data xuu(  2) /   0.100000000000000D+00 /
  data xuu(  3) /   0.200000000000000D+00 /
  data xuu(  4) /   0.300000000000000D+00 /
  data xuu(  5) /   0.400000000000000D+00 /
  data xuu(  6) /   0.500000000000000D+00 /
  data xuu(  7) /   0.600000000000000D+00 /
  data xuu(  8) /   0.700000000000000D+00 /
  data xuu(  9) /   0.800000000000000D+00 /
  data xuu( 10) /   0.900000000000000D+00 /
  data xuu( 11) /   0.100000000000000D+01 /
  data xuu( 12) /   0.110000000000000D+01 /
  data xuu( 13) /   0.120000000000000D+01 /
  data yuu(  1) /   0.000000000000000D+00 /
  data yuu(  2) /  -0.113953324143752D+01 /
  data yuu(  3) /  -0.145709859805864D+01 /
  data yuu(  4) /  -0.174913308002738D+01 /
  data yuu(  5) /  -0.202960322136630D+01 /
  data yuu(  6) /  -0.225202324967546D+01 /
  data yuu(  7) /  -0.242723053979436D+01 /
  data yuu(  8) /  -0.255171976467357D+01 /
  data yuu(  9) /  -0.260521638832322D+01 /
  data yuu( 10) /  -0.264397894381693D+01 /
  data yuu( 11) /  -0.265707884842034D+01 /
  data yuu( 12) /  -0.264564149400021D+01 /
  data yuu( 13) /  -0.260870604452106D+01 /
  data buu(  1) /  -0.183757286015853D+02 /
  data buu(  2) /  -0.574233124410516D+01 /
  data buu(  3) /  -0.236790436375322D+01 /
  data buu(  4) /  -0.307404645857774D+01 /
  data buu(  5) /  -0.251104850116555D+01 /
  data buu(  6) /  -0.196846462620234D+01 /
  data buu(  7) /  -0.154391254686695D+01 /
  data buu(  8) /  -0.846780636273251D+00 /
  data buu(  9) /  -0.408540363905760D+00 /
  data buu( 10) /  -0.286833282404628D+00 /
  data buu( 11) /  -0.309389414590161D-06 /
  data buu( 12) /   0.236958014464143D+00 /
  data buu( 13) /   0.503352368511243D+00 /
  data cuu(  1) /   0.830779120415016D+02 /
  data cuu(  2) /   0.432560615333001D+02 /
  data cuu(  3) /  -0.951179272978074D+01 /
  data cuu(  4) /   0.245037178153561D+01 /
  data cuu(  5) /   0.317960779258630D+01 /
  data cuu(  6) /   0.224623095704576D+01 /
  data cuu(  7) /   0.199928983630817D+01 /
  data cuu(  8) /   0.497202926962879D+01 /
  data cuu(  9) /  -0.589626545953876D+00 /
  data cuu( 10) /   0.180669736096520D+01 /
  data cuu( 11) /   0.106163236918694D+01 /
  data cuu( 12) /   0.130795086934864D+01 /
  data cuu( 13) /   0.135599267112235D+01 /
  data duu(  1) /  -0.132739501694005D+03 /
  data duu(  2) /  -0.175892847543603D+03 /
  data duu(  3) /   0.398738817043878D+02 /
  data duu(  4) /   0.243078670350231D+01 /
  data duu(  5) /  -0.311125611846847D+01 /
  data duu(  6) /  -0.823137069125319D+00 /
  data duu(  7) /   0.990913144440207D+01 /
  data duu(  8) /  -0.185388527186089D+02 /
  data duu(  9) /   0.798774635639692D+01 /
  data duu( 10) /  -0.248354997259420D+01 /
  data duu( 11) /   0.821061667205675D+00 /
  data duu( 12) /   0.160139339245701D+00 /
  data duu( 13) /   0.160139339245701D+00 /

  if ( arg <= argmin ) then
    func =   0.000000000000000D+00
    dfunc  = 0.0D+00
    d2func = 0.0D+00
    return
  end if

  call seval ( nuu, arg, xuu, yuu, buu, cuu, duu, func, dfunc, d2func )

  return
end
subroutine seval ( n, u, x, y, b, c, d, f, df, d2f )

!*******************************************************************************
!
!! SEVAL evaluates a cubic spline function.
!
!  Discussion:
!
!    This subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
!  Modified:
!
!    28 October 2005
!
  implicit none

  integer n

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) d2f
  real ( kind = 8 ) df
  real ( kind = 8 ) dx
  real ( kind = 8 ) f
  integer, save :: i = 1
  integer j
  integer k
  real ( kind = 8 ) u
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  if ( i .ge. n ) i = 1
  if ( u .lt. x(i) ) go to 10
  if ( u .le. x(i+1) ) go to 30
!
!  binary search
!
10 continue

  i = 1
  j = n+1

20 continue

  k = (i+j)/2
  if ( u .lt. x(k) ) j = k
  if ( u .ge. x(k) ) i = k
  if ( j .gt. i+1 ) go to 20
!
!  Evaluate the spline.
!
   30 continue

  dx = u - x(i)
  f = y(i) + dx * ( b(i) + dx * ( c(i) + dx * d(i) ) )
  df = b(i) + dx * ( 2.0D+00 * c(i) + 3.0D+00 * dx * d(i) )
  d2f = 2.0D+00 * c(i) + 6.0D+00 * dx * d(i)

  return
end
