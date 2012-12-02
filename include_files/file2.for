c
c=======================================================================
c
c     12/16/04 - Modified so that nq(1D) increases by 1 for
c                each increase in iqp?
c
c     02/16/04 - Allow for unequally spaced p-nodes.
c
c                Deleted functions f_rho & f_drho (not used).
c
c     12/15/03 - Fixed iorder computation in josh1d
c
c=======================================================================
c
c     Routine to be provided by Josh.
c
c     2D version.
c
c     Input:
c
c       itype - Element type (3 or 4)
c       iqp   - Quadrature indicator (1 - standard order, 2 - higher).
c       lp    - P-level (0 - (bi)linear, etc.)
c       nn    - Number of nodes on input. Standard order.
c
c     Output:
c
c       nq    - Number of quadrature nodes.
c       xq,yq - Quadrature node locations.
c       wq    - Quadrature weights.
c
c       phinq - Values of nodal basis functions at Q-nodes.
c       cqdx  - x-derivatives at quadrature nodes.
c       cqdy  - x-derivatives at quadrature nodes.
c
c     Simplest to call separate routines for information for
c     reference elements:
c
c     rquad   - Given itype, lp, nn, give nq,sqr,tqr,wqr.
c     dpdst   - Given itype, nq,..., return cqds, cqdt.
c
c     On entry, nq is the number of nodes
c
c     Assumes nq points can integrate a polynomial of order 2*nq-1.
c
c       n = nn-1 - Order of polynomials.
c
c       iqp = 1  - Suitable for integrating a product of two
c                  n'th order polynomials under a linear map.
c                  Order is 2*n.
c
c       iqp = 2  - Suitable for integrating a product of two n'th
c                  order polynomials under an isoparameteric map.
c                  Order is 5*n-1.
c
c       nq is taken to be order/2+1 (integer arithmetic)..
c
c     General notation:
c
c       rho - 1D standard nodal basis functions.
c       chi - 1D nodal basis functions defined at quadrature nodes.
c       phi - 2D nodal basis functions.
c
c=======================================================================
c
      subroutine josh1d(iqp,lp,nn,x,y,nq,xq,yq,wq,rhonq,cqdt,alen)
c
c=======================================================================
c
c     Find Gaussian Quadrature points/weights on an element edge.
c
c     This involves line integrals.
c
c     Let t(x(s),y(s)) be an arclength parameter. I need
c     int f(t) dt = int f(s)*dt/ds ds
c
c     dt/ds = sqrt(dx/ds^2+dy/ds^2).
c
c     This is because t(s) = int sqrt(dx/ds^2+dy/ds^2) ds,
c     so t'(s)
c
c=======================================================================
c
      implicit real*8(a-h,o-z)
c
c
c     02/27/04 - Added values:
c
c       nuc_max - Max number of user-coded functions.
c       ntl_max - Max number of table functions. (already defined)
c       npl_max - Max number of polynomial functions.
c       nsf_max - Max number of segmented functions.
c       nrg_max - Max number of registry values.
c
c
c     03/27/03 - Added nev_max (number of element vertices).
c                nne_max was used in itn, itb, etc. This was
c                too large, since at most 4 are used.
c
c     03/09/03 - Added some parameters for higher order.
c                Also arranged parameters into groups.
c
c     05/23/02 - Added table bounds.
c
c     These parameters (all of the form xxxx_max) are used for defining
c     the maximum numbers of equations, boundary conditions, zones, etc.
c     allowed in the code (particularly in readop and seta). This file
c     must be included (using the line "include 'params.max'") before
c     including any of the following files:
c
c       common.eqs   (storage for interpreted equations)
c       common.bcs   (storage for boundary conditions)
c       common.tre   (storage for trees)
c       common.ins   (storage for instruction list)
c       common.stc   (storage for stored constants)
c       common.all   List of names for base functions.
c       common.ho    Storage of higher-order info.
c
c     Order of finite element spaces.
c     -------------------------------
c
c       iqp_max  - Order indicator for finite elements.
c       lp_kmax  - Maximum "kick" order.
c       lp_lmax  - Maximum number of p-refinement levels.
c
c     Note: Using maximum lp_kmax & lp_lmax gives
c
c             nns_max = 1+(lp_kmax+1)*2^(lp_lmax-1)
c
c           Then, iqp_max = 1 gives
c
c             2*nqs_max-1 >= 2*(nns_max-1)
c
c     Any combination of k & lp with 1+(k+1)*2^(lp-1) <= nns_max is ok.
c
      parameter(iqp_max =4 )
      parameter(lp_kmax =2 )
      parameter(lp_lmax =5 )
c
c     Geometry
c     --------
c
c       nz_max  - Maximum number of zones allowed.
c       ns_max  - Maximum number of boundary numbers in domain.
c
      parameter(nz_max =10)
      parameter(ns_max =20)
c
c       nev_max - Maximum number of vertices per element.
c       nns_max - Maximum number of nodes/element side.
c       nqs_max - Maximum number of quadrature points/element side.
c     * nqb_max - Maximum number of quadrature points (edge).
c
      parameter(nev_max=4 )
      parameter(nns_max=10)
      parameter(nqs_max=10)
c     parameter(nqb_max=2 )  Obsolete. Same as nqs_max.
c
c       nne_max - Maximum number of nodes/element.
c       nne_max - Maximum number of nodes/element.
c       nqp_max - Maximum number of quadrature points (element).
c
        parameter(nne_max=nns_max*nns_max)
        parameter(nqp_max=nqs_max*nqs_max)
c
c     Equation reading, processing & storage.
c     ---------------------------------------
c
c       ncl_max - Maximum number of continuation lines allowed.
c       ntr_max - Maximum number of storage for trees.
c       nin_max - Maximum size of instruction list.
c
      parameter(ncl_max=5 )
      parameter(ntr_max=20000)
      parameter(nin_max=20000)
c
c       ncn_max - Maximum number of named constants.
c       nbf_max - Maximum number of base functions.
c       nfn_max - Maximum number of defined functions.
c       nun_max - Maximum number of unknowns allowed
c       nmt_max - Maximum number of matrices.
c       nsc_max - Maximum number of stored constants.
c
c       nis_max - Maximum number of initialization specs.
c
      parameter(ncn_max=50)
      parameter(nbf_max=50)
      parameter(nfn_max=100)
      parameter(nun_max=10)
      parameter(nmt_max=50)
      parameter(nsc_max=5000)
c
      parameter(nis_max=50)
c
c       ne_max  - Average number of entries per matrix.
c       nme_max - Maximum total number of matrix entries.
c
      parameter(ne_max=4)
        parameter(nme_max=nmt_max*ne_max)
c
c       nam_max - Maximum number of names allowed.
c
        parameter(nam_max=ncn_max+nbf_max+nun_max+nfn_max+nmt_max)
c
c       neq_max - Maximum number of equations allowed
c       nbc_max - Maximum number of boundary conditions.
c       nt_max  - Maximum average number of terms per equation.
c       ntx_max - Maximum number of terms per equation.
c       ntb_max - Maximum number of terms allowed in bc's.
c
      parameter(neq_max=50)
      parameter(nbc_max=50)
      parameter(nt_max =50)
      parameter(ntx_max=500)
        parameter(ntq_max=neq_max*nt_max)
        parameter(ntb_max=nt_max*nbc_max)
c
c       ntl_max - Maximum number of tables.
c       ntd_max - Maximum space for table data.
c
      parameter(ntl_max=10)
      parameter(ntd_max=100*ntl_max)
c
c     02/27/04
c
c       nuc_max - Max number of user-coded functions.
c       npl_max - Max number of polynomial functions.
c       nsf_max - Max number of segmented functions.
c       nrg_max - Max number of registry values.
c
c       ntl_max - Maximum number of tables.
c       ntd_max - Maximum space for table data.
c
      parameter(nuc_max=50)
      parameter(npl_max=50)
      parameter(npd_max=10*npl_max)
      parameter(nsf_max=50)
      parameter(nsd_max=10*nsf_max)
      parameter(nrg_max=100)
c
      dimension x (nn)
      dimension y (nn)
c
      dimension xq(nqs_max)
      dimension yq(nqs_max)
      dimension wq(nqs_max)
c
      dimension sq(nqs_max)
c
      dimension rhonq(nns_max,nqs_max)
      dimension cqdt (nqs_max,nqs_max)
c
c     Internal storage.
c
      dimension drho(nns_max,nqs_max)
      dimension dchi(nqs_max,nqs_max)
c
      save
c
c=======================================================================
c
c     Determine order of polynomials to be integrated exactly.
c
      iorder = 2*(nn+iqp-2)
c     if(iqp.eq.1) then
c       iorder = 2*(nn-1)
c     elseif(iqp.eq.2) then
c       iorder = 4*(nn-1)-1
c     elseif(iqp.eq.3) then
c       iorder = 8*(nn-1)-1
c     endif
c
c     Get 1D quadrature nodes.
c
      call quad1d(iorder,nq,sq,wq)
      if(nq.gt.nqs_max) stop 'nq > nqs_max'
c
c     Get derivatives of chi's at quadrature nodes.
c
      call getdchi(nq,sq,dchi,nqs_max)
c
c     Get 1D basis functions.
c
      call getrho (nn,nq,sq,rhonq,nns_max)
c
c     Get derivatives of 1D basis functions at q-nodes.
c
      call getdrho(nn,nq,sq,drho,nns_max)
c
c     Compute dx/ds and dy/ds at quad nodes.
c
      alen = 0.0d0
      do iq=1,nq
        xq(iq) = 0.0d0
        yq(iq) = 0.0d0
        do in=1,nn
          xq(iq) = xq(iq)+x(in)*rhonq(in,iq)
          yq(iq) = yq(iq)+y(in)*rhonq(in,iq)
        enddo
        dxds = 0.0d0
        dyds = 0.0d0
        do in=1,nn
          dxds = dxds+x(in)*drho(in,iq)
          dyds = dyds+y(in)*drho(in,iq)
        enddo
        dtds = dsqrt(dxds*dxds+dyds*dyds)
        wq(iq) = wq(iq)*dtds
        alen = alen+wq(iq)
c
        dsdt = 1.0d0/dtds
c
        do jq=1,nq
          cqdt(iq,jq) = dchi(jq,iq)*dsdt
        enddo
      enddo
c
c     area = 0.0d0
c     do iq=1,nq
c       area = area+wq(iq)
c     enddo
c
      return
      end
c
c=======================================================================
c
      subroutine josh2d(itype,iqp,lp,nn,x,y,nq,xq,yq,wq,phinq,
     *                  cqdx,cqdy,area)
c
c=======================================================================
c
c     Find Gaussian Quadrature points/weights on the element.
c
c=======================================================================
c
      implicit real*8(a-h,o-z)
c
c
c     02/27/04 - Added values:
c
c       nuc_max - Max number of user-coded functions.
c       ntl_max - Max number of table functions. (already defined)
c       npl_max - Max number of polynomial functions.
c       nsf_max - Max number of segmented functions.
c       nrg_max - Max number of registry values.
c
c
c     03/27/03 - Added nev_max (number of element vertices).
c                nne_max was used in itn, itb, etc. This was
c                too large, since at most 4 are used.
c
c     03/09/03 - Added some parameters for higher order.
c                Also arranged parameters into groups.
c
c     05/23/02 - Added table bounds.
c
c     These parameters (all of the form xxxx_max) are used for defining
c     the maximum numbers of equations, boundary conditions, zones, etc.
c     allowed in the code (particularly in readop and seta). This file
c     must be included (using the line "include 'params.max'") before
c     including any of the following files:
c
c       common.eqs   (storage for interpreted equations)
c       common.bcs   (storage for boundary conditions)
c       common.tre   (storage for trees)
c       common.ins   (storage for instruction list)
c       common.stc   (storage for stored constants)
c       common.all   List of names for base functions.
c       common.ho    Storage of higher-order info.
c
c     Order of finite element spaces.
c     -------------------------------
c
c       iqp_max  - Order indicator for finite elements.
c       lp_kmax  - Maximum "kick" order.
c       lp_lmax  - Maximum number of p-refinement levels.
c
c     Note: Using maximum lp_kmax & lp_lmax gives
c
c             nns_max = 1+(lp_kmax+1)*2^(lp_lmax-1)
c
c           Then, iqp_max = 1 gives
c
c             2*nqs_max-1 >= 2*(nns_max-1)
c
c     Any combination of k & lp with 1+(k+1)*2^(lp-1) <= nns_max is ok.
c
      parameter(iqp_max =4 )
      parameter(lp_kmax =2 )
      parameter(lp_lmax =5 )
c
c     Geometry
c     --------
c
c       nz_max  - Maximum number of zones allowed.
c       ns_max  - Maximum number of boundary numbers in domain.
c
      parameter(nz_max =10)
      parameter(ns_max =20)
c
c       nev_max - Maximum number of vertices per element.
c       nns_max - Maximum number of nodes/element side.
c       nqs_max - Maximum number of quadrature points/element side.
c     * nqb_max - Maximum number of quadrature points (edge).
c
      parameter(nev_max=4 )
      parameter(nns_max=10)
      parameter(nqs_max=10)
c     parameter(nqb_max=2 )  Obsolete. Same as nqs_max.
c
c       nne_max - Maximum number of nodes/element.
c       nne_max - Maximum number of nodes/element.
c       nqp_max - Maximum number of quadrature points (element).
c
        parameter(nne_max=nns_max*nns_max)
        parameter(nqp_max=nqs_max*nqs_max)
c
c     Equation reading, processing & storage.
c     ---------------------------------------
c
c       ncl_max - Maximum number of continuation lines allowed.
c       ntr_max - Maximum number of storage for trees.
c       nin_max - Maximum size of instruction list.
c
      parameter(ncl_max=5 )
      parameter(ntr_max=20000)
      parameter(nin_max=20000)
c
c       ncn_max - Maximum number of named constants.
c       nbf_max - Maximum number of base functions.
c       nfn_max - Maximum number of defined functions.
c       nun_max - Maximum number of unknowns allowed
c       nmt_max - Maximum number of matrices.
c       nsc_max - Maximum number of stored constants.
c
c       nis_max - Maximum number of initialization specs.
c
      parameter(ncn_max=50)
      parameter(nbf_max=50)
      parameter(nfn_max=100)
      parameter(nun_max=10)
      parameter(nmt_max=50)
      parameter(nsc_max=5000)
c
      parameter(nis_max=50)
c
c       ne_max  - Average number of entries per matrix.
c       nme_max - Maximum total number of matrix entries.
c
      parameter(ne_max=4)
        parameter(nme_max=nmt_max*ne_max)
c
c       nam_max - Maximum number of names allowed.
c
        parameter(nam_max=ncn_max+nbf_max+nun_max+nfn_max+nmt_max)
c
c       neq_max - Maximum number of equations allowed
c       nbc_max - Maximum number of boundary conditions.
c       nt_max  - Maximum average number of terms per equation.
c       ntx_max - Maximum number of terms per equation.
c       ntb_max - Maximum number of terms allowed in bc's.
c
      parameter(neq_max=50)
      parameter(nbc_max=50)
      parameter(nt_max =50)
      parameter(ntx_max=500)
        parameter(ntq_max=neq_max*nt_max)
        parameter(ntb_max=nt_max*nbc_max)
c
c       ntl_max - Maximum number of tables.
c       ntd_max - Maximum space for table data.
c
      parameter(ntl_max=10)
      parameter(ntd_max=100*ntl_max)
c
c     02/27/04
c
c       nuc_max - Max number of user-coded functions.
c       npl_max - Max number of polynomial functions.
c       nsf_max - Max number of segmented functions.
c       nrg_max - Max number of registry values.
c
c       ntl_max - Maximum number of tables.
c       ntd_max - Maximum space for table data.
c
      parameter(nuc_max=50)
      parameter(npl_max=50)
      parameter(npd_max=10*npl_max)
      parameter(nsf_max=50)
      parameter(nsd_max=10*nsf_max)
      parameter(nrg_max=100)
c
      dimension x (nn)
      dimension y (nn)
c
      dimension xq(*)
      dimension yq(*)
      dimension wq(*)
c
      dimension phinq(nne_max,nqp_max)
      dimension cqdx (nqp_max,nqp_max)
      dimension cqdy (nqp_max,nqp_max)
c
c     1D data.
c
      dimension wq1(nqs_max)
c
c     Internal storage.
c
      dimension rho (nns_max,nqs_max)
      dimension drho(nns_max,nqs_max)
      dimension dchi(nqs_max,nqs_max)
c
      dimension sq(1000)
      dimension tq(1000)
c
      data iicall /0/
c
      save
c
c=======================================================================
c
      if(iicall.eq.0) then
        print *,'Josh2d - iqp = ',iqp
        iqp0 = iqp
        iicall = 1
      endif
      if(iqp.ne.iqp0) print *,'New iqp = ',iqp
c
c     Triangles.
c     ----------
c
      if(itype.eq.3) then
        ax = x(1)
        bx = x(3)-x(1)
        cx = x(2)-x(1)
        ay = y(1)
        by = y(3)-y(1)
        cy = y(2)-y(1)
        det = x(1)*(y(3)-y(2))+x(2)*(y(1)-y(3))+x(3)*(y(2)-y(1))
        if(iqp.eq.0) then
          nq = 3
          xc = 0.0d0
          yc = 0.0d0
          do i=1,3
            xc = xc+x(i)
            yc = yc+y(i)
          enddo
          xc = xc/3.0d0
          yc = yc/3.0d0
          do i=1,3
            xq(i) = (0.999d0*x(i)+0.001d0*xc)
            yq(i) = (0.999d0*y(i)+0.001d0*yc)
            do j=1,3
              phinq(j,i) = 0.0d0
            enddo
            phinq(i,i) = 1.0d0
          enddo
          wq(1) = dabs(det)/6.0d0
          wq(2) = dabs(det)/6.0d0
          wq(3) = dabs(det)/6.0d0
        else
          call quadt(iqp,lp,nn,nq,sq,tq,wq)
          nq = 4
          do iiq=1,nq
            xq(iiq) = ax+bx*sq(iiq)+cx*tq(iiq)
            yq(iiq) = ay+by*sq(iiq)+cy*tq(iiq)
            wq(iiq) = wq(iiq)*dabs(det)
c
            phinq(1,iiq) = 1.0d0-sq(iiq)-tq(iiq)
            phinq(2,iiq) = tq(iiq)
            phinq(3,iiq) = sq(iiq)
            phinq(4,iiq) = 0.0d0
          enddo
        endif
c
c     Quadrilaterals.
c     ---------------
c
      elseif(itype.eq.4) then
c
c       Get data for 1D
c
c       Determine order of polynomials to be integrated exactly.
c
        n1 = sqrt(float(nn))
        iorder = 2*(n1+iqp-2)
c       if(iqp.eq.1) then
c         iorder = 2*(n1-1)
c       elseif(iqp.eq.2) then
c         iorder = 4*(n1-1)-1
c       elseif(iqp.eq.3) then
c         iorder = 8*(n1-1)-1
c       endif
c
c       Get 1D quadrature nodes.
c
        call quad1d(iorder,nq1,sq,wq1)
c       print *,'  josh2d - iorder,n1,nq1=',iorder,n1,nq1
        if(nq1.gt.nqs_max) stop 'nq > nqs_max'
c
c       Get derivatives of chi's at quadrature nodes (1D).
c
        call getdchi(nq1,sq,dchi,nqs_max)
c
c       Get 1D basis functions.
c
        call getrho (n1,nq1,sq,rho,nns_max)
c
c       Get derivatives of 1D basis functions at q-nodes.
c
        call getdrho(n1,nq1,sq,drho,nns_max)
c
c       Get 2D basis functions.
c
        nn = n1*n1
        nq = nq1*nq1
        call getphinq(rho,nns_max,n1,nq1,phinq,nne_max,nn,nq)
c
c       Compute dx/ds, dx/dt, dy/ds and dy/dt at quad nodes.
c
        iq1 = 0
        do j1=1,nq1
          do i1=1,nq1
            iq1 = iq1+1
            dxds = 0.0d0
            dxdt = 0.0d0
            dyds = 0.0d0
            dydt = 0.0d0
            xq(iq1) = 0.0d0
            yq(iq1) = 0.0d0
            in = 0
            do j2=1,n1
              do i2=1,n1
                in = in+1
                dphids = drho(i2,i1)*rho(j2,j1)
                dphidt = rho(i2,i1)*drho(j2,j1)
                dxds = dxds+x(in)*dphids
                dxdt = dxdt+x(in)*dphidt
                dyds = dyds+y(in)*dphids
                dydt = dydt+y(in)*dphidt
                xq(iq1) = xq(iq1)+x(in)*phinq(in,iq1)
                yq(iq1) = yq(iq1)+y(in)*phinq(in,iq1)
              enddo
            enddo
            det = dxds*dydt-dxdt*dyds
c           write(9,'(''  J det = '',i3,2x,1p,d12.5)') iq1,det
            wq(iq1) = wq1(i1)*wq1(j1)*det
c
            dsdx =  dydt/det
            dsdy = -dxdt/det
            dtdx = -dyds/det
            dtdy =  dxds/det
c
c           Get coefficients for derivatives.
c
            do iq2=1,nq
              cqdx(iq1,iq2) = 0.0d0
              cqdy(iq1,iq2) = 0.0d0
            enddo
c
            cqdx(iq1,iq1) = dchi(i1,i1)*dsdx+dchi(j1,j1)*dtdx
            cqdy(iq1,iq1) = dchi(i1,i1)*dsdy+dchi(j1,j1)*dtdy
c
c           Location (i2,j1). s-line.
c
            j2 = j1
            do i2=1,nq1
              if(i2.ne.i1) then
                iq2 = (j2-1)*nq1+i2
                cqdx(iq1,iq2) = dchi(i2,i1)*dsdx
                cqdy(iq1,iq2) = dchi(i2,i1)*dsdy
              endif
            enddo
c
c           Location (i1,j2). t-line.
c
            i2 = i1
            do j2=1,nq1
              if(j2.ne.j1) then
                iq2 = (j2-1)*nq1+i2
                cqdx(iq1,iq2) = dchi(j2,j1)*dtdx
                cqdy(iq1,iq2) = dchi(j2,j1)*dtdy
              endif
            enddo
          enddo
        enddo
      endif
c
      area = 0.0d0
      do iq=1,nq
        area = area+wq(iq)
      enddo
c
      return
      end
c
c=======================================================================
c
      subroutine quadt(iqp,lp,nn,nq,sq,tq,wq)
c
c=======================================================================
c
c     Return 2D quadrature nodes & weights on reference triangle.
c
c     Currently only 4 & 5-node quadrature.
c
c=======================================================================
c
      implicit real*8(a-h,o-z)
c
      dimension sq(*)
      dimension tq(*)
      dimension wq(*)
c
c     Triangle - 4 node quadrature
c
      dimension s34(4)
      dimension t34(4)
      dimension w34(4)
c
      data w34 /1.718630568649268D-01,
     *          1.089432529485410D-01,
     *          1.089432529485410D-01,
     *          1.102504372379912D-01/
      data s34 /1.728105738621062D-01,
     *          9.078239356710568D-02,
     *          7.092176064328946D-01,
     *          4.518105512614261D-01/
      data t34 /1.728105738621062D-01,
     *          7.092176064328946D-01,
     *          9.078239356710568D-02,
     *          4.518105512614261D-01/
c
c     Triangle - 5 node quadrature
c
      dimension s35(5)
      dimension t35(5)
      dimension w35(5)
c
      data w35 /8.897730939166018D-02,
     *          1.017073718944203D-01,
     *          1.017073718944203D-01,
     *          8.897730939166018D-02,
     *          1.186306374278391D-01/
      data s35 /1.103471248136220D-01,
     *          9.435493353697824D-02,
     *          7.226288329447310D-01,
     *          2.410452222695360D-01,
     *          4.409269851976064D-01/
      data t35 /2.410452222695360D-01,
     *          7.226288329447310D-01,
     *          9.435493353697824D-02,
     *          1.103471248136220D-01,
     *          4.409269851976064D-01/
c
      save
c
c=======================================================================
c
c     Set weights, points for nodal quadrature.
c
c     (Set just inside elements for discontinuous coefficients.)
c
c     12/04/01 - Fixed to set quadrature nodes with same order
c                as with iqp > 0. This requires a change in wqn,
c                xq, and yq. Initially, the Q point ordering was
c                the same as the node ordering. Differentiation
c                coefficients were defined assuming the "usual"
c                Q point ordering, however.
c
c                      2----3             +----+
c                      |    |             |3  4|
c                Nodes |    |    Q Points |1  2|
c                      1----4             +----+
c
c                See notes on reversal of r & s at top of file.
c
      if(iqp.eq.0) then
        nq = nn
        stop 'iqp = 0 not in effect in quad2d'
      endif
c
c     Triangles.
c     ----------
c
c       f = a + b*r + c*s + d*t
c
      if(iqp.eq.1) then
        nq = 4
        do iiq=1,nq
          sq(iiq) = s34(iiq)
          tq(iiq) = t34(iiq)
          wq(iiq) = w34(iiq)
        enddo
c
      elseif(iqp.eq.2) then
        nq = 5
        do iiq=1,nq
          sq(iiq) = s35(iiq)
          tq(iiq) = t35(iiq)
          wq(iiq) = w35(iiq)
        enddo
      endif
c
      return
      end
c
c=======================================================================
c
      subroutine quad1d(iorder,nq,sq,wq)
c
c=======================================================================
c
c     Determine 1D quadrature nodes & weights.
c
c     Reference element is [-1,1].
c
c     iorder is the polynomial order to be integrated exactly.
c
c=======================================================================
c
      implicit real*8(a-h,o-z)
c
c
c     02/27/04 - Added values:
c
c       nuc_max - Max number of user-coded functions.
c       ntl_max - Max number of table functions. (already defined)
c       npl_max - Max number of polynomial functions.
c       nsf_max - Max number of segmented functions.
c       nrg_max - Max number of registry values.
c
c
c     03/27/03 - Added nev_max (number of element vertices).
c                nne_max was used in itn, itb, etc. This was
c                too large, since at most 4 are used.
c
c     03/09/03 - Added some parameters for higher order.
c                Also arranged parameters into groups.
c
c     05/23/02 - Added table bounds.
c
c     These parameters (all of the form xxxx_max) are used for defining
c     the maximum numbers of equations, boundary conditions, zones, etc.
c     allowed in the code (particularly in readop and seta). This file
c     must be included (using the line "include 'params.max'") before
c     including any of the following files:
c
c       common.eqs   (storage for interpreted equations)
c       common.bcs   (storage for boundary conditions)
c       common.tre   (storage for trees)
c       common.ins   (storage for instruction list)
c       common.stc   (storage for stored constants)
c       common.all   List of names for base functions.
c       common.ho    Storage of higher-order info.
c
c     Order of finite element spaces.
c     -------------------------------
c
c       iqp_max  - Order indicator for finite elements.
c       lp_kmax  - Maximum "kick" order.
c       lp_lmax  - Maximum number of p-refinement levels.
c
c     Note: Using maximum lp_kmax & lp_lmax gives
c
c             nns_max = 1+(lp_kmax+1)*2^(lp_lmax-1)
c
c           Then, iqp_max = 1 gives
c
c             2*nqs_max-1 >= 2*(nns_max-1)
c
c     Any combination of k & lp with 1+(k+1)*2^(lp-1) <= nns_max is ok.
c
      parameter(iqp_max =4 )
      parameter(lp_kmax =2 )
      parameter(lp_lmax =5 )
c
c     Geometry
c     --------
c
c       nz_max  - Maximum number of zones allowed.
c       ns_max  - Maximum number of boundary numbers in domain.
c
      parameter(nz_max =10)
      parameter(ns_max =20)
c
c       nev_max - Maximum number of vertices per element.
c       nns_max - Maximum number of nodes/element side.
c       nqs_max - Maximum number of quadrature points/element side.
c     * nqb_max - Maximum number of quadrature points (edge).
c
      parameter(nev_max=4 )
      parameter(nns_max=10)
      parameter(nqs_max=10)
c     parameter(nqb_max=2 )  Obsolete. Same as nqs_max.
c
c       nne_max - Maximum number of nodes/element.
c       nne_max - Maximum number of nodes/element.
c       nqp_max - Maximum number of quadrature points (element).
c
        parameter(nne_max=nns_max*nns_max)
        parameter(nqp_max=nqs_max*nqs_max)
c
c     Equation reading, processing & storage.
c     ---------------------------------------
c
c       ncl_max - Maximum number of continuation lines allowed.
c       ntr_max - Maximum number of storage for trees.
c       nin_max - Maximum size of instruction list.
c
      parameter(ncl_max=5 )
      parameter(ntr_max=20000)
      parameter(nin_max=20000)
c
c       ncn_max - Maximum number of named constants.
c       nbf_max - Maximum number of base functions.
c       nfn_max - Maximum number of defined functions.
c       nun_max - Maximum number of unknowns allowed
c       nmt_max - Maximum number of matrices.
c       nsc_max - Maximum number of stored constants.
c
c       nis_max - Maximum number of initialization specs.
c
      parameter(ncn_max=50)
      parameter(nbf_max=50)
      parameter(nfn_max=100)
      parameter(nun_max=10)
      parameter(nmt_max=50)
      parameter(nsc_max=5000)
c
      parameter(nis_max=50)
c
c       ne_max  - Average number of entries per matrix.
c       nme_max - Maximum total number of matrix entries.
c
      parameter(ne_max=4)
        parameter(nme_max=nmt_max*ne_max)
c
c       nam_max - Maximum number of names allowed.
c
        parameter(nam_max=ncn_max+nbf_max+nun_max+nfn_max+nmt_max)
c
c       neq_max - Maximum number of equations allowed
c       nbc_max - Maximum number of boundary conditions.
c       nt_max  - Maximum average number of terms per equation.
c       ntx_max - Maximum number of terms per equation.
c       ntb_max - Maximum number of terms allowed in bc's.
c
      parameter(neq_max=50)
      parameter(nbc_max=50)
      parameter(nt_max =50)
      parameter(ntx_max=500)
        parameter(ntq_max=neq_max*nt_max)
        parameter(ntb_max=nt_max*nbc_max)
c
c       ntl_max - Maximum number of tables.
c       ntd_max - Maximum space for table data.
c
      parameter(ntl_max=10)
      parameter(ntd_max=100*ntl_max)
c
c     02/27/04
c
c       nuc_max - Max number of user-coded functions.
c       npl_max - Max number of polynomial functions.
c       nsf_max - Max number of segmented functions.
c       nrg_max - Max number of registry values.
c
c       ntl_max - Maximum number of tables.
c       ntd_max - Maximum space for table data.
c
      parameter(nuc_max=50)
      parameter(npl_max=50)
      parameter(npd_max=10*npl_max)
      parameter(nsf_max=50)
      parameter(nsd_max=10*nsf_max)
      parameter(nrg_max=100)
c
      dimension sq(*)
      dimension wq(*)
c
      dimension q2(2)
      data q2  /-5.773502691896257E-01, 5.773502691896257E-01/
c
      dimension w2(2)
      data w2  / 1.000000000000000E+00, 1.000000000000000E+00/
c
      dimension q3(3)
      data q3  /-7.745966692414834E-01, 0.000000000000000E+00,
     *           7.745966692414834E-01/
c
      dimension w3(3)
      data w3  / 5.555555555555556E-01, 8.888888888888890E-01,
     *           5.555555555555556E-01/
c
      dimension q4(4)
      data q4  /-8.611363115940527E-01,-3.399810435848563E-01,
     *           3.399810435848563E-01, 8.611363115940527E-01/
c
      dimension w4(4)
      data w4  / 3.478548451374538E-01, 6.521451548625462E-01,
     *           6.521451548625462E-01, 3.478548451374538E-01/
c
      dimension q5(5)
      data q5  /-9.061798459386643E-01,-5.384693101056841E-01,
     *           0.000000000000000E+00, 5.384693101056841E-01,
     *           9.061798459386643E-01/
c
      dimension w5(5)
      data w5  / 2.369268850561884E-01, 4.786286704993661E-01,
     *           5.688888888888909E-01, 4.786286704993661E-01,
     *           2.369268850561884E-01/
c
      dimension q6(6)
      data q6  /-9.324695142031521E-01,-6.612093864662638E-01,
     *          -2.386191860831955E-01, 2.386191860831955E-01,
     *           6.612093864662638E-01, 9.324695142031521E-01/
c
      dimension w6(6)
      data w6  / 1.713244923791704E-01, 3.607615730481403E-01,
     *           4.679139345726893E-01, 4.679139345726893E-01,
     *           3.607615730481403E-01, 1.713244923791704E-01/
c
      dimension q7(7)
      data q7  /-9.491079123427586E-01,-7.415311855993956E-01,
     *          -4.058451513774004E-01, 0.000000000000000E+00,
     *           4.058451513774004E-01, 7.415311855993956E-01,
     *           9.491079123427586E-01/
c
      dimension w7(7)
      data w7  / 1.294849661688693E-01, 2.797053914892749E-01,
     *           3.818300505051175E-01, 4.179591836734764E-01,
     *           3.818300505051175E-01, 2.797053914892749E-01,
     *           1.294849661688693E-01/
c
      dimension q8(8)
      data q8  /-9.602898564975422E-01,-7.966664774136561E-01,
     *          -5.255324099163844E-01,-1.834346424956806E-01,
     *           1.834346424956806E-01, 5.255324099163844E-01,
     *           7.966664774136561E-01, 9.602898564975422E-01/
c
      dimension w8(8)
      data w8  / 1.012285362903613E-01, 2.223810344533440E-01,
     *           3.137066458778776E-01, 3.626837833784173E-01,
     *           3.626837833784173E-01, 3.137066458778776E-01,
     *           2.223810344533440E-01, 1.012285362903613E-01/
c
      dimension q9(9)
      data q9  /-9.681602395076134E-01,-8.360311073265632E-01,
     *          -6.133714327004148E-01,-3.242534234036133E-01,
     *           0.000000000000000E+00, 3.242534234036133E-01,
     *           6.133714327004148E-01, 8.360311073265632E-01,
     *           9.681602395076134E-01/
c
      dimension w9(9)
      data w9  / 8.127438836160783E-02, 1.806481606949452E-01,
     *           2.606106964030327E-01, 3.123470770399117E-01,
     *           3.302393550010050E-01, 3.123470770399117E-01,
     *           2.606106964030327E-01, 1.806481606949452E-01,
     *           8.127438836160783E-02/
c
      dimension q10(10)
      data q10 /-9.739065285171347E-01,-8.650633666888045E-01,
     *          -6.794095682986596E-01,-4.333953941288043E-01,
     *          -1.488743389814015E-01, 1.488743389814015E-01,
     *           4.333953941288043E-01, 6.794095682986596E-01,
     *           8.650633666888045E-01, 9.739065285171347E-01/
c
      dimension w10(10)
      data w10 / 6.667134430878152E-02, 1.494513491507625E-01,
     *           2.190863625161436E-01, 2.692667193099595E-01,
     *           2.955242247143529E-01, 2.955242247143529E-01,
     *           2.692667193099595E-01, 2.190863625161436E-01,
     *           1.494513491507625E-01, 6.667134430878152E-02/
c
      dimension q11(11)
      data q11 /-9.782286581466345E-01,-8.870625997708689E-01,
     *          -7.301520055796078E-01,-5.190961292137052E-01,
     *          -2.695431559572268E-01, 0.000000000000000E+00,
     *           2.695431559572268E-01, 5.190961292137052E-01,
     *           7.301520055796078E-01, 8.870625997708689E-01,
     *           9.782286581466345E-01/

      dimension w11(11)
      data w11 / 5.566856711472156E-02, 1.255803694621458E-01,
     *           1.862902109252768E-01, 2.331937645921465E-01,
     *           2.628045445140306E-01, 2.729250867833574E-01,
     *           2.628045445140306E-01, 2.331937645921465E-01,
     *           1.862902109252768E-01, 1.255803694621458E-01,
     *           5.566856711472156E-02/
c
      dimension q12(12)
      data q12 /-9.815606342448744E-01,-9.041172563614603E-01,
     *          -7.699026741754944E-01,-5.873179542609674E-01,
     *          -3.678314989748968E-01,-1.252334085018861E-01,
     *           1.252334085018861E-01, 3.678314989748968E-01,
     *           5.873179542609674E-01, 7.699026741754944E-01,
     *           9.041172563614603E-01, 9.815606342448744E-01/
c
      dimension w12(12)
      data w12 / 4.717533639116828E-02, 1.069393260045217E-01,
     *           1.600783285527755E-01, 2.031674267261821E-01,
     *           2.334925365300782E-01, 2.491470457952742E-01,
     *           2.491470457952742E-01, 2.334925365300782E-01,
     *           2.031674267261821E-01, 1.600783285527755E-01,
     *           1.069393260045217E-01, 4.717533639116828E-02/
c
      dimension q13(13)
      data q13 /-9.841830547302339E-01,-9.175983992813748E-01,
     *          -8.015780908620191E-01,-6.423493396334936E-01,
     *          -4.484927512426413E-01,-2.304583160899322E-01,
     *           0.000000000000000E+00, 2.304583160899322E-01,
     *           4.484927512426413E-01, 6.423493396334936E-01,
     *           8.015780908620191E-01, 9.175983992813748E-01,
     *           9.841830547302339E-01/
c
      dimension w13(13)
      data w13 / 4.048400473575120E-02, 9.212149977586202E-02,
     *           1.388735101458644E-01, 1.781459807150328E-01,
     *           2.078160475642222E-01, 2.262831803744856E-01,
     *           2.325515533775633E-01, 2.262831803744856E-01,
     *           2.078160475642222E-01, 1.781459807150328E-01,
     *           1.388735101458644E-01, 9.212149977586202E-02,
     *           4.048400473575120E-02/
c
      dimension q14(14)
      data q14 /-9.862838086912444E-01,-9.284348836333490E-01,
     *          -8.272013149923210E-01,-6.872929046647064E-01,
     *          -5.152486361390562E-01,-3.191123686989603E-01,
     *          -1.080549486040908E-01, 1.080549486040908E-01,
     *           3.191123686989603E-01, 5.152486361390562E-01,
     *           6.872929046647064E-01, 8.272013149923210E-01,
     *           9.284348836333490E-01, 9.862838086912444E-01/
c
      dimension w14(14)
      data w14 / 3.511946034614696E-02, 8.015808719520255E-02,
     *           1.215185707471110E-01, 1.572031672353473E-01,
     *           1.855383975339648E-01, 2.051984636695097E-01,
     *           2.152638532727177E-01, 2.152638532727177E-01,
     *           2.051984636695097E-01, 1.855383975339648E-01,
     *           1.572031672353473E-01, 1.215185707471110E-01,
     *           8.015808719520255E-02, 3.511946034614696E-02/
c
      dimension q15(15)
      data q15 /-9.879925180289779E-01,-9.372733924338538E-01,
     *          -8.482065834429273E-01,-7.244177313068012E-01,
     *          -5.709721723633163E-01,-3.941513466458632E-01,
     *          -2.011940936203377E-01, 0.000000000000000E+00,
     *           2.011940936203377E-01, 3.941513466458632E-01,
     *           5.709721723633163E-01, 7.244177313068012E-01,
     *           8.482065834429273E-01, 9.372733924338538E-01,
     *           9.879925180289779E-01/
c
      dimension w15(15)
      data w15 / 3.075324197565902E-02, 7.036604746656214E-02,
     *           1.071592205008836E-01, 1.395706780695285E-01,
     *           1.662692060392694E-01, 1.861610001218824E-01,
     *           1.984314850909571E-01, 2.025782414705157E-01,
     *           1.984314850909571E-01, 1.861610001218824E-01,
     *           1.662692060392694E-01, 1.395706780695285E-01,
     *           1.071592205008836E-01, 7.036604746656214E-02,
     *           3.075324197565902E-02/
c
      dimension q16(16)
      data q16 /-9.894009348053501E-01,-9.445750221571226E-01,
     *          -8.656312004236659E-01,-7.554044054509828E-01,
     *          -6.178762411475158E-01,-4.580167748879537E-01,
     *          -2.816035490858031E-01,-9.501250930110497E-02,
     *           9.501250930110497E-02, 2.816035490858031E-01,
     *           4.580167748879537E-01, 6.178762411475158E-01,
     *           7.554044054509828E-01, 8.656312004236659E-01,
     *           9.445750221571226E-01, 9.894009348053501E-01/
c
      dimension w16(16)
      data w16 / 2.715245988238494E-02, 6.225352488578796E-02,
     *           9.515851275863560E-02, 1.246289719739357E-01,
     *           1.495959887529875E-01, 1.691565185319576E-01,
     *           1.826034138545381E-01, 1.894506093597724E-01,
     *           1.894506093597724E-01, 1.826034138545381E-01,
     *           1.691565185319576E-01, 1.495959887529875E-01,
     *           1.246289719739357E-01, 9.515851275863560E-02,
     *           6.225352488578796E-02, 2.715245988238494E-02/
c
      dimension q17(17)
      data q17 /-9.905754731988323E-01,-9.506755109286769E-01,
     *          -8.802391285190447E-01,-7.815139612433595E-01,
     *          -6.576711003169748E-01,-5.126904693661325E-01,
     *          -3.512317014746474E-01,-1.784841437241518E-01,
     *           0.000000000000000E+00, 1.784841437241518E-01,
     *           3.512317014746474E-01, 5.126904693661325E-01,
     *           6.576711003169748E-01, 7.815139612433595E-01,
     *           8.802391285190447E-01, 9.506755109286769E-01,
     *           9.905754731988323E-01/
c
      dimension w17(17)
      data w17 / 2.414830826353617E-02, 5.545954122636184E-02,
     *           8.503616477988524E-02, 1.118838649150996E-01,
     *           1.351363821906063E-01, 1.540457637361635E-01,
     *           1.680040872484053E-01, 1.765626726681942E-01,
     *           1.794464299434957E-01, 1.765626726681942E-01,
     *           1.680040872484053E-01, 1.540457637361635E-01,
     *           1.351363821906063E-01, 1.118838649150996E-01,
     *           8.503616477988524E-02, 5.545954122636184E-02,
     *           2.414830826353617E-02/
c
      dimension q18(18)
      data q18 /-9.915651552644213E-01,-9.558238817219664E-01,
     *          -8.926023066698663E-01,-8.037046825563874E-01,
     *          -6.916866475115471E-01,-5.597703484081864E-01,
     *          -4.117506698550099E-01,-2.518858458532465E-01,
     *          -8.477486827585889E-02, 8.477486827585889E-02,
     *           2.518858458532465E-01, 4.117506698550099E-01,
     *           5.597703484081864E-01, 6.916866475115471E-01,
     *           8.037046825563874E-01, 8.926023066698663E-01,
     *           9.558238817219664E-01, 9.915651552644213E-01/
c
      dimension w18(18)
      data w18 / 2.161604712623478E-02, 4.971462371675944E-02,
     *           7.642583732330538E-02, 1.009421665972655E-01,
     *           1.225553167380302E-01, 1.406429709609593E-01,
     *           1.546846289761835E-01, 1.642763052039739E-01,
     *           1.691421033572879E-01, 1.691421033572879E-01,
     *           1.642763052039739E-01, 1.546846289761835E-01,
     *           1.406429709609593E-01, 1.225553167380302E-01,
     *           1.009421665972655E-01, 7.642583732330538E-02,
     *           4.971462371675944E-02, 2.161604712623478E-02/
c
      dimension q19(19)
      data q19 /-9.924069252504020E-01,-9.602085692181328E-01,
     *          -9.031568746217613E-01,-8.227163080368496E-01,
     *          -7.209684928513420E-01,-6.005480751831690E-01,
     *          -4.645735505499827E-01,-3.165664003208389E-01,
     *          -1.603599369336122E-01, 0.000000000000000E+00,
     *           1.603599369336122E-01, 3.165664003208389E-01,
     *           4.645735505499827E-01, 6.005480751831690E-01,
     *           7.209684928513420E-01, 8.227163080368496E-01,
     *           9.031568746217613E-01, 9.602085692181328E-01,
     *           9.924069252504020E-01/
c
      dimension w19(19)
      data w19 / 1.946158064917772E-02, 4.481377064401640E-02,
     *           6.904390616355083E-02, 9.148932102419895E-02,
     *           1.115660506917955E-01, 1.287536835335080E-01,
     *           1.426069278403044E-01, 1.527668260608829E-01,
     *           1.589700388089103E-01, 1.610557891673102E-01,
     *           1.589700388089103E-01, 1.527668260608829E-01,
     *           1.426069278403044E-01, 1.287536835335080E-01,
     *           1.115660506917955E-01, 9.148932102419895E-02,
     *           6.904390616355083E-02, 4.481377064401640E-02,
     *           1.946158064917772E-02/
c
      dimension q20(20)
      data q20 /-9.931288558459677E-01,-9.639732440533934E-01,
     *          -9.122375019660727E-01,-8.391222238213100E-01,
     *          -7.463393329480198E-01,-6.360627202301004E-01,
     *          -5.108765051928753E-01,-3.737145022118465E-01,
     *          -2.277916295223091E-01,-7.652857473953303E-02,
     *           7.652857473953303E-02, 2.277916295223091E-01,
     *           3.737145022118465E-01, 5.108765051928753E-01,
     *           6.360627202301004E-01, 7.463393329480198E-01,
     *           8.391222238213100E-01, 9.122375019660727E-01,
     *           9.639732440533934E-01, 9.931288558459677E-01/
c
      dimension w20(20)
      data w20 / 1.761335247052325E-02, 4.059998721060333E-02,
     *           6.267002263771275E-02, 8.327448170371499E-02,
     *           1.019281252567558E-01, 1.181934022237969E-01,
     *           1.316889125381483E-01, 1.420980141551057E-01,
     *           1.491762718044848E-01, 1.527574299991544E-01,
     *           1.527574299991544E-01, 1.491762718044848E-01,
     *           1.420980141551057E-01, 1.316889125381483E-01,
     *           1.181934022237969E-01, 1.019281252567558E-01,
     *           8.327448170371499E-02, 6.267002263771275E-02,
     *           4.059998721060333E-02, 1.761335247052325E-02/
c
      save
c
c=======================================================================
c
c     Determine number of quadrature nodes needed.
c
c     For 2*nq-1 >= iorder, nq = iorder/2+1.
c
      nq = iorder/2+1
      nq = max0(nq,2)
      nq = min0(nq,20)
c     print *,'  quad1d - iorder,nq = ',iorder,nq
c
      if(nq.gt.nqs_max) then
        print *,'*** ERROR ***'
        print *,'    Not enough storage for quadrature points ',nq
        print *,'    Try increasing nqs_max or decreasing iqp.'
        stop
      endif
c
      if(nq.le.2) then
        do i=1,nq
          sq(i) = q2(i)
          wq(i) = w2(i)
        enddo
      elseif(nq.eq.3) then
        do i=1,nq
          sq(i) = q3(i)
          wq(i) = w3(i)
        enddo
      elseif(nq.eq.4) then
        do i=1,nq
          sq(i) = q4(i)
          wq(i) = w4(i)
        enddo
      elseif(nq.eq.5) then
        do i=1,nq
          sq(i) = q5(i)
          wq(i) = w5(i)
        enddo
      elseif(nq.eq.6) then
        do i=1,nq
          sq(i) = q6(i)
          wq(i) = w6(i)
        enddo
      elseif(nq.eq.7) then
        do i=1,nq
          sq(i) = q7(i)
          wq(i) = w7(i)
        enddo
      elseif(nq.eq.8) then
        do i=1,nq
          sq(i) = q8(i)
          wq(i) = w8(i)
        enddo
      elseif(nq.eq.9) then
        do i=1,nq
          sq(i) = q9(i)
          wq(i) = w9(i)
        enddo
      elseif(nq.eq.10) then
        do i=1,nq
          sq(i) = q10(i)
          wq(i) = w10(i)
        enddo
      elseif(nq.eq.11) then
        do i=1,nq
          sq(i) = q11(i)
          wq(i) = w11(i)
        enddo
      elseif(nq.eq.12) then
        do i=1,nq
          sq(i) = q12(i)
          wq(i) = w12(i)
        enddo
      elseif(nq.eq.13) then
        do i=1,nq
          sq(i) = q13(i)
          wq(i) = w13(i)
        enddo
      elseif(nq.eq.14) then
        do i=1,nq
          sq(i) = q14(i)
          wq(i) = w14(i)
        enddo
      elseif(nq.eq.15) then
        do i=1,nq
          sq(i) = q15(i)
          wq(i) = w15(i)
        enddo
      elseif(nq.eq.16) then
        do i=1,nq
          sq(i) = q16(i)
          wq(i) = w16(i)
        enddo
      elseif(nq.eq.17) then
        do i=1,nq
          sq(i) = q17(i)
          wq(i) = w17(i)
        enddo
      elseif(nq.eq.18) then
        do i=1,nq
          sq(i) = q18(i)
          wq(i) = w18(i)
        enddo
      elseif(nq.eq.19) then
        do i=1,nq
          sq(i) = q19(i)
          wq(i) = w19(i)
        enddo
      elseif(nq.ge.20) then
        do i=1,nq
          sq(i) = q20(i)
          wq(i) = w20(i)
        enddo
      endif
c
      return
      end
c
c=======================================================================
c
c     Routines for nodal basis functions on standard nodes.
c
c=======================================================================
c
c     Routines for computing basis functions & derivatives.
c
c=======================================================================
c
      subroutine getrho(n,nq,sq,rho,nd)
c
c=======================================================================
c
c     Compute standard nodal basis functions at quadrature nodes.
c
c=======================================================================
c
      implicit real*8(a-h,o-z)
c
      dimension sq (*)
      dimension rho(nd,*)
c
      dimension s  (20)
c
c     Storage for data set (computed values).
c
      parameter(nddim=10000)
      dimension dat(nddim)
      dimension i_pt(20)
      dimension i_nn(20)
      dimension i_nq(20)
      data      nds /0/
      data      kds /1/
      save
c
c-----------------------------------------------------------------------
c
c     Take derivatives of interpolating functions.
c
      if(nq.gt.20) stop 'nq > 20 in dpds'
c
c     Test for stored values.
c
      call dbfind2(i_nn,n,i_nq,nq,nds,i_pt,iptr)
c
c     Load existing values.
c
      if(iptr.gt.0) then
        call dbload2(dat,iptr,rho,nd,n,nq)
c
c     Compute values & store.
c
      else
c       h = 2.0d0/dfloat(n-1)
c
c       Set node locations.
c
c       do i=1,n
c         s(i) = -1.0d0+(i-1)*h
c       enddo
c
        call getpn1d(n,s)
c
        do i=1,n
          dd = 1.0d0
          do j=1,n
            if(j.ne.i) dd = dd*(s(i)-s(j))
          enddo
c
          do k=1,nq
            rho(i,k) = 1.0d0
            do j=1,n
              if(j.ne.i) rho(i,k) = rho(i,k)*(sq(k)-s(j))
            enddo
            rho(i,k) = rho(i,k)/dd
          enddo
        enddo
c
c       Save data.
c
        if(kds+n*nq.lt.nddim) then
          nds = nds+1
          i_nq(nds) = nq
          i_nn(nds) = n
          i_pt(nds) = kds
          call dbsave2(dat,kds,rho,nd,n,nq)
        else
          print *,'No room to store data in dphiq'
        endif
c
      endif
c
      return
      end
c
c=======================================================================
c
      subroutine getdrho(n,nq,sq,drho,nd)
c
c=======================================================================
c
c     Compute standard nodal basis functions at quadrature nodes.
c
c=======================================================================
c
      implicit real*8(a-h,o-z)
c
      dimension sq  (*)
      dimension drho(nd,*)
c
      dimension s  (20)
      dimension pds(20)
c
c     Storage for data set (computed values).
c
      parameter(nddim=10000)
      dimension dat(nddim)
      dimension i_pt(20)
      dimension i_nn(20)
      dimension i_nq(20)
      data      nds /0/
      data      kds /1/
      save
c
c-----------------------------------------------------------------------
c
c     Take derivatives of interpolating functions.
c
      if(nq.gt.20) stop 'nq > 20 in dpds'
c
c     Test for stored values.
c
      call dbfind2(i_nq,nq,i_nn,n,nds,i_pt,iptr)
c
c     Load existing values.
c
      if(iptr.gt.0) then
        call dbload2(dat,iptr,drho,nd,n,nq)
c
c     Compute values & store.
c
      else
        h = 2.0d0/dfloat(n-1)
c
c       Set node locations.
c
        do i=1,n
          s(i) = -1.0d0+(i-1)*h
          pds(i) = 1.0d0
        enddo
c
        do i=1,n
        enddo
c
c       Compute partial products.
c
        do i=1,n
          do j=i+1,n
            ds = s(i)-s(j)
            pds(i) =  pds(i)*ds
            pds(j) = -pds(j)*ds
          enddo
        enddo
c
c       Compute derivative coefficients.
c
        do i=1,n
          do j=1,nq
            drho(i,j) = 0.0d0
            do k=1,n
              if(k.ne.i) then
                dt = 1.0d0
                do l=1,n
                  if(l.ne.i.and.l.ne.k) then
                    dt = dt*(sq(j)-s(l))
                  endif
                enddo
                drho(i,j) = drho(i,j)+dt
              endif
            enddo
            drho(i,j) = drho(i,j)/pds(i)
          enddo
        enddo
c
c       Save data.
c
        if(kds+n*nq.lt.nddim) then
          nds = nds+1
          i_nq(nds) = nq
          i_nn(nds) = n
          i_pt(nds) = kds
          call dbsave2(dat,kds,drho,nd,n,nq)
        else
          print *,'No room to store data in dphiq'
        endif
      endif
c
      return
      end
c
c=======================================================================
c
      subroutine getdchi(nq,sq,dchi,nd)
c
c=======================================================================
c
c     Compute derivatives dchi(i,j) = dchi_i(s_j)/ds
c
c     These give difference formulas on quadrature nodes:
c
c     df(s_i)/fs = sum f_j*dchi(s_i)/ds = sum f_j*dchi(j,i)
c                   j                         j
c
c=======================================================================
c
      implicit real*8(a-h,o-z)
c
      dimension sq (*)
      dimension dchi(nd,*)
      dimension ds (20,20)
      dimension pds(20)
c
c     Storage for data set (computed values).
c
      parameter(nddim=10000)
      dimension dat(nddim)
      dimension i_pt(20)
      dimension i_nq(20)
      data      nds /0/
      data      kds /1/
      save
c
c-----------------------------------------------------------------------
c
c     Take derivatives of interpolating functions.
c
      if(nq.gt.20) stop 'nq > 20 in dpds'
c
c     Test for stored values.
c
      call dbfind1(i_nq,nq,nds,i_pt,iptr)
c
c     Load existing values.
c
      if(iptr.gt.0) then
        call dbload2(dat,iptr,dchi,nd,nq,nq)
c
c     Compute values & store.
c
      else
        do i=1,nq
          pds(i) = 1.0d0
        enddo
c
c       Compute partial products.
c
        do i=1,nq
          do j=i+1,nq
            ds(i,j) = sq(i)-sq(j)
            ds(j,i) = -ds(i,j)
            pds(i) = pds(i)*ds(i,j)
            pds(j) = pds(j)*ds(j,i)
          enddo
        enddo
c
c       Compute derivative coefficients.
c
        do i=1,nq
          dchi(i,i) = 0.0d0
          do j=1,nq
            ii = ii+1
            if(j.ne.i) then
              dchi(i,j) = pds(j)/(pds(i)*ds(j,i))
              dchi(i,i) = dchi(i,i)+1.0d0/ds(i,j)
            endif
          enddo
        enddo
c
c       Save data.
c
        if(kds+nq*nq.lt.nddim) then
          nds = nds+1
          i_nq(nds) = nq
          i_pt(nds) = kds
          call dbsave2(dat,kds,dchi,nd,nq,nq)
        else
          print *,'No room to store data in dphiq'
        endif
c
      endif
c
      return
      end
c
c=======================================================================
c
      subroutine getphinq(rho,ndr,n1,nq1,phinq,ndp,nn,nq)
c
c=======================================================================
c
c     Compute values of 2D nodal basis functions at quad points.
c
c     Store indexed by nn & nq.
c
c=======================================================================
c
      implicit real*8(a-h,o-z)
c
      dimension rho (ndr,*)
      dimension phinq(ndp,*)
c
c     Storage for data set (computed values).
c
      parameter(nddim=50000)
      dimension dat(nddim)
      dimension i_pt(20)
      dimension i_nn(20)
      dimension i_nq(20)
      data      nds /0/
      data      kds /1/
      save
c
c-----------------------------------------------------------------------
c
c     Test for stored values.
c
      call dbfind2(i_nq,nq,i_nn,nn,nds,i_pt,iptr)
c
c     Load existing values.
c
      if(iptr.gt.0) then
        call dbload2(dat,iptr,phinq,ndp,nn,nq)
c
c     Compute values & store.
c
      else
        in = 0
        do j1=1,n1
          do i1=1,n1
            in = in+1
            iq = 0
            do j2=1,nq1
              do i2=1,nq1
                iq = iq+1
                phinq(in,iq) = rho(i1,i2)*rho(j1,j2)
              enddo
            enddo
          enddo
        enddo
c
c       Save data.
c
        if(kds+nn*nq.lt.nddim) then
          nds = nds+1
          i_nq(nds) = nq
          i_nn(nds) = nn
          i_pt(nds) = kds
          call dbsave2(dat,kds,phinq,ndp,nn,nq)
        else
          print *,'No room to store data in getphinq'
        endif
      endif
c
      return
      end
