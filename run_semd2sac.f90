!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
!
subroutine run_semd2sac(ievt,simu_type)

  use specfem_par, only: CUSTOM_REAL, MAX_STRING_LEN, IIN, network_name, station_name,nrec, nrec_local, &
       seismograms_d, ispec_selected_rec, number_receiver_global, myrank, t0
  use specfem_par_coupling, only: phi_FK
  use shared_input_parameters, only: NSTEP, DT, SUPPRESS_UTM_PROJECTION

  use specfem_par_elastic, only: ispec_is_elastic, ELASTIC_SIMULATION
  use specfem_par_acoustic, only: ispec_is_acoustic, ACOUSTIC_SIMULATION
  use fullwave_adjoint_tomo_par
  use my_mpi
  use measure_adj_mod
  use telestf_mod
  
  implicit none
  character(len=MAX_STRING_LEN)                             :: simu_type
  integer                                                   :: irec, irec_local, ispec
  integer                                                   :: ievt,icomp,ier,igaus
  character(len=MAX_STRING_LEN)                             :: datafile
  integer                                                   :: npt1
  double precision, dimension(NDIM)                         :: datarray
  real(kind=4), dimension(NSTEP)                            :: stf_array
  double precision, dimension(NSTEP)                        :: uin, win, rfi
  real(kind=4), dimension(:), allocatable                   :: tmpl
  double precision                                          :: t01,dt1
  ! for rotation
  real(kind=CUSTOM_REAL), dimension(3,NSTEP)                :: seismo_syn
  ! for reading elon,elat,stlon,stlat from CMTSOLUTION and STATIONS
  double precision                                          :: elat,elon,edep
  double precision                                          :: geo_src_r,geo_src_lon,geo_src_lat
  double precision                                          :: geo_sta_r,geo_sta_lon,geo_sta_lat
  double precision                                          :: src_x,src_y,src_z,sta_x,sta_y,sta_z
  double precision, dimension(nrec)                         :: stlat,stlon,stele,stbur                                           
  character(len=MAX_STRING_LEN)                             :: src_filename, rec_filename,string,bandname
  double precision                                          :: gcarc, dist, baz, az
  ! for covert xyz coords to geo coords
  double precision                                          :: lon_center_chunk, lat_center_chunk, chunk_azi
  double precision, dimension(3,3)                          :: rotation_matrix

  lon_center_chunk=0.
  lat_center_chunk=0.
  chunk_azi=0.

  !
  if (myrank==0) write(OUT_FWAT_LOG,*) 'This is run_semd2sac ...' 
  if (myrank==0 ) then
    if (simu_type=='tele' .or. simu_type=='leq') then
       if (simu_type=='leq') then
          src_filename='src_rec/CMTSOLUTION_'//trim(evtid_names(ievt))
       else
          src_filename='src_rec/CMTSOLUTION'
       endif
       write(*,*)'src_filename=',trim(src_filename)
       open(unit = 123,file=trim(src_filename),status='old',action='read',iostat=ier)
       if (ier /= 0) then
          print *,'Error opening file: ',trim(src_filename)
          stop 'Error opening CMTSOLUTION file'
       endif
       read(123,"(a256)",iostat=ier) string
       ! skips empty lines
       do while (len_trim(string) == 0)
         read(123,"(a256)",iostat=ier) string
         if (ier /= 0) then
           print *, 'Error reading header line in source ',ievt
           stop 'Error reading header line in station in CMTSOLUTION file'
         endif
       enddo
       ! ignore line with event name
       read(123,"(a)",iostat=ier) string
       ! read time shift
       read(123,"(a)",iostat=ier) string
       ! read half duration
       read(123,"(a)",iostat=ier) string
       ! read latitude
       read(123,"(a)",iostat=ier) string
       read(string(10:len_trim(string)),*) elat
       ! read longitude
       read(123,"(a)",iostat=ier) string
       read(string(11:len_trim(string)),*) elon
       ! read depth
       read(123,"(a)",iostat=ier) string
       read(string(7:len_trim(string)),*) edep
       close(123)
    elseif (simu_type == 'noise') then
        src_filename='src_rec/FORCESOLUTION_'//trim(evtid_names(ievt))
       write(*,*)'src_filename=',trim(src_filename)
       open(unit = 123,file=trim(src_filename),status='old',action='read',iostat=ier)
       if (ier /= 0) then
          print *,'Error opening file: ',trim(src_filename)
          stop 'Error opening FORCESOLUTION file'
       endif
       ! ignore line with event name
       read(123,"(a)",iostat=ier) string
       ! read time shift
       read(123,"(a)",iostat=ier) string
       ! read half duration
       read(123,"(a)",iostat=ier) string
       ! read latitude
       read(123,"(a)",iostat=ier) string
       read(string(10:len_trim(string)),*) elat
       ! read longitude
       read(123,"(a)",iostat=ier) string
       read(string(11:len_trim(string)),*) elon
       ! read depth
       read(123,"(a)",iostat=ier) string
       read(string(7:len_trim(string)),*) edep
       close(123)
    else
      elat = -12345.
      elon = -12345.
      edep = -12345.
    endif 
   !  write(*,*)'elon,elat,edep=',elon,elat,edep 

    ! opens STATIONS or STATIONS_ADJOINT file
    rec_filename='src_rec/STATIONS_'//trim(evtid_names(ievt))//'_FILTERED'
    open(unit=IIN,file=trim(rec_filename),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))
    ! reads all stations
    do irec = 1,nrec
      read(IIN,*,iostat=ier) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
      if (ier /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))
    enddo
    ! close receiver file
    close(IIN)
    
    if (simu_type=='tele') then
      ! open STF file
      datafile='src_rec/STF_'//trim(evtid_names(ievt))//'.sac'
      call drsac1(trim(datafile),datarray,npt1,t01,dt1) 
      if(npt1.ne.NSTEP) then
        write(*,*) 'npts of STF not match NSTEP of specfem3d!!!'
        stop
      endif
      stf_array(1:NSTEP)=datarray(1:NSTEP)
    endif
  endif  ! end if of myrank
  ! broadcast 
  call bcast_all_singledp(elat)
  call bcast_all_singledp(elon)
  call bcast_all_singledp(edep)
  call bcast_all_dp(stlat,nrec)
  call bcast_all_dp(stlon,nrec)
  call bcast_all_dp(stele,nrec)
  call bcast_all_dp(stbur,nrec)
  if (simu_type=='tele') then
    call bcast_all_cr(stf_array,NSTEP)
  endif
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%p%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(nrec_local > 0) then
    !write(*,*) "Rank ",myrank,' has ',nrec_local,' receivers:'
    do irec_local = 1, nrec_local
       irec = number_receiver_global(irec_local)
       ispec = ispec_selected_rec(irec)

       if (ELASTIC_SIMULATION) then
          if (ispec_is_elastic(ispec)) then
            !====================================================================================================
            seismo_syn(:,:)=seismograms_d(:,irec_local,1:NSTEP) 
            if (SUPPRESS_UTM_PROJECTION) then
               ! call calc_dist_baz_cart(elon,elat,stlon(irec),stlat(irec),dist,baz)
               ! gcarc=dist/110000.
               ! write(*,*)'Cartesian: net.stnm,elat,elon,edep,stlat,stlon,gcarc,dist,baz=',trim(network_name(irec))//'.'//trim(station_name(irec)),&
               ! elat,elon,edep,stlat(irec),stlon(irec),gcarc,dist,baz
               src_x=elon
               src_y=elat
               src_z=6371000.-edep*1000.
               sta_x=stlon(irec)
               sta_y=stlat(irec)
               sta_z=6371000. + stele(irec) - stbur(irec)   

               call compute_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)
               call cart2geogr(src_x,src_y,src_z,rotation_matrix,geo_src_lon,geo_src_lat,geo_src_r)
               call cart2geogr(sta_x,sta_y,sta_z,rotation_matrix,geo_sta_lon,geo_sta_lat,geo_sta_r)
            else
               geo_src_lon=elon
               geo_src_lat=elat
               geo_sta_lon=stlon(irec)
               geo_sta_lat=stlat(irec)     
            endif
            call calc_delta_dist_baz(geo_src_lat,geo_src_lon,geo_sta_lat,geo_sta_lon,gcarc,dist,baz)
            dist=dist/1000.
            if (baz>0) then
               az=360.-baz !! az is not used in measure_adj, only baz used
            else
               az=360.+baz 
            endif
            !write(*,*)'net.stnm,elat,elon,edep,stlat,stlon,gcarc,dist,baz=',trim(network_name(irec))//'.'//trim(station_name(irec)),&
            !geo_src_lat,geo_src_lon,edep,geo_sta_lat,geo_sta_lon,gcarc,dist,baz
            do icomp=1,NRCOMP
               if (icomp.eq.1 ) then
                  if (trim(dat_coord)=='ZRT') then
                     call rotate_ZNE_to_ZRT(seismograms_d(3,irec_local,:),seismograms_d(2,irec_local,:),&
                     seismograms_d(1,irec_local,:),seismo_syn(1,:),seismo_syn(2,:),seismo_syn(3,:),NSTEP,real(baz)) 
                  else
                     seismo_syn(1,:)=seismograms_d(3,irec_local,:)
                     seismo_syn(2,:)=seismograms_d(2,irec_local,:)
                     seismo_syn(3,:)=seismograms_d(1,irec_local,:)
                  endif
               endif
               if (simu_type=='tele') then
                  ! convolve with STF for TeleFWI
                  allocate(tmpl(NSTEP))
                  call myconvolution(seismo_syn(icomp,:),stf_array(:),NSTEP,NSTEP,tmpl,0) 
                  seismo_syn(icomp,:)=tmpl*DT
               endif
               datafile='./'//trim(in_dat_path(ievt))//'/'//trim(network_name(irec))//'.'&
                        //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))//'.sac' 
               call dwsac0(trim(datafile),dble(seismo_syn(icomp,:)),NSTEP,-dble(t0),dble(DT),&
                           trim(network_name(irec)),trim(station_name(irec)),trim(CH_CODE)//trim(RCOMPS(icomp)),&
                           dist,az,baz,geo_sta_lat,geo_sta_lon,geo_src_lat,geo_src_lon,edep)

               if (simu_type=='tele') then
                  deallocate(tmpl)
               endif
            enddo ! end icomp  
          !====================================================================================================
          endif
       endif


       if (ACOUSTIC_SIMULATION) then
         if (ispec_is_acoustic(ispec)) then

          !! compute adjoint source according to cost L2 function
            write(*,*) 'To be added '
            stop
         endif
       endif
    enddo ! end loop irec
  endif ! end nrec >0

  call synchronize_all()
  if(myrank==0) write(*,*) 'Finished semd2sac here.'
  
end subroutine run_semd2sac


!================================================================================
! subroutines from specfem3D

  subroutine calc_delta_dist_baz(slati,sloni,rlati,rloni,gcarc,dist,baz)

   double precision, parameter :: deg2rad = 0.017453292519943
   double precision, parameter :: rad2deg = 57.295779513082323

    double precision, intent(in)  :: slati, sloni, rlati, rloni
    double precision, intent(out) :: gcarc, dist, baz

    double precision :: slat, slon, rlat, rlon, den, num, cosarc
    double precision :: r=6371000.

    !*** Pass in radian
    slat = deg2rad * slati
    slon = deg2rad * sloni
    rlat = deg2rad * rlati
    rlon = deg2rad * rloni

    !*** Compute distance in degree and in km
    cosarc = sin(rlat)*sin(slat) + cos(rlat)*cos(slat)*cos(rlon-slon)
    gcarc  = acos(cosarc)
    dist   = r * gcarc

    !*** Compute back-azimuth
    num = sin(rlon-slon) * cos(slat) * cos(rlat)
    den =     sin(slat) - cos(gcarc) * sin(rlat)
    baz= -(rad2deg * atan2(num,den)) - 360.
    if (baz >= 360.) then
       baz = baz -360.
    else if (baz < 0.) then
       baz = baz +360.
    else
       baz = baz
    endif

  end subroutine calc_delta_dist_baz

  ! ================================================================================
  ! Compute offset and baz from stations/source informations
  subroutine calc_dist_baz_cart(src_x,src_y,sta_x,sta_y,dist,baz)

   double precision, parameter :: deg2rad = 0.017453292519943
   double precision, parameter :: rad2deg = 57.295779513082323

   double precision, intent(in)  :: src_x, src_y, sta_x, sta_y
   double precision, intent(out) :: dist, baz

   double precision :: dx, dy

   !*** Compute distance in degree and in km
   dx   = src_x-sta_x
   dy   = src_y-sta_y
   dist = sqrt(dx*dx + dy*dy)

   !*** Compute back-azimuth
   baz  = rad2deg * atan2(dx,dy)
   if (baz >= 360.) then
      baz = baz -360.
   else if (baz < 0.) then
      baz = baz +360.
   else
      baz = baz
   endif

 end subroutine calc_dist_baz_cart

 !=======================================================================================================
 !!!!!!!!!!!!!!!! subroutines from axisem_couple_specfem3d/UTILS_COUPLING_SpecFEM !!!!!!!
 !=======================================================================================================
!!!
!!
 subroutine cart2geogr(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius)
   !
       implicit none
   
       integer NDIM !,NGLLX,NGLLY,NGLLZ
       double precision xstore,ystore,zstore
       double precision longitud,latitud,radius
       double precision rotation_matrix(3,3)
       double precision vector_ori(3),vector_rotated(3)
       double precision rayon,x,y,z,deg2rad,long,lati
       integer i,j !,k
   
       deg2rad=3.141592653589793d0/180.d0
   
       NDIM=3
   
       vector_ori(1)=xstore
       vector_ori(2)=ystore
       vector_ori(3)=zstore
   
       do i = 1,NDIM  !
          vector_rotated(i) = 0.d0
          do j = 1,NDIM
   
             vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
   
          enddo
       enddo
   
       x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
       rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)
   
       long=atan2(y,x)
       lati=asin(z/rayon)
   
       longitud = long/deg2rad
       latitud  = lati/deg2rad
       radius   = rayon/1000.d0
   
   
   
   
     end subroutine cart2geogr
   
   
   !=======================================================================================================
   !!!
   !!
     subroutine geogr2cart(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius)
   !
       implicit none
   
       integer NDIM
       double precision xstore,ystore,zstore
       double precision longitud,latitud,radius
       double precision rotation_matrix(3,3)
       double precision vector_ori(3),vector_rotated(3)
       double precision deg2rad !rayon,x,y,z,deg2rad,long,lati
       integer i,j !,k
   
       deg2rad=3.141592653589793d0/180.d0
       NDIM=3
   
       vector_ori(1)=1000.d0*radius*cos(deg2rad*longitud)*sin(deg2rad*(90.d0-latitud))
       vector_ori(2)=1000.d0*radius*sin(deg2rad*longitud)*sin(deg2rad*(90.d0-latitud))
       vector_ori(3)=1000.d0*radius*cos(deg2rad*(90.d0-latitud))
   
       do i = 1,NDIM
          vector_rotated(i) = 0.d0
          do j = 1,NDIM
   
             vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
   
          enddo
       enddo
   
       xstore=vector_rotated(1)
       ystore=vector_rotated(2)
       zstore=vector_rotated(3)
   
     end subroutine geogr2cart
   
   
   !=======================================================================================================
   !
   !=======================================================================================================
   
     subroutine  compute_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)
   
     implicit none
   
     double precision rotation_matrix(3,3),lon_center_chunk,lat_center_chunk, chunk_azi
     double precision R0(3,3),R1(3,3),R2(3,3),axe_rotation(3),R00(3,3)
   
     ! je met le chunk en 0,0
     axe_rotation(1)=0.d0; axe_rotation(2)=1.d0; axe_rotation(3)=0.d0
     call rotation_matrix_axe(R00,axe_rotation,90.d0)  ! je ramene le chunk en (0,0)
     ! rotation de l'azimuth du chunk
     axe_rotation(1)=1.d0; axe_rotation(2)=0.d0; axe_rotation(3)=0.d0
     call rotation_matrix_axe(R0,axe_rotation,90.-chunk_azi)
     ! on met le chunk a la bonne latitude
     axe_rotation(1)=0.d0; axe_rotation(2)=-1.d0; axe_rotation(3)=0.d0
     call rotation_matrix_axe(R1,axe_rotation,lat_center_chunk)
     ! on met le chunk a la bonne longitude
     axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=1.d0
     call rotation_matrix_axe(R2,axe_rotation, lon_center_chunk)
     ! rotation resultante
     call compose4matrix(rotation_matrix,R00,R0,R1,R2)
   
     end subroutine compute_rotation_matrix
   
   !=======================================================================================================
   
     subroutine  compute_inv_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)
   
     implicit none
   
     double precision rotation_matrix(3,3),lon_center_chunk,lat_center_chunk, chunk_azi
     double precision R0(3,3),R1(3,3),R2(3,3),axe_rotation(3),R00(3,3)
   
     ! on met le chunk a la bonne longitude
     axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=-1.d0
     call rotation_matrix_axe(R00,axe_rotation, lon_center_chunk)
     ! on met le chunk a la bonne latitude
     axe_rotation(1)=0.d0; axe_rotation(2)=1.d0; axe_rotation(3)=0.d0
     call rotation_matrix_axe(R0,axe_rotation,lat_center_chunk)
     ! rotation de l'azimuth du chunk
     axe_rotation(1)=-1.d0; axe_rotation(2)=0.d0; axe_rotation(3)=0.d0
     call rotation_matrix_axe(R1,axe_rotation,90.-chunk_azi)
   ! je met le chunk en 0,0
     axe_rotation(1)=0.d0; axe_rotation(2)=-1.d0; axe_rotation(3)=0.d0
     call rotation_matrix_axe(R2,axe_rotation,90.d0)  ! je ramene le chunk en (0,0)
     ! rotation resultante
     call compose4matrix(rotation_matrix,R00,R0,R1,R2)
   
   end subroutine compute_inv_rotation_matrix
   
   
   !=======================================================================================================
   !
   !   ROUTINES POUR FAIRE DES ROTATIONS 3D ET DIVERS CHANGEMENTS DE REPERES
   !
   ! Vadim Monteiller Mars 2013
   !
   !-------------------------------------------------------------------------------
   ! matrice de rotation 3D d'axe "axe" et d'angle theta (en degres)
   ! cette matrice est en complexe
   !
   !=======================================================================================================
   !
     subroutine rotation_matrix_axe(R,axe,theta)
   
     implicit none
   
     double precision axe(3),theta,pi,deg2rad
     double precision R(3,3)
     double precision c,s,ux,uy,uz,norme_axe
   
     pi=3.1415926535897932d0
     deg2rad = pi / 180.d0
     ! on normalise l'axe
     norme_axe=dsqrt(axe(1)**2 + axe(2)**2 + axe(3)**2)
   
     ! composantes de l'axe
     ux=axe(1)/norme_axe
     uy=axe(2)/norme_axe
     uz=axe(3)/norme_axe
   
     ! on calcule le cos et sin
     c=dcos(deg2rad * theta);s=dsin(deg2rad * theta)
   
     ! matrice de rotation complexe
     R(1,1)=(ux**2 + (1.d0-ux**2)*c)
     R(1,2)=(ux*uy*(1.d0-c)-uz*s)
     R(1,3)=(ux*uy*(1.d0-c)+uy*s)
   
     R(2,1)=(ux*uy*(1.d0-c)+uz*s)
     R(2,2)=(uy**2+(1.d0-uy**2)*c)
     R(2,3)=(uy*uz*(1.d0-c)-ux*s)
   
     R(3,1)=(ux*uz*(1.d0-c)-uy*s)
     R(3,2)=(uy*uz*(1.d0-c)+ux*s)
     R(3,3)=(uz**2+(1.d0-uz**2)*c)
   
     !write(49,*) ' MATRICE ROTATION '
     !write(49,*) R(1,:)
     !write(49,*) R(2,:)
     !write(49,*) R(3,:)
     !write(49,*)
   
     end subroutine rotation_matrix_axe
   
   !=======================================================================================================
   
   !=======================================================================================================
   !
   ! R=R2*R1*R0
   !
   !=======================================================================================================
   
     subroutine compose4matrix(R,R00,R0,R1,R2)
   
     implicit none
   
     double precision R(3,3),R0(3,3),R1(3,3),R2(3,3),R00(3,3),Rtmp(3,3)
     integer i,j,k
   
   
     R(:,:)=0.d0
     ! multiplication R=R0*R00
     do j=1,3
        do i=1,3
           do k=1,3
              R(i,j)=R(i,j) + R0(i,k)*R00(k,j)
           enddo
        enddo
     enddo
   
     ! multiplication R=R1*R
     Rtmp=R
     R(:,:)=0.d0
     do j=1,3
        do i=1,3
           do k=1,3
              R(i,j)=R(i,j) + R1(i,k)*Rtmp(k,j)
           enddo
        enddo
     enddo
   
     ! multiplication R=R2*R
     Rtmp=R
     R(:,:)=0.d0
     do j=1,3
        do i=1,3
           do k=1,3
              R(i,j)=R(i,j) + R2(i,k)*Rtmp(k,j)
           enddo
        enddo
     enddo
   
     !write(49,*) ' MATRICE ROTATION COMPLETE '
     !write(49,*) R(1,:)
     !write(49,*) R(2,:)
     !write(49,*) R(3,:)
     !write(49,*)
   
     end subroutine compose4matrix
