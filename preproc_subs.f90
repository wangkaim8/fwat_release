module preproc_subs 
    use specfem_par, only: cr => CUSTOM_REAL, MAX_LEN => MAX_STRING_LEN,&
                                network_name, station_name,nrec,myrank,&
                                t0,NSTEP,DT,IIN,OUTPUT_FILES
    use fullwave_adjoint_tomo_par
    use measure_adj_mod
    use sacio

    implicit none
    integer                                                :: ier,irec
contains

subroutine read_fktimes(ievt, ttp, tb, te)
    character(len=MAX_LEN)                             :: datafile, dummystring
    character(len=256),dimension(:), allocatable       :: fk_netwk,fk_stnm
    real(kind=4) , dimension(:), allocatable           :: fk_ttp,fk_tb,fk_te 
    integer                                            :: fk_irec,fk_nrec,ievt
    real                                               :: dummy
    real(kind=4)                                       :: ttp(nrec),tb(nrec),te(nrec)

    if (myrank == 0) then
        datafile='src_rec/FKtimes_'//trim(evtid_names(ievt))
        open(unit=IIN,file=trim(datafile),status='old',action='read',iostat=ier)
        if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(datafile))
            ! reads all stations
        fk_nrec=0
        do
            read(IIN,*,iostat=ier) dummystring,dummystring,dummy
            if (ier /= 0) exit
            fk_nrec=fk_nrec+1
        enddo
        ! close receiver file
        close(IIN)
        allocate(fk_netwk(fk_nrec)) 
        allocate(fk_stnm(fk_nrec)) 
        allocate(fk_ttp(fk_nrec)) 
        allocate(fk_tb(fk_nrec)) 
        allocate(fk_te(fk_nrec)) 
        open(unit=IIN,file=trim(datafile),status='old',action='read',iostat=ier)
        do irec=1,fk_nrec
            if (TW_AFTER.eq.0.) then
                if (irec==1) write(*,*) 'Read TW_BEFORE and TW_AFTER from FKtimes...'
                read(IIN,*,iostat=ier) fk_netwk(irec),fk_stnm(irec),fk_ttp(irec),fk_tb(irec),fk_te(irec)
            else
                if (irec==1) write(*,*) 'Read TW_BEFORE and TW_AFTER from FWAT.PAR'
                read(IIN,*,iostat=ier) fk_netwk(irec),fk_stnm(irec),fk_ttp(irec)
            endif
        enddo
        close(IIN)
        ! find the right ttp among all stations (Note STATIONS_FILTERED might be
        ! less than STATIONS/FKtimes, thus we need to search for the right ttp )
        do irec=1,nrec
            do fk_irec=1,fk_nrec
                if (trim(network_name(irec))==trim(fk_netwk(fk_irec)) .and. &
                    trim(station_name(irec))==trim(fk_stnm(fk_irec)) ) then
                    ttp(irec)=fk_ttp(fk_irec)
                    if (TW_AFTER.eq.0.) then
                        tb(irec)=fk_tb(fk_irec)
                        te(irec)=fk_te(fk_irec)
                    else
                        tb(irec)=TW_BEFORE
                        te(irec)=TW_AFTER
                    endif
                    write(*,*)'netwk,stnm,ttp,tb,te=',trim(network_name(irec)),'.',trim(station_name(irec)), &
                                                        ttp(irec),tb(irec),te(irec)
                    if( (ttp(irec)-t0+te(irec)) > (-t0+(NSTEP-1)*DT) ) then
                        write(*,*)'ttp exceed data range'
                        stop
                    endif
                endif
            enddo
        enddo
        deallocate(fk_netwk)
        deallocate(fk_stnm)
        deallocate(fk_ttp)
        deallocate(fk_tb)
        deallocate(fk_te)
    endif
    call bcast_all_cr(ttp,nrec)
    call bcast_all_cr(tb,nrec)
    call bcast_all_cr(te,nrec)
end subroutine

subroutine pre_proc_tele_elastic(ievt, irec, glob_sem_disp, win_tb,win_te, fstart0, fend0, bandname, &
                                 baz_all, glob_stnm, glob_dat_tw, glob_syn_tw, glob_ff, ttp)

    use telestf_mod

    integer                                          :: icomp, ievt, irec, NDIM_CUT
    logical                                          :: findfile
    double precision                                 :: t01,dt1, fstart0, fend0
    integer                                          :: yr,jda,ho,mi, npt1
    double precision                                 :: t0_inp,t1_inp,dt_inp
    double precision, dimension(NDIM)                :: datarray
    double precision, dimension(NDIM)                :: dat_inp, syn_inp
    double precision, dimension(NDIM)                :: dat_inp_bp, syn_inp_bp
    real(kind=cr), dimension(nrec,NSTEP,3)           :: glob_sem_disp
    real(kind=cr)                                    :: seismo_syn(3,NSTEP)                           
    real(kind=4), dimension(NSTEP)                   :: one_seismo_dat,one_seismo_syn
    double precision                                 :: sec,dist,az,baz,slat,slon
    character(len=10)                                :: net,sta,chan_dat
    character(len=MAX_LEN)                           :: bandname,adjfile,datafile
    double precision                                 :: baz_all(nrec)
    real                                             :: win_tb,win_te
    character(len=256), dimension(nrec)              :: glob_stnm
    real(kind=4), dimension(nrec,NSTEP,NRCOMP), intent(inout)       :: glob_dat_tw,glob_syn_tw
    real(kind=4), dimension(nrec,NSTEP)              :: glob_ff
    real(kind=4)                                     :: ttp(nrec)

    seismo_syn(:,:)=0.d0 !glob_sem_disp(irec,:,:)

    do icomp=1,NRCOMP
        datafile='./'//trim(in_dat_path(ievt))//'/'//trim(network_name(irec))//'.'&
                //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))//'.sac' 
        !write(*,*)'myrank,datafile=',myrank,trim(datafile)
        inquire(file=trim(datafile),exist=findfile)
        if ( .not. findfile ) then
            stop 'No such sac file in the fwat_data'
        else
            call drsac1(trim(datafile),datarray,npt1,t01,dt1)
            call get_sacfile_header(trim(datafile),yr,jda,ho,mi,sec,net,sta, &
                                    chan_dat,dist,az,baz,slat,slon)
            
            baz_all(irec)=baz
            !write(*,*)'myrank,datafile,npt1,t01,dt1,NSTEP,t0,DT,dist,ttp= ',myrank,trim(datafile),npt1,t01,dt1,&
            !                                   NSTEP,t0,DT,dist,ttp(irec)
            if(abs(dt1-dT)>0.0001) stop 'delta of data does NOT match DT of SEM synthetics'
            ! only do the rotation for one time
            if (icomp.eq.1 ) then
                if (trim(dat_coord)=='ZRT') then
                    call rotate_ZNE_to_ZRT(glob_sem_disp(irec,:,3),glob_sem_disp(irec,:,2),&
                    glob_sem_disp(irec,:,1),seismo_syn(1,:),seismo_syn(2,:),seismo_syn(3,:),NSTEP,real(baz)) 
                else
                    seismo_syn(1,:)=glob_sem_disp(irec,:,3)
                    seismo_syn(2,:)=glob_sem_disp(irec,:,2)
                    seismo_syn(3,:)=glob_sem_disp(irec,:,1)
                endif
            endif
            !****************  pre-process dat and syn ****************************
            dt_inp=DT !
            t0_inp=-t0 
            t1_inp=-t0+(NSTEP-1)*DT
            NDIM_CUT=NSTEP !(t1_inp-t0_inp)/dt_inp+1
            dat_inp(:)=0.
            syn_inp(:)=0.
            dat_inp(1:npt1)=datarray(1:npt1)
            syn_inp(1:NSTEP)=dble(seismo_syn(icomp,1:NSTEP))
            !!! Filter
            ! rtrend
            !call detrend(dat_inp,NDIM_CUT,t0_inp,dt_inp)
            !call detrend(dat_inp,NDIM_CUT,t0_inp,dt_inp)
            call detrend(dat_inp,NDIM_CUT)
            call detrend(syn_inp,NDIM_CUT)
            ! remean
            dat_inp(1:NDIM_CUT)=dat_inp(1:NDIM_CUT)-sum(dat_inp(1:NDIM_CUT))/NDIM_CUT
            syn_inp(1:NDIM_CUT)=syn_inp(1:NDIM_CUT)-sum(syn_inp(1:NDIM_CUT))/NDIM_CUT
            dat_inp_bp=dat_inp
            syn_inp_bp=syn_inp
            call bandpass(dat_inp_bp,NDIM_CUT,dt_inp,fstart0,fend0)
            call bandpass(syn_inp_bp,NDIM_CUT,dt_inp,fstart0,fend0)
            ! rtrend
            !call detrend(dat_inp_bp,NDIM_CUT,t0_inp,dt_inp)
            !call detrend(syn_inp_bp,NDIM_CUT,t0_inp,dt_inp)
            call detrend(dat_inp_bp,NDIM_CUT)
            call detrend(syn_inp_bp,NDIM_CUT)
            ! remean
            dat_inp_bp(1:NDIM_CUT)=dat_inp_bp(1:NDIM_CUT)-sum(dat_inp_bp(1:NDIM_CUT))/NDIM_CUT
            syn_inp_bp(1:NDIM_CUT)=syn_inp_bp(1:NDIM_CUT)-sum(syn_inp_bp(1:NDIM_CUT))/NDIM_CUT

            if (VERBOSE_MODE) then
                adjfile=trim(OUTPUT_FILES)//'/syn.'//trim(network_name(irec))//'.'&
                        //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
                        //'.sac'//'.'//trim(bandname) 
                call dwsac1(trim(adjfile),syn_inp_bp,NDIM_CUT,t0_inp,dt_inp)
                adjfile=trim(OUTPUT_FILES)//'/dat.'//trim(network_name(irec))//'.'&
                        //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
                        //'.sac'//'.'//trim(bandname) 
                call dwsac1(trim(adjfile),dat_inp_bp,NDIM_CUT,t0_inp,dt_inp)
            endif

            one_seismo_dat(1:NDIM_CUT)=dat_inp_bp(1:NDIM_CUT)
            glob_dat_tw(irec,1:NDIM_CUT,icomp)=dat_inp_bp(1:NDIM_CUT)
            one_seismo_syn(1:NDIM_CUT)=syn_inp_bp(1:NDIM_CUT)
            glob_syn_tw(irec,1:NDIM_CUT,icomp)=syn_inp_bp(1:NDIM_CUT)
            glob_stnm(irec)=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
                            //trim(station_name(irec))//'.'//trim(CH_CODE)& 
                            //'.sac'//'.'//trim(bandname)
            if (icomp==1)  then !!! get only vertical deconved traces
            write(*,*)'run time_deconv for ',trim(glob_stnm(irec))//' on window: ',win_tb,win_te
            call time_iterdeconv(glob_stnm,one_seismo_dat,one_seismo_syn,glob_dat_tw,glob_syn_tw,glob_ff,&
                                 ttp(irec),win_tb,win_te,irec,nrec,NSTEP,real(-t0),real(DT),fstart0,fend0)
            endif
        endif ! end findfile
    enddo ! end icomp  
end subroutine pre_proc_tele_elastic


subroutine average_amp_scale(glob_dat_tw, icomp, avgamp)
    real, intent(inout)                              :: avgamp
    real                                             :: avgamp0
    integer                                          :: igood, icomp
    real(kind=4), dimension(nrec,NSTEP,NRCOMP)       :: glob_dat_tw

    ! use only Z component for amplitude scale
    avgamp0=0.
    do irec =1 ,nrec
        avgamp0=avgamp0+maxval(abs(glob_dat_tw(irec,:,icomp))) 
    enddo
    avgamp0=avgamp0/nrec
    avgamp=0
    igood=0
    do irec =1, nrec
        if ((maxval(abs(glob_dat_tw(irec,:,icomp)))-avgamp0)<0.2*avgamp0) then
            avgamp=avgamp+maxval(abs(glob_dat_tw(irec,:,icomp)))
            igood=igood+1
        endif
    enddo
    avgamp=avgamp/igood
end subroutine average_amp_scale

end module
