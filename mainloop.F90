!> \file mainloop.F90

!> \brief the main program loop
!! 
!! Loops through the snapshots calling the raytracing and output routines 
!<

module mainloop_mod
  use myf03_mod

  ! routines
  use gadget_general_class, only: gadget_constants_type
  use main_input_mod, only: readin_snapshot
  use oct_tree_mod, only: buildtree, setparticleorder
  use ray_mod
  use raylist_mod
  use ion_temperature_update, only: update_raylist
  use mt19937_mod, only: genrand_real1
  use output_mod, only: output_total_snap, ion_frac_out
  
  ! variables
  use global_mod, only: psys
  use global_mod, only: tree
  use global_mod, only: GV
  use global_mod, only: PLAN
  use global_mod, only: active_rays


  
  implicit none
  
  integer(i8b), parameter :: one = 1

  integer(i4b), parameter :: rays_per_leaf = 2
  integer(i4b), parameter :: rays_per_dt = 1
contains
  
  !> this is the main driver of SPHRAY
  !======================================
  subroutine mainloop()
    implicit none
    integer tid, omp_get_thread_num
    type(raylist_type) :: raylist       !< ray/particle intersections
    character(clen), parameter :: myname="mainloop"
    logical, parameter :: crash=.true.
    integer, parameter :: verb=1
    character(clen) :: str,fmt

    type(raystat_type) :: raystats(raystatbuffsize)
    integer(i8b) :: raystatcnt
    type(gadget_constants_type) :: gconst

    !  local counters 
    !-----------------
    
    integer(i8b) :: snapn !< snapshot counter
    integer(i8b) :: rayn  !< ray counter
    integer(i8b) :: raym  !< ray counter, inner loop
    integer(i8b) :: srcn  !< source counter
    integer(i8b) :: rayoops !< count of lasthit > ray%itime
    
    real(r8b) :: rn       !< random number
    real(r8b) :: MB       !< MBs for memory consumption tracking
    
    ! work variables
    !-----------------
    real(r8b) :: outmark
    
#ifdef incHrec
    integer(i8b) :: rindx
    integer(i8b) :: pindx
    integer(i8b) :: lasthitcheck
#endif
    
    
    raystatcnt = 0
    
    ! loop over the snapshots 
    !=========================
    snaps: do snapn = GV%StartSnapNum, GV%EndSnapNum

       !  read in particle and source snapshot
       !----------------------------------------------------------------      
       call readin_snapshot()
       psys%src%lastemit = GV%itime
              
       !  build oct tree.  only need to do this once per snap (for now)
       !----------------------------------------------------------------
       call buildtree(psys,tree,MB,GV%PartPerCell)
       GV%MB = GV%MB + MB
       call setparticleorder(psys, tree)             
       
      
       if (GV%raystats) then
          write(GV%raystatlun) PLAN%snap(snapn)%SrcRays, raystatbuffsize 
       endif
       

       if(GV%JustInit) then
          write(str,"(A,F10.2)") "total memory allocated [MB] = ", GV%MB
          call mywrite(str,verb)
          call mywrite("just initializing", verb)
          call mywrite("",verb)
          stop
       end if

       
       if (GV%DoInitialOutput) then
          GV%OutputIndx = 0
          call output_total_snap(psys)      
          GV%OutputIndx = 1
       end if
       
       
       ! begin ray tracing 
       !------------------------- 
       allocate(active_rays(GV%IonFracOutRays))
       GV%rayoops = 0
       GV%totalhits = 0

       src_rays: do rayn = one, PLAN%snap(snapn)%SrcRays, GV%IonFracOutRays

          
          do raym = 1, GV%IonFracOutRays
            ! begin creation of a ray
            GV%TotalSourceRaysCast = GV%TotalSourceRaysCast + 1                
            GV%itime = GV%itime + 1
          
            !  select a source randomly (weighted by their luminosity)
            rn = genrand_real1() * psys%src(size(psys%src))%Lcdf
            srcn=1
            do while(psys%src(srcn)%Lcdf.LT.rn)
               srcn=srcn+1
               if(srcn.GT.size(psys%src)) then
                  write(*,*) srcn, rn, psys%src(size(psys%src))%Lcdf, size(psys%src)
                  stop "src num > number of sources in mainloop.f90"
               endif
            enddo
                    
            !  create a source ray and calc the impacts
            call src_ray_make( active_rays(raym), psys%src(srcn), GV%itime, GV%dt_s, GV%Lunit, psys%box )

            ! begin stat
            if (GV%raystats) then
             
               raystatcnt = raystatcnt + 1
             
               raystats(raystatcnt)%srcn  = srcn
               raystats(raystatcnt)%start = active_rays(raym)%start  
               raystats(raystatcnt)%ryd   = active_rays(raym)%freq
             
               if (raystatcnt == raystatbuffsize) then
                  write(GV%raystatlun) raystats
                  flush(GV%raystatlun)
                  raystatcnt = 0
               end if
                          
            end if
          ! done stat 
          
            GV%TotalPhotonsCast = GV%TotalPhotonsCast + active_rays(raym)%pini
          ! done creation of a ray
          enddo

         !$OMP PARALLEL FIRSTPRIVATE(raym, raylist, rayoops, TID)
         TID = OMP_GET_THREAD_NUM()
         PRINT *, 'Hello from thread', TID
         !$OMP DO SCHEDULE(DYNAMIC, 1)
          do raym = 1, GV%IonFracOutRays
            ! begin ray tracing and updating 
            call prepare_raysearch(psys, raylist, active_rays(raym))
            call trace_ray(raylist, psys, tree) 
            call update_raylist(raylist,psys%par,psys%box)
            !$OMP ATOMIC
            GV%rayoops = GV%rayoops + raylist%rayoops
            !$OMP ATOMIC
            GV%totalhits = GV%totalhits + raylist%lastnnb
            ! done ray tracing and updating
            ! free up the memory from the globalraylist.
            call kill_raylist(raylist)
          enddo
         !$OMP END DO
         !$OMP END PARALLEL
          ! update some really unused global variables only before output
          ! yfeng1
          GV%IonizingPhotonsPerSec = GV%TotalPhotonsCast / (GV%itime * GV%dt_s)

          !        output routines
          !------------------------
          ! a patch of IonFracOutRays has been processed, write output
          call ion_frac_out(psys, tree )

          ! check if this time step requires a full output
          if ( GV%OutputIndx <= GV%NumTotOuts ) then
             
             ! set correct time marker unit
             if ( trim(GV%OutputTiming) == "standard" ) then
                
                outmark = GV%start_time_code + GV%itime * GV%dt_code
                
             else if ( trim(GV%OutputTiming) == "forced" ) then
                
                if (trim(GV%ForcedUnits) == "mwionfrac") then
                   outmark = GV%mwionfrac
                else 
                   outmark = GV%itime * GV%dt_code
                end if
                
             else
                
                write(*,*) "output type ", trim(GV%OutputTiming), "not recognized"
                stop 
                
             end if


             
             ! check outmark against tabulated output "times"
             if ( outmark >= PLAN%OutputTimes(GV%OutputIndx) ) then
                call output_total_snap(psys)
                GV%OutputIndx = GV%OutputIndx + 1 
             end if
             
          end if
          
          
          ! if we are on the last ray and we havent gotten to the last 
          ! output, do a full output.
          if ( snapn == GV%EndSnapNum ) then
             if ( GV%OutputIndx <= GV%NumTotOuts ) then
                if ( GV%TotalSourceRaysCast==PLAN%snap(snapn)%SrcRays ) then
                   write(*,*) "doing an output on the last ray"
                   call output_total_snap(psys)
                end if
             end if
          end if
          
          
          
          
       end do src_rays
       deallocate(active_rays)


       
       
    end do snaps
    
    close(GV%ionlun)
    
    
    if (GV%raystats) then
       close(GV%raystatlun)
    end if
    
    
  end subroutine mainloop
  
end module mainloop_mod
