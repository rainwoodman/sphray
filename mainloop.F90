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
  use ion_temperature_update, only: update_intersection
  use mt19937_mod, only: genrand_real1
  use output_mod, only: output_total_snap, ion_frac_out
  
  ! variables
  use global_mod, only: psys
  use global_mod, only: tree
  use global_mod, only: GV
  use accounting_mod, only: AV
  use accounting_mod
  use global_mod, only: PLAN
  use global_mod, only: active_rays
  use config_mod, only: CV
  use resolve_mod, only: resolution_type, resolve_more, resolve_all, resolution_get_resolved_intersection, prepare_resolution, kill_resolution
  
  implicit none
  
  integer(i8b), parameter :: one = 1

  integer(i4b), parameter :: rays_per_leaf = 2
  integer(i4b), parameter :: rays_per_dt = 1
contains
  
  !> this is the main driver of SPHRAY
  !======================================
  subroutine mainloop()
    implicit none
    integer TID, OMP_GET_THREAD_NUM, NTRD, OMP_GET_MAX_THREADS
    real(r8b) OMP_GET_WTIME
    type(raylist_type),allocatable :: raylists(:)       !< ray/particle intersections
    type(accounting_variables_type),allocatable :: localAVs(:)       !< ray/particle intersections
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
    integer(i8b) :: raybatch !< ray counter, at integer batch startings
    integer(i4b) :: rayn  !< ray counter, inner loop
    integer(i8b) :: srcn  !< source counter
    integer(i8b) :: j !< resolution encoded_impacts counter
    integer(i8b) :: good_nnb_sum, iterations, small_count
    
    type(resolution_type):: resolution !< thread racing and causality resolution
    type(intersection_type):: intersection!< intersection
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
    snaps: do snapn = CV%StartSnapNum, CV%EndSnapNum

       !  read in particle and source snapshot
       !----------------------------------------------------------------      
       call readin_snapshot()
       psys%src%lastemit = GV%itime
              
       !  build oct tree.  only need to do this once per snap (for now)
       !----------------------------------------------------------------
       call buildtree(psys,tree,MB,CV%PartPerCell)
       GV%MB = GV%MB + MB
       call setparticleorder(psys, tree)             
       
      
       if (CV%raystats) then
          write(GV%raystatlun) PLAN%snap(snapn)%SrcRays, raystatbuffsize 
       endif
       

       if(CV%JustInit) then
          write(str,"(A,F10.2)") "total memory allocated [MB] = ", GV%MB
          call mywrite(str,verb)
          call mywrite("just initializing", verb)
          call mywrite("",verb)
          stop
       end if

       
       if (CV%DoInitialOutput) then
          GV%OutputIndx = 0
          call output_total_snap(psys)      
          GV%OutputIndx = 1
       end if
       
       
       ! begin ray tracing 
       !------------------------- 
       allocate(active_rays(CV%IonFracOutRays))
       call ion_frac_out(psys, tree )
       call timer_init(AV, 1, "tracing", running = .False.)
       call timer_init(AV, 2, 'resolve', running=.False.)
       call timer_init(AV, 3, 'good', running=.False.)
       call timer_init(AV, 4, 'bad', running=.False.)
       call timer_init(AV, 5, 'update', running=.False.)
       call timer_init(AV, 6, 'preRes', running=.False.)

       src_rays: do raybatch = one, PLAN%snap(snapn)%SrcRays, CV%IonFracOutRays

          
          do rayn = 1, CV%IonFracOutRays
            ! begin creation of a ray
            AV%TotalSourceRaysCast = AV%TotalSourceRaysCast + 1                
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
            call src_ray_make( active_rays(rayn), psys%src(srcn), GV%itime, GV%dt_s, GV%Lunit, psys%box )

            ! begin stat
            if (CV%raystats) then
             
               raystatcnt = raystatcnt + 1
             
               raystats(raystatcnt)%srcn  = srcn
               raystats(raystatcnt)%start = active_rays(rayn)%start  
               raystats(raystatcnt)%ryd   = active_rays(rayn)%freq
             
               if (raystatcnt == raystatbuffsize) then
                  write(GV%raystatlun) raystats
                  flush(GV%raystatlun)
                  raystatcnt = 0
               end if
                          
            end if
          ! done stat 
          
            AV%TotalPhotonsCast = AV%TotalPhotonsCast + active_rays(rayn)%pini
          ! done creation of a ray
          enddo

         allocate(raylists(CV%IonFracOutRays))
         NTRD = OMP_GET_MAX_THREADS()

         allocate(localAVs(0:NTRD-1))
         do tid = 0, NTRD - 1
           call clear_accounting_variables(localAVs(tid))
         enddo

         !$OMP PARALLEL FIRSTPRIVATE(rayn, TID, j, intersection)

         TID = OMP_GET_THREAD_NUM()
         !PRINT *, 'Hello from thread', TID, NTRD
         ! begin ray tracing and updating 
         call timer_resume(AV, 1)
         !$OMP DO SCHEDULE(DYNAMIC, 1)
          do rayn = 1, CV%IonFracOutRays
            call prepare_raysearch(raylists(rayn))
            call trace_ray(rayn, raylists(rayn), psys, tree) 
          enddo
         !$OMP END DO
         !$OMP END PARALLEL
         call timer_pause(AV, 1)

         ! this section resolves the races and causalities
         call timer_resume(AV, 6)
         call prepare_resolution(resolution, raylists, active_rays, psys%par)
         call timer_pause(AV, 6)
           small_count = 0
           good_nnb_sum = 0
           iterations = 0

         do while(resolution%remaining_nnb > 0)
           if (NTRD == 1) then
             call resolve_all(resolution, raylists)
           else
             call timer_resume(AV, 2)
             call resolve_more(resolution, raylists)
             call timer_pause(AV, 2)
           endif
           ! this section updates the intersections
           call timer_resume(AV,5)
           if (resolution%good_nnb < 5) then
              small_count = small_count + 1
           endif
           !$OMP PARALLEL FIRSTPRIVATE(j, TID, intersection) IF(resolution%good_nnb > 1)
           TID = OMP_GET_THREAD_NUM()
           !$OMP DO SCHEDULE(DYNAMIC, 1)
           do j = 1, resolution%good_nnb
             call resolution_get_resolved_intersection(resolution, raylists, j, intersection)
!             print *, resolution%good_nnb, 'inte', intersection%pindx, intersection%rayn, intersection%t
!             print *, active_rays(intersection%rayn)%emit_time, psys%par(intersection%pindx)%lasthit
             call update_intersection(intersection, psys%par,psys%box, localAVs(TID))
           enddo 
           !$OMP END DO
           !$OMP END PARALLEL
           call timer_pause(AV,5)
           good_nnb_sum = good_nnb_sum + resolution%good_nnb
           iterations = iterations + 1
         end do

         ! free up the memory from the globalraylist.
         ! done ray tracing and updating
         call kill_resolution(resolution)
         do tid = 0, NTRD - 1
           call reduce_accounting_variables(localAVs(tid))
         enddo
         deallocate(localAVs)
         do rayn = 1, CV%IonFracOutRays
           call kill_raylist(raylists(rayn))
         enddo
         deallocate(raylists)

         ! if vacuum BCs and exiting box, claim the leftovers
         if(psys%box%tbound(1)==0) then
           do rayn = 1, CV%IonFracOutRays
             if(.not. active_rays(rayn)%exhausted) then
               AV%PhotonsLeavingBox = AV%PhotonsLeavingBox + active_rays(rayn)%pcnt
             end if
           enddo
         end if

          ! update some really unused global variables only before output
          ! yfeng1
          AV%IonizingPhotonsPerSec = AV%TotalPhotonsCast / (GV%itime * GV%dt_s)

          !        output routines
          !------------------------
          ! a patch of IonFracOutRays has been processed, write output
          call ion_frac_out(psys, tree )
         call timer_print(AV, 1)
         call timer_print(AV, 2)
         call timer_print(AV, 3)
         call timer_print(AV, 4)
         call timer_print(AV, 5)
         call timer_print(AV, 6)
         print *, real(good_nnb_sum) / iterations, small_count, iterations
          ! check if this time step requires a full output
          if ( GV%OutputIndx <= GV%NumTotOuts ) then
             
             ! set correct time marker unit
             if ( trim(CV%OutputTiming) == "standard" ) then
                
                outmark = GV%start_time_code + GV%itime * GV%dt_code
                
             else if ( trim(CV%OutputTiming) == "forced" ) then
                
                if (trim(CV%ForcedUnits) == "mwionfrac") then
                   outmark = GV%mwionfrac
                else 
                   outmark = GV%itime * GV%dt_code
                end if
                
             else
                
                write(*,*) "output type ", trim(CV%OutputTiming), "not recognized"
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
          if ( snapn == CV%EndSnapNum ) then
             if ( GV%OutputIndx <= GV%NumTotOuts ) then
                if ( AV%TotalSourceRaysCast==PLAN%snap(snapn)%SrcRays ) then
                   write(*,*) "doing an output on the last ray"
                   call output_total_snap(psys)
                end if
             end if
          end if
          
          
          
          
       end do src_rays
       deallocate(active_rays)


       
       
    end do snaps
    
    close(GV%ionlun)
    
    
    if (CV%raystats) then
       close(GV%raystatlun)
    end if
    
    
  end subroutine mainloop
  

end module mainloop_mod
