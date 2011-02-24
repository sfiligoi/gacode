pro user_control, SKIP_LARGE = do_skip_large, REMOTE_DIR=use_remote_dir

   common GLOBAL

  ;;----------------------------------------------------
  ;; REMOTE SIMULATION DIRECTORY (0=off,1=on)
  ;;
  DEFAULT, use_remote_dir, 1
  remotedir_flag = use_remote_dir
  remotedir = '/global/scratch/sd/bassem/VZ_runs_142111/newruns'
  ;;----------------------------------------------------

  ;;----------------------------------------------------
  ;; READ LARGEST FILES (0=read,1=skip,2=only read u.out)
  ;; Do it this way b/c skip_large is a global variable (CH)
  DEFAULT, do_skip_large, 0
  skip_large = do_skip_large
  ;;----------------------------------------------------
  
  ;;----------------------------------------------------
  ;; QUIET MODE (affects meanstddev)
  ;;
  quiet = 1
  ;;----------------------------------------------------

  ;;----------------------------------------------------
  ;; WINDOW DIMENSIONS
  ;;
  sx = 720
  sy = 480  
  ;;----------------------------------------------------

  ;;----------------------------------------------------
  ;; POSTSCRIPT
  ;;
  ;; Default postscript plot size
  ;;
  xsize = 16
  ysize = 10
  
  ;; 0=ps, 1=eps
  ps_val = 0

  ;; 0=b/w, 1=color
  ps_color = 1

  ;; 0=label on, 1=label off
  ps_label = 1

  ;; 0=thin, 1=thick
  ps_thick = 1
  ;;----------------------------------------------------

end

