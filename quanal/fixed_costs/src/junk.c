      /*}
      else
	{
	  uint fail_cnt=0;
	  fprintf(logfile,KBLU "One-shot calibration for exit rate of 45%...\n" RESET);

	  double grid[homotopy_times];

	  linspace(0.1,0.2,homotopy_times,grid);
	  int h;
	  for(h=homotopy_times-1;h>=0;h--)
	    {
	      exit_rate_target = grid[h];
	      fprintf(logfile,KNRM "\n\tCalibrating to exporter exit rate of %0.3f" RESET,exit_rate_target);
	      status = find_root_deriv_mkl(&f);
	      if(status && status != 99)
		{
		  fail_cnt +=1;
		  //fprintf(logfile,KRED "\nError solving calibration function!\n" RESET);
		  //break;
		}
		}

	  if(status == 0 || status == 99)
	    {
	      linspace(0.05,0.1,homotopy_times,grid);
	      for(h=homotopy_times-1;h>=0;h--)
		{
		  exit_rate_target = grid[h];
		  fprintf(logfile,KNRM "\n\tCalibrating to exporter exit rate of %0.3f" RESET,exit_rate_target);
		  status = find_root_deriv_mkl(&f);
		  if(status && status != 99)
		    {
		      fail_cnt+=1;
		      //fprintf(logfile,KRED "\nError solving calibration function!\n" RESET);
		      //break;
		    }
		}
	    }

	  if(status == 0 || status == 99)
	    {
	      linspace(0.025,0.05,homotopy_times,grid);
	      for(h=homotopy_times-1;h>=0;h--)
		{
		  exit_rate_target = grid[h];
		  fprintf(logfile,KNRM "\n\tCalibrating to exporter exit rate of %0.3f" RESET,exit_rate_target);
		  status = find_root_deriv_mkl(&f);
		  if(status && status != 99)
		    {
		      fail_cnt+=1;
		      //fprintf(logfile,KRED "\nError solving calibration function!\n" RESET);
		      //break;
		    }
		}
	    }

	  if(status == 0 || status == 99)
	    {
	      linspace(0.02,0.025,homotopy_times,grid);
	      for(h=homotopy_times-1;h>=0;h--)
		{
		  exit_rate_target = grid[h];
		  fprintf(logfile,KNRM "\n\tCalibrating to exporter exit rate of %0.3f" RESET,exit_rate_target);
		  status = find_root_deriv_mkl(&f);
		  if(status && status != 99)
		    {
		      fail_cnt+=1;
		      //fprintf(logfile,KRED "\nError solving calibration function!\n" RESET);
		      //break;
		    }
		}
	    }

	  if(fail_cnt>20)
	    {
	      fprintf(logfile,KRED "\nToo many failures in homotopy process!\n" RESET);
	      status=1;
	      }*/
