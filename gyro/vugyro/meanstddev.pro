FUNCTION NEAR,ARRAY,VALUE
;       WRITTEN BY:     R. E. BELL, 1998
;       MODIFIED:       R. E. BELL, OCT. 2000
;-
  diff=abs(array-value)
;  k=(where(diff eq min(diff)))[0]
;  return,k
  min_val=min(diff)
  k=where(diff eq min_val)
; DRM October 2005: Problems on Mac with spurious returns; check here:
  if (abs(array[k[0]]-value) ne min_val) then begin
     print," Problem inside function NEAR:"
     print," array[k[0]], value, min_val:",array[k[0]], value, min_val
  endif
  return,k[0]
  end

FUNCTION meanstddev, Filename,FuncLabel, values,times,tstart, tend, $
    PlotStr=PlotStr,quiet=quiet
;
; Purpose:
;
;	This function returns structure containing the time average
;	and its standard deviation.
;  Input: 
;      Filename (string) Same as GYRO sim_dir name
;      Values   (array) The flux time series to be analyzed
;      Times    (array) The time vector associated with "Values"
;      Tstart	Analysis window starting time
;      Tend	Analysis window ending time
;  Not used in vugyro version:
;      PlotStr  Controls which plots are drawn
;      quiet    switch to shut off print statements
; The RETURNed result contains: 
;  Avg	Time average of the array Values
;       Estimates of the standard deviation of Avg are based on the
;       lag=1 autocorrelation of a series of sub-interval averages:
;  Trend (Avg over last half) - (Avg over first half)
;  Hump  (Avg over middle half) - (Avg over first and last quarters)
;  Sig1  A_corr falls to 1/sqrt(2*N_interval); the noise level
;  Sig2  A_corr falls to zero  (For general automatic use)
;  Sig3  A_corr falls to -1/sqrt(2*N_interval); the 'far-side' noise level
; Tavg	  The duration of the sub-intervals used for:
; Sig_avg The estimates of the standard deviation (using Tavg sub-intervals)
; Acorr_avg The lag=1 autocorrelation (using sub-interval durations in Tavg)
; Nsig2   The total number of sub-intervals used for the Sig2 estimate.
; Tsig1   The sub-interval duration used for the Sig1 estimate.
; Tsig2   The sub-interval duration used for the Sig2 estimate.
; Tlagu	  The uniform time base used for a 'standard' autocorrelation
; Acorru  The lag=1 autocorrelation based on uniform time base.
;
; > 1.0 get the desired window of data
n_times=n_elements(times)
if (tstart le 0.) then tstart=times[0]
if (tend le 0.) then tend=times[n_times-1]

ind_range=where( times gt tstart and times le tend, count )
if (count gt 1) then begin
  gridvals=times[ind_range]
  funcVals=values[ind_range]
  npoints=n_elements(gridvals)
endif else begin
  print,' GKV_Lag: no times between tstart,tend:',tstart,tend
  output=0.
  return, output
endelse

; > 1.1 set up time intervals; correct large intervals
;help,gridvals
; back_diff is the backward time difference (except for the first element, for which that is undefined).
back_diff=gridVals-shift(gridVals,1)
back_diff[0]=back_diff[1] ; a plausible guess.

;  Show largest back_diff as a check:
ind_wsort=sort(back_diff)
if (not keyword_set(quiet)) then print,' Largest back_diff:',back_diff[ind_wsort[npoints-6:npoints-1]]
; Check for all back_diff == 1.0; probably had a problem reading time variable
if (ind_wsort[0] eq 1.) and  (ind_wsort[npoints-1] eq 1.) then begin
   print,"." ; a blank line
   print,"  ********* Warning: all time differences are 1.00  ************"
   print,"   there may have been a problem reading the time variable "
endif

; Remove large intervals that arise from missing data
ind_missing=where (back_diff gt 10, count)
if(count gt 0) then begin
if (not keyword_set(quiet)) then begin
  for imiss=0,n_elements(ind_missing)-1 do begin
    print,' Data missing between t= ',$
            gridvals[ind_missing[imiss]-1:ind_missing[imiss]]
  endfor
endif
  back_diff[ind_missing]=0. ; missing sections don't count now.
  ; Next, look for realistic time steps to use instead of back_diff=0.
  for  ig=0,n_elements(ind_missing)-1 do begin
    for ishft=1,ind_missing[ig] do begin
      ; look earlier and earlier for a finite weight.
      iuse=ind_missing[ig]-ishft
      if back_diff[iuse] gt 0. then begin
        back_diff[ind_missing[ig]] =back_diff[iuse]  ; gap has a reasonable time step now.
        break
      endif
    endfor
  endfor
endif
; Set a time grid that closes up the jumps:
newgrid=total(back_diff,/cumulative)
newgrid=newgrid+gridvals[0] ; sometimes use this instead of gridvals
 if (not keyword_set(quiet)) then print,' newgrid start and end:'
 if (not keyword_set(quiet)) then print,newgrid[0:3]
n_ng=n_elements(newgrid)
 if (not keyword_set(quiet)) then print,newgrid[n_ng-4:n_ng-1]

; New weighting 11 Oct 2005:
; > 1.15 Set the weights for time averaging
weights=(back_diff+shift(back_diff,-1))/2. ; half of back- and forward-difference
weights[0]=back_diff[1]/2. ; half of the initial forward difference
weights[npoints-1]=back_diff[npoints-1]/2. ; half of the last back difference

; > 1.2 Get average flux
Total_Fluence=total(weights*funcVals)
Avg_Flux=Total_Fluence/total(weights)
 if (not keyword_set(quiet)) then print,' condensed time range is ',newgrid[0],newgrid[npoints-1]
; if (not keyword_set(quiet)) then print,' Time average flux = ',Avg_Flux


; Get the variance of N_corr equal intervals:
;  The minimum value of Ncorr_arr should be >> 1).
;    see notes below about A_CORRELATE.
Ncorr_arr=[7,8,9,10,11,12,14,17,20,24,28,34,40,48,56,68,80,96,112,136,160,192,224,320,448,640,896]
Max_Narr=n_elements(Ncorr_arr)
Max_Ncorr=Max_Narr-1 ; the last index
if ((newgrid[npoints-1]-newgrid[0])/Ncorr_arr[Max_Narr-1] lt 3.) then begin
  ; Use only part of the array
  Max_Ncorr=near(Ncorr_arr,(newgrid[npoints-1]-newgrid[0])/3.)
endif
T_intrvl=(newgrid[npoints-1]-newgrid[0])/Ncorr_arr
dt_min=max(weights) > 3.
ind_long_intrvl=where(t_intrvl gt dt_min)
Max_Ncorr=n_elements(ind_long_intrvl)-1

Nuse_arr=reverse(Ncorr_arr[0:Max_Ncorr])
Max_Ncorr=n_elements(Nuse_arr) ; now represents total number, not the index

Tcorr_arr=fltarr(Max_Ncorr)
Var_corr=fltarr(Max_Ncorr)
sig_corr=fltarr(Max_Ncorr)
A_corr=fltarr(Max_Ncorr)
Acorr_flag1=0
Acorr_flag2=0
Acorr_flag3=0
isig1=0
isig2=0
isig3=0
Sig_max1=0.
Sig_max2=0.
Sig_max3=0.

for itcorr=0,Max_Ncorr-1 do begin

  T_corr=(newgrid[npoints-1]-newgrid[0])/Nuse_arr[itcorr]
  Tcorr_arr[itcorr]=T_corr
  N_corr=Nuse_arr[itcorr]
  Tcen_corr=fltarr(N_corr)
  Avg_corr=fltarr(N_corr)
  Wgt_corr=fltarr(N_corr)
  nd_corr=lonarr(N_corr+1) ; intarr is limited to 32k
  nd_corr[0]=0
  nd_corr[N_corr]=npoints  ; the -1 is included in the "ip" loop below:
  target=T_corr*findgen(N_corr)+newgrid[0]
  for intrvl=1,N_corr-1 do begin
    nd_corr[intrvl]=near(newgrid,target[intrvl])
; The previous line occasionally has the wrong result, so check for an error
;   and correct it.
    if (nd_corr[intrvl] le nd_corr[intrvl-1]) then begin
       print, " After first try:"
       print,'  GKV_Lag problem with nd_corr[',intrvl,']'
       print,intrvl,target[intrvl],newgrid[nd_corr[intrvl-1:intrvl]]
; Try resetting the previous value:
       nd_corr[intrvl-1]=near(newgrid,target[intrvl-1])
       if (nd_corr[intrvl] le nd_corr[intrvl-1]) then begin
         print, " After second try:"
         print,'  GKV_Lag problem with nd_corr[',intrvl,']'
         print,intrvl,target[intrvl],newgrid[nd_corr[intrvl-1:intrvl]]
       endif
    endif
  endfor

  for ip=0,N_corr-1 do begin
    Wgt_corr[ip]=total(weights[nd_corr[ip]:nd_corr[ip+1]-1])
    Avg_corr[ip]=total(weights[nd_corr[ip]:nd_corr[ip+1]-1]*$
          funcVals[nd_corr[ip]:nd_corr[ip+1]-1])/Wgt_corr[ip]
    Tcen_corr[ip]=0.5*(newgrid[nd_corr[ip+1]-1]+newgrid[nd_corr[ip]])
  endfor; End of ip loop

  Avg_of_avgs=total(Wgt_corr*Avg_corr)/total(Wgt_corr)
  if (abs(Avg_of_avgs-Avg_Flux) ge 0.001*abs(Avg_Flux)) then begin
    print,'  ********* Averge flux calculations disagree ****************'
    print, 'Avg_of_avgs,Avg_Flux:',Avg_of_avgs,Avg_Flux
  endif

  Var_corr[itcorr]=sqrt(total(Wgt_corr*(Avg_corr-Avg_Flux)^2)/total(Wgt_corr))
  sig_corr[itcorr]=Var_corr[itcorr]/sqrt(N_corr-1)

; The lagged auto-correlation calculation is not reliable if
;  the number of data points (length of Avg_corr) is not large compared
;  to the maximum of the absolute magnitudes of the elements of lag_array.
;  So N_corr (and the minimum value of Ncorr_arr) should be >> 1).
; Use a better definition for small sample sizes:
; Don't use Wgt_corr in A_corr: which element would be relevant?:
  denom=total((Avg_corr-Avg_Flux)^2)/float(N_corr)
  sh_avg=shift(Avg_corr,-1)
  A_corr[itcorr]=total((Avg_corr[0:N_corr-2]-Avg_Flux)*$
    (sh_avg[0:N_corr-2]-Avg_Flux))/(float(N_corr-1)*denom)

  ; Save maxima for RETURN
  ; test Acorr_flags before setting them so the maximum follows A_corr thresh.
  if (Acorr_flag1 ge 1) then begin
    if (Sig_max1 le 0.) and (sig_corr[itcorr] le sig_corr[itcorr-1]) then begin
      Sig_max1=sig_corr[itcorr-1]
      isig1=itcorr-1
    endif
    if (Acorr_flag2 ge 1) and (Sig_max2 le 0.) then begin
      if (sig_corr[itcorr] le sig_corr[itcorr-1]) then begin
        Sig_max2=sig_corr[itcorr-1]
        isig2=itcorr-1
      endif
    endif
    if (Acorr_flag3 ge 1) and (Sig_max3 le 0.) then begin
      if (sig_corr[itcorr] le sig_corr[itcorr-1]) then begin
        Sig_max3=sig_corr[itcorr-1]
        isig3=itcorr-1
      endif
    endif
  endif
  Acorr_test=1./sqrt(2.*N_corr) ; Avg. abs(A_corr) for uncorrelated sample.
  if (Acorr_flag1 le 0) and (A_corr[itcorr] le Acorr_test) then Acorr_flag1=itcorr
  if (Acorr_flag2 le 0) and (A_corr[itcorr] le 0.) then Acorr_flag2=itcorr
  if (Acorr_flag3 le 0) and (A_corr[itcorr] le -Acorr_test) then Acorr_flag3=itcorr
endfor ; End of itcorr loop

After_loop:

;print,' Sig_max1,Sig_max2=',Sig_max1,Sig_max2
;print,' isig1,isig2=',isig1,isig2
if (Acorr_flag1 le 0) then begin
  if (not keyword_set(quiet)) then $
    print,' GKV_Lag problem: Auto_correlation coeff. always .GE. 1./sqrt(2*N_corr)'
  Sig_max1=max(sig_corr)
  isig1=near(sig_corr,Sig_max1)
  Acorr_flag1=Max_Ncorr-1 ; To be used in plots below
endif else begin
  if (Acorr_flag1 eq Max_Ncorr-1) then begin
    Sig_max1=max(sig_corr)
    isig1=near(sig_corr,Sig_max1)
  endif
  if (Sig_max1 le 0.) then begin
    Sig_max1=sig_corr[Acorr_flag1]
    isig1=Acorr_flag1
  endif
endelse
if (Acorr_flag2 le 0) then begin
  if (not keyword_set(quiet)) then $
     print,' GKV_Lag: Auto_correlation coeff. always positive'
  Sig_max2=max(sig_corr)
  isig2=near(sig_corr,Sig_max2)
  Acorr_flag2=Max_Ncorr-1 ; To be used in plots below
endif else begin
  if (Acorr_flag2 eq Max_Ncorr-1) then begin
    Sig_max2=max(sig_corr)
    isig2=near(sig_corr,Sig_max2)
  endif
  if (Sig_max2 le 0.) then begin
    Sig_max2=sig_corr[Acorr_flag2]
    isig2=Acorr_flag2
  endif
endelse
if (Acorr_flag3 le 0) then begin
  if (not keyword_set(quiet)) then $
    print,' GKV_Lag: Auto_correlation coeff. never very negative'
  Sig_max3=max(sig_corr)
  isig3=near(sig_corr,Sig_max3)
  Acorr_flag3=Max_Ncorr-1 ; To be used in plots below
endif else begin
  if (Acorr_flag3 eq Max_Ncorr-1) then begin
    Sig_max3=max(sig_corr)
    isig3=near(sig_corr,Sig_max3)
  endif
  if (Sig_max3 le 0.) then begin
    Sig_max3=sig_corr[Acorr_flag3]
    isig3=Acorr_flag3
  endif
endelse
;print,' Sig_max1,Sig_max2=',Sig_max1,Sig_max2
;print,' isig1,isig2=',isig1,isig2


tavg_variance=sqrt(total(weights*(FuncVals-Avg_Flux)^2)/total(weights))

imid=near(newgrid,(newgrid[npoints-1]+newgrid[0])/2.)
Avg_1=total(weights[0:imid-1]*funcVals[0:imid-1])/$
          total(weights[0:imid-1])
Avg_2=total(weights[imid:npoints-1]*funcVals[imid:npoints-1])/$
          total(weights[imid:npoints-1])
diff_12=Avg_2-Avg_1
i1qrt=near(newgrid,newgrid[0]+(newgrid[npoints-1]-newgrid[0])/4.)
i3qrt=near(newgrid,newgrid[0]+3.*(newgrid[npoints-1]-newgrid[0])/4.)
Avg_ends=(total(weights[0:i1qrt-1]*funcVals[0:i1qrt-1])$
          + total(weights[i3qrt:npoints-1]*funcVals[i3qrt:npoints-1]) ) $
          /(total(weights[0:i1qrt-1])+total(weights[i3qrt:npoints-1]) )
Avg_center=total(weights[i1qrt:i3qrt-1]*funcVals[i1qrt:i3qrt-1])/$
          total(weights[i1qrt:i3qrt-1])
diff_ctr_ends= Avg_center-Avg_ends
if (not keyword_set(quiet)) then print,' Time average flux = ',Avg_Flux
if (not keyword_set(quiet)) then begin
   maxval=abs(diff_12) > Sig_max1 > Sig_max2 > Sig_max3
   mag=fix(alog(maxval)/alog(10.)) > 0
   ndec=7-mag ; number of decimal places in 10 character field
   fmtstr=string("(4F10.",ndec,")",format='(A,I1,A)')
   print, "  diff_12 , Sig_max1, Sig_max2, Sig_max3:"
   print,diff_12,Sig_max1,Sig_max2,Sig_max3,format=fmtstr
endif

;  Interpolate on uniform grid to calculate A_corr for unaveraged f(t)
;   different spacings look similar : delt_arr=[0.5,1.,2.,5.,10.,20.]
delt=1.
Nint_times= fix((newgrid[npoints-1]-newgrid[0])/delt)
Tuni=delt*indgen(Nint_times) + newgrid[0]
Funi=interpol(FuncVals,newgrid,Tuni)
;;;; Tmax_uni=Tcorr_arr[isig2] < (newgrid[npoints-1]-newgrid[0])/4.
Tmax_uni=Tcorr_arr[Max_Ncorr-1]
Max_lag=fix(Tmax_uni/delt)
Lagu_arr=indgen(Max_lag) ; let first element be zero to get A_corr(0)=1.00
Tlag=delt*indgen(Max_lag)
Auni_corr=A_correlate(Funi,Lagu_arr)

Returning:
output={avg:Avg_Flux, Sig1:Sig_max1, Sig2:Sig_max2, Sig3:Sig_max3, $
  Tlagu:Tlag, Acorru:Auni_corr, Trend:diff_12, Hump:diff_ctr_ends, $
  Tavg:Tcorr_arr, Acorr_avg:A_corr, Sig_avg:Sig_corr,$
  Nsig2:Nuse_arr[isig2], Tsig1:Tcorr_arr[isig1], Tsig2:Tcorr_arr[isig2] }
RETURN, output

END ; ****** GKVs1D::MeanStdDev ****** ;
