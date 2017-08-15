;@fit_functions
pro conf_error2d,x,y,z,zerr,p0,perror0,model,perror,bestfit,pvarunit,bestchisq,zfit,psave,doplot=doplot,pmin0=pmin0,log=log,d2=d2,parinfo=parinfo,weights=weights,wpar=wpar,delchi0=delchi0
  
; z=image, zerr = error image 
; p0 = parameter structure/array 
; 


 if n_params() eq 0 then begin
     print,'syntax - conf_error,x,y,z,zerr,p0,perror0,model,perror'
     return
  endif 
  
  ;;default 90% confidence delchi=2.706
  if n_elements(delchi0) eq 0 then delchi0=2.706
  ;;1-sigma
;  delchi0=1.
  
  npar=n_elements(p0) ;number of model parameters
;  !p.multi=[0,round(npar/2),2]

  perror=dblarr(2,npar) ;upper and lower errors for each parameter
  
  com='zfit0='+model+'(x,y,p0)' 
  
  tmp=execute(com) ;return model image from initial parameters
  
  nvar=20. ;number of iterations (variations)             
;  pmin=0.9*p0
;  pmax=1.1*p0

  w0=where(perror0 eq 0.,nw0) ;if initial parameter error = 0, 
                              ; replace with 0.1*parameter
  if nw0 gt 0 then perror0[w0]=0.1*p0[w0]
  w0=where(p0 eq 0.,nw0) ; if initial parameter = 0, 
  if nw0 gt 0 then perror0[w0]=1.  ; set parameter error = 1

  mult=[replicate(2.,npar)] ;array with one element per parameter (all set to 2)
  pmin=p0-perror0*mult ;create array to hold parameter-2*error values
  if n_elements(pmin0) gt 0 then begin ;check if pmin0 argument passed,
     wmin=where(pmin-pmin0 lt 0.,nwmin)
     if nwmin gt 0 then pmin[wmin]=pmin0[wmin] ;replaces all pmin values that
                                ; are out of allowed range with
                                ; minimum allowed
  endif 
  pmax=p0+perror0*mult ;create array to hold parameter+error values
  if n_elements(doplot) eq 0 then doplot=1
  delchi=dblarr(nvar) & chisq=dblarr(nvar) ;create arrays to hold delchi

                              ; and chi2 from each iteration
  bestfit=dblarr(npar) & pvarunit=bestfit & bestchisq=bestfit ;array to hold
                                ;best fit model parameters (bestchisq
                                ;and pvarunit not used)
  psave=dblarr(n_elements(p0),nvar) ;save best fit model params from each iteration (overwritten each parameter loop like delchi)
              
  wno0=where(z ne 0 and zfit0 ne 0 and zerr ne 0,nwno0) ;flag zero pixels
  ;save initial chi2
  chisq0=total((z[wno0]-zfit0[wno0])^2/zerr[wno0]^2)*n_elements(zerr)/nwno0  
  bestfit=p0
  doag=0
  multlow=mult ;initialize arrays (one element per parameter)
  multup=mult     
  if n_elements(log) eq 0 then log=intarr(npar)
  ppix=intarr(35) ;why 35?
  ppix[[2,3,4,5,9,13,17,21]]=1 ;numbers correspond to x0,y0,sigr,and sig0-3
                               ;  indices in parameter array
  once=0
  go=1
  ; wpar=array of indices corresponding to certain parameters? new_sn1987a_fit
  ;   passes [4,5], which would correspond to sigr and r0 - 
  ;   only calculates errors for these parameters?
  if n_elements(wpar) eq 0 then wpar=indgen(npar) else npar=n_elements(wpar)
  print, 'n wpar = ',npar
  for i=0,npar-1 do begin
     if n_elements(parinfo) gt 0 then begin ;make sure parinfo is defined
        if parinfo[wpar[i]].fixed eq 1 then go=0 else go=1 ;skip fixed params
        if go eq 0 then print, 'skipping parameter ',i
     endif 
     if go then begin 
        print, 'beginning loop over parameter '+string(i)
        p=p0 
;     pvarunit[i]=0.01*(pmax[i]-pmin[i])
    ;'    pvar=dblarr(nvar)
        pvar=dblarr(floor(nvar/2d)*2d) ;-> EAH ;array with one element per iteration
        pvar[nvar/2.:*]=(findgen(nvar/2.)+1.5)/(nvar/2.+1)*(pmax[wpar[i]]-p0[wpar[i]])+p0[wpar[i]] ;fill last half with (multiples of uppererrorbar) + paramvalue -- steps through parameter values within error range
        pvar[0:nvar/2.-1]=(findgen(nvar/2.)+1.5)/(nvar/2.+1)*(p0[wpar[i]]-pmin[wpar[i]])+pmin[wpar[i]] ;fill first half with (multiples of uppererrorbar) + paramvalue
        olddelchi=0
        doagain:
        nvar=n_elements(pvar)
;        print,'p0 = ',p0
;        print,'pmin = ',pmin
;        print,'pmax = ',pmax
        for j=0,nvar-1 do begin  ;iteration loop
           if delchi[j] eq 0. then begin ;if iteration not already done
              ;; but delchi is not reset after parameter loop! so
              ;; after first parameter, this will always fail...fixed
              p[wpar[i]]=pvar[j] ; shift parameter
;           if npar eq 1 then begin 
;              if not d2 then com='yfit='+model+'(x,p)' else com='yfit='+model+'(x,x,p)'
;              tmp=execute(com)
;           endif else begin 
              ;;need to fit free params
              if n_elements(parinfo) gt 0 then nparinfo=parinfo else $
                 nparinfo = parinfo_struct(npar) ;set up parameter array
              nparinfo[wpar[i]].fixed=1          ;  for this iteration
              nparinfo[wpar[i]].value=pvar[j]    ;fix current param to shifted 
                                                 ; value, fit others
              newp=mpfit2dfun(model,x,y,z,zerr,p,parinfo=nparinfo,yfit=zfit,/quiet, ftol=0.0001)
;           endelse 

              wno0=where(z ne 0 and zfit ne 0 and zerr ne 0,nwno0)
              psave[*,j]=newp ;save fit parameters from this iteration
              chisq[j]=total((z[wno0]-zfit[wno0])^2/zerr[wno0]^2)*n_elements(zerr)/nwno0
              print, 'iteration, value: ', j,pvar[j]
              print, 'new, original chi2: ',chisq[j], chisq0
              delchi[j]=chisq[j]-chisq0
           endif 
;           print,j

        endfor ;end loop over iterations for parameter wpar[i]

        ; compare delchi values from this iteration (from all 
        ;  parameters in grid) to saved delchis
        wlow=where(pvar le p0[wpar[i]],nlow) ;iteratins with paramvalue < initial
        wupp=where(pvar ge p0[wpar[i]],nupp) ;iterations with paramvalue > initial
        if max(delchi[wlow]) lt delchi0 then begin ;check low param value delchi
           multlow[wpar[i]]=multlow[wpar[i]]*2.    ; if all pass, set multlow=4
           pmin[wpar[i]]=p0[wpar[i]]-perror0[wpar[i]]*multlow[wpar[i]] ;pmin[wpar[i]] ; set param low error bound to param-initialerror*4 (expand range)
           if n_elements(pmin0) gt 0 then begin ;fails unless pmin0 was passed
              if pmin[wpar[i]] le pmin0[wpar[i]] then begin 
                 pmin[wpar[i]]=pmin0[wpar[i]]   
                 if not once then begin 
                    once=1 
                    doag=1
                 endif else begin
                    doag=0
                    once=0
                 endelse 
              endif else doag=1 ;do again if new pmin > previous pmin
           endif else doag=1 ;if pmin0 not set and lower bound expanded, then do again
        endif ;endif need to expand lower parameter bound

        if max(delchi[wupp]) lt delchi0 then begin ;check if upper bound needs to be expanded to reach the desired conf. threshold
;           print,'wupp delchi ',max(delchi[wupp])
           multup[wpar[i]]=multup[wpar[i]]*2. ;if all new high params pass
           pmax[wpar[i]]=p0[wpar[i]]+perror0[wpar[i]]*multup[wpar[i]] ;expand upper bound
;        pvar=findgen(nvar+1)/(nvar)*(pmax[wpar[i]]-pmin[wpar[i]])+pmin[wpar[i]]
;        pvar=[pvar[wlow],indgen(nupp)/(nupp-1.)*(pmax[wpar[i]]-p0[wpar[i]])+p0[wpar[i]]]
;        pvar=indgen(nvar)/(nvar-1)*(pmax[wpar[i]]*2.-pmin[wpar[i]])+pmin[wpar[i]]
;        pvar=[pvar[wlow],pvar[wupp]*2.]
           if max(delchi[wupp]) eq olddelchi then doag=0 else doag=1 ;do again if upper bound has changed (max delchi is different than last time around)
           olddelchi=max(delchi[wupp]);set olddelchi equal to biggest new upper delchi

        endif ;endif expand upper bound
        if doag then begin 
;           pvar=dblarr(nvar)
;           pvar[nvar/2.+0.5:*]=(findgen(nvar/2.)+1.5)/(nvar/2.+1)*(pmax[wpar[i]]-p0[wpar[i]])+p0[wpar[i]]
;           pvar[0:nvar/2.-1]=(findgen(nvar/2.)+1.5)/(nvar/2.+1)*(p0[wpar[i]]-pmin[wpar[i]])+pmin[wpar[i]]
           nolow=0 & nohigh=0
           if parinfo[wpar[i]].limited[0] eq 1 then begin
              if pmin[wpar[i]] lt parinfo[wpar[i]].limits[0] or max(delchi[wlow]) gt delchi0 then nolow=1 ;decide if should do again the lower range of par
           endif 
           if parinfo[wpar[i]].limited[1] eq 1 then begin
              if pmax[wpar[i]] gt parinfo[wpar[i]].limits[1] or max(delchi[wupp]) gt delchi0 then nohigh=1 ;decide if should do again the upper range of par
           endif 
           nnew=nvar/4
           ppmin=min(pvar)-perror0[wpar[i]]*multlow[wpar[i]] ;expand ranges
           ppmax=max(pvar)+perror0[wpar[i]]*multup[wpar[i]]
           if ppmin ne parinfo[wpar[i]].limits[0] and ppmax ne parinfo[wpar[i]].limits[1] then begin  ;check new range against limits
              if ppmin lt parinfo[wpar[i]].limits[0] then ppmin=parinfo[wpar[i]].limits[0]
              if ppmax gt parinfo[wpar[i]].limits[1] and parinfo[wpar[i]].limits[1] ne 0 then ppmax=parinfo[wpar[i]].limits[1]
              plow=(findgen(nvar/4.)+1.5)/(nvar/4.+1)*(min(pvar)-ppmin)+ppmin
              phigh=(findgen(nvar/4.)+1.5)/(nvar/4.+1)*(ppmax-max(pvar))+max(pvar) ;create upper and lower extensions for pvar array

              if parinfo[wpar[i]].limited[0] eq 1 then begin ;if lower limit 
;              winlim=where(pvar ge parinfo[wpar[i]].limits[0] and pvar le parinfo[wpar[i]].limits[1])
;              pvar=pvar[winlim]
                 woutlim=where(pvar gt parinfo[wpar[i]].limits[1] or pvar lt parinfo[wpar[i]].limits[0],nwout) ;number of points in pvar outside the limits
                 nps=n_elements(psave[0,*]) ;original number of iterations
                 if nwout eq 0 then begin ;if no points outside limits
                    doag=0
                    if not nolow then begin ;expand arrays to include new iterations
                       pvar=[plow,pvar] ;reset pvar grid
                       delchi=[dblarr(nnew),delchi] ;expand delchi array
                       chisq=[dblarr(nnew),chisq] ;expand chisq array
                       ;pstmp=dblarr(35,nps+nnew)
                        pstmp=dblarr(n_elements(p0),nps+nnew) ;change made by EAH
                       pstmp[*,nnew:*]=psave
                       psave=pstmp ;expand the psave array
;                       psave=[dblarr(35,nnew),psave]
                    endif 
                    if not nohigh then begin ;same as previous loop
                       nps=n_elements(psave[0,*])
                       pvar=[pvar,phigh]
                       delchi=[delchi,dblarr(nnew)]
                       chisq=[chisq,dblarr(nnew)]
 ;                      pstmp=dblarr(35,nps+nnew)
                       pstmp=dblarr(n_elements(p0),nps+nnew) ;change made by EAH
                       pstmp[*,0:nps-1]=psave
                       psave=pstmp
;                       psave=[psave,dblarr(35,nnew)]
                    endif 
;                    pvar=[plow,pvar,phigh]
;                    delchi=[dblarr(nnew),delchi,dblarr(nnew)]
;                    chisq=[dblarr(nnew),chisq,dblarr(nnew)]
;                    stop
                    if not nohigh or not nolow then goto,doagain ;repeat iteration loop
                 endif 
              endif else begin ;do again if no lower limit
                 doag=0
                 nps=n_elements(psave[0,*]) ;expand arrays for new iterations
                 pvar=[plow,pvar,phigh] ;update pvar grid
                 delchi=[dblarr(nnew),delchi,dblarr(nnew)]
                 chisq=[dblarr(nnew),chisq,dblarr(nnew)]
 ;                pstmp=dblarr(35,nps+nnew*2) ;JR
                 pstmp=dblarr(n_elements(p0),nps+nnew*2) ;EAH -> this appears to work
 ;                stop
                 pstmp[*,nnew:nps+nnew-1]=psave
;                 pstmp[*,nnew:nps-1]=psave
                 psave=pstmp
;                 psave=[dblarr(35,nnew),psave,dblarr(35,nnew)]
                 goto,doagain
              endelse 
           endif ;endif ppmin,ppmax not equal to limits
           
        endif ;endif doag
;     polint,delchi[wlow],pvar[wlow],delchi0,intlow
;     intlow=base_interp(pvar[wlow],delchi[wlow],delchi0,order=2,minsig=100.)

;     polint,delchi[wupp],pvar[wupp],delchi0,intupp
;     intupp=base_interp(pvar[wupp],delchi[wupp],delchi0,order=2,minsig=100.)

;     intlow=interpol(pvar[wlow],delchi[wlow],delchi0)
;     intupp=interpol(pvar[wupp],delchi[wupp],delchi0)

        intlow=findval(delchi[wlow],pvar[wlow],delchi0) ;linearly interpolate to get parameter value that would be associated with desired conf. threshold
        print, 'pvar[wlow],delchi[wlow] = '
        print, pvar[wlow]
        print, delchi[wlow]
       intupp=findval(delchi[wupp],pvar[wupp],delchi0)
        if n_elements(intlow) gt 1 then intlow=min(intlow) ;in case multiple values would give delchi=delchi0
        if n_elements(intupp) gt 1 then intupp=max(intupp)
        ;if findval failed, use interpol instead
        if intlow eq -1 then intlow=interpol(pvar[wlow],delchi[wlow],delchi0)
        if intupp eq -1 then intupp=interpol(pvar[wupp],delchi[wupp],delchi0)
;        print, 'p0[wpar[i]],intlow = ',intlow,p0[wpar[i]]
        lowerr=p0[wpar[i]]-intlow ;get lower 90% error bar     
        upperr=intupp-p0[wpar[i]] ;get upper 90% error bar
;        print, 'max delchi[wlow], delchi0 = ',max(delchi[wlow]),delchi0
        ;set lower error bar = best fit parameter (upper limit) if necessary
        if max(delchi[wlow]) lt delchi0 and lowerr lt 0. then lowerr=p0[wpar[i]]
        if keyword_set(doplot) then begin ;plot stuff
  ;         stop
;        w=where(delchi lt 5)
           plot,pvar,delchi,xlog=log[wpar[i]],yrange=[min(delchi)-2,10]
           oplot,[p0[wpar[i]]-lowerr,p0[wpar[i]]-lowerr],[min(delchi),max(delchi)*2.],line=2
           oplot,[p0[wpar[i]]+upperr,p0[wpar[i]]+upperr],[min(delchi),max(delchi)*2.],line=2
           oplot,[p0[wpar[i]],p0[wpar[i]]],[min(delchi),max(delchi)*2],line=0
           oplot,[min(pvar),max(pvar)],[delchi0,delchi0],line=1
;        oplot,[bestfit[wpar[i]],bestfit[wpar[i]]],[min(delchi),max(delchi)*2],line=1
;        k=get_kbrd(10)   
;        if k eq 's' then stop
        endif ;end plot stuff
;     print,bestfit[wpar[i]],lowerr,upperr,perror0[wpar[i]]

;        if keyword_set(dosys) then begin
        ;;;interpolate other parameters to get effect on systematics
        whpar=0
        ns=n_elements(psave[0,*]) ;(final) number of iterations
        range=fltarr(3,n_elements(psave[*,0]))
        for k=0,n_elements(psave[*,0])-1 do begin ;loop over all p0 params
           ps=psave[k,*] ;parameter k values from all iterations
           s=indgen(n_elements(ps)) ;array of indices for ps
           if ps[0] ne min(ps) then s=reverse(s) ;reverse array (why?)
           ps=ps[s]
           if k eq 5 then print,'ps = ', ps
           mpar=min(delchi[s],m) ;find iteration with min delchi (equivalent to min chisq)

           if k eq 5 then print, 'mpar, m, ps[m] = ',mpar,m,ps[m]
           if m lt ns-1 and m gt 0 then begin ;if not 0th or last element
              lowpar=where(ps[0:m-1] lt ps[m],nlow) ; or ps[0:m-1] gt ps[m],nlow) ;find iterations where parameter k is < value associated with min delchi
              highpar=where(ps[m+1:*] gt ps[m],nhigh) ; or ps[m+1:*] lt ps[m],nhigh)
              if abs(nlow-m) lt 5 and abs(nhigh-m) lt 5 then begin ;choose if the parameter should be included in systematics check? if fewer than 5 on each side are greater than ps[m], then continue
                 whpar=[whpar,k] ;add parameter k to whpar list
                 ;interpolate parameter k that gives delchi0
                 intlow=interpol(ps[0:m],delchi[s[0:m]],delchi0)
                 intupp=interpol(ps[m:*],delchi[s[m:*]],delchi0)
                 ;save interpolated 90% upper and lower
                 ;bounds and the best fit value for
                 ;parameter k
                 range[*,k]=[intlow,p0[k],intupp]
              endif ;endif 5
           endif ;endif m not first or last
        endfor ;endfor loop over parameters
;     endif ;endif dosys

        if n_elements(whpar) gt 1 then begin 
           print, 'whpar loop'
           whpar=whpar[1:*] ;remove leading 0
           print,whpar
           r=where(range[0,*] ne 0) ;find ranges only for used parameters, is equal to whpar
;           range=range[*,r]
;           wr=where(range[0,*] lt 10 and ppix eq 1)
           wr=[4,5] ;add par 26 if also getting flux errors?
           ;print, 'range[0,wr] = ',range[0,wr]
           syserrlow=sqrt(total((range[0,wr]-range[1,wr])^2)) ;get quadratic sum of error bars on only specified (wr) parameters
           ;should these really be the same for every parameter? (e.g. sigR and r0) maybe, it is the systematic error
           syserrupp=sqrt(total((range[2,wr]-range[1,wr])^2))
           w5=where(r eq 5,nw5) ;check if par5=r0 is included
           if nw5 eq 0 then begin ;if r0 not included, add syserr+err
              print, 'nw5 loop'
              syserrlow=sqrt(syserrlow^2+lowerr^2)
              syserrupp=sqrt(syserrupp^2+upperr^2)
           endif 
           
        endif else begin ;if no systematic terms (whpar empty)
           syserrlow=lowerr ;set systematic error = 0 (err = regular err only)
           syserrupp=upperr
        endelse 
        print, 'lowerr,upperr = ',lowerr,upperr        
        
;        stop
        perror[0,wpar[i]]=syserrlow;lowerr
        perror[1,wpar[i]]=syserrupp;upperr
     endif ;endif go


     ;reset the delchi array so it doesn't erroneously skip the iteration loop for the next parameter -> KAF
     nvar = 20.0
     delchi=dblarr(nvar) & chisq=dblarr(nvar)
     psave=dblarr(n_elements(p0),nvar) ;reset psave and chisq array in case they changed size -> KAF
     
  endfor ;endfor loop over wpar (getting error on wpar[i])

;  print, 'perror = ',perror

  com='zfit='+model+'(x,y,bestfit)'
  tmp=execute(com)
;  !p.multi=0
;;;;implement new loop for if better chisq is found then start over  
;  stop
  return
end 
