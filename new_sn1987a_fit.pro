pro new_sn1987a_fit,im,a,sigma,chisq,dof,iter,fitim,perror,lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat, lobsrc=lobsrc,doplot=doplot,manual=manual,fixedsigr=fixedsigr

  if keyword_set(doplot) then doplot=1 else doplot=0
  if n_elements(fixedsigr) eq 0 then fixedsigr=0

  a0=0.
  a1=1.
  x0=0.
  y0=0.
;  y0=50.
;  sigr=10.
;  sigr=13.
  sigr=6.0
;  sigr=2.
;  sigr=8.
;  r0=50.
;  r0=40.
  if fixedsigr eq 1 then sigr = 8.45;16.9 ;~effective psf width in deconvolved images
;  r0=56.3 
  r0=53.
;  r0 = 45.0
  s0=1.
  sigt0=20.*!dtor
  sigr0=2.
  s1=1.
  sigt1=20.*!dtor
  sigr1=2.
  s2=1.
  sigt2=20.*!dtor
  sigr2=2.
  s3=1.
  sigt3=20.*!dtor
  sigr3=2.
 
  
  theta=[45,135,225,315]*!dtor
  if keyword_set(bilat) then theta=[0,180,225,315]*!dtor
  theta0=theta[0]
  theta1=theta[1]
  theta2=theta[2]
  theta3=theta[3]
  
  xsrc = 0.
  ysrc = 0.
  wsrc = 6.78
;  wsrc = 5.
  fsrc = 1.
;  fsrc = .35
  ; point source counts = 2pi*f*wsrc^2
  
  counts=total(im)
  
;  a =[a0,a1,x0,y0,sigx,sigy,b1,r0,theta0,sigr,sigt,c0,c1,phi0,counts,mu0,sigt2,b2]
  a=[a0,a1,x0,y0,sigr,r0,s0,theta0,sigt0,sigr0,s1,theta1,sigt1,sigr1,$
     s2,theta2,sigt2,sigr2,s3,theta3,sigt3,sigr3,counts,$
     xsrc, ysrc, wsrc, fsrc]*1d
  
;  new_sn1987a_model,x,a,fim
  
;  stop
  parname=['a0','a1','x0','y0','sigr','r0',$
           's0','theta0','sigt0','sigr0',$
           's1','theta1','sigt1','sigr1',$
           's2','theta2','sigt2','sigr2',$
           's3','theta3','sigt3','sigr3','counts',$
           'xsrc', 'ysrc', 'wsrc', 'fsrc']
  
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0], tied:'',$
                       parname:''}, n_elements(a))
  
 if not keyword_set(manual) then begin ;default initial values and limits
  ; setting parameter limits
  parinfo.parname=parname
  parinfo.value=a
  parinfo[22].fixed=1 ;;counts
  parinfo[[2,3,5]].limited=[1,1] ;x0,y0,r0 
  if fixedsigr eq 0 then parinfo[4].limited=[1,1] else parinfo[4].fixed=1;sigr (KAF)
  
;  parinfo[[2,3,5,7,11,15,19]].limited=[1,1]
  
;  parinfo[[0,1]].limited=[1,0] ;a0 & a1
;  parinfo[[0,1]].limits=[0,1]
  
;  parinfo[[6,10,14,18]].limited=[1,0]
;  parinfo[[6,10,14,18]].limits=[0,1e100]
  parinfo[[2,3]].limits=[-30,30] ;x0&y0 ;original
;  parinfo[[2,3]].limits=[-60,60] ;x0&y0
  parinfo[5].limits=[30,80]                            ; r0
  if fixedsigr eq 0 then parinfo[4].limits=[1.0,100.0] ;sigr (KAF)
  parinfo[[8,12,16,20]].limited=[1,1] ;;sigt
;  parinfo[[8,12,16,20]].limits=[20,90]*!dtor ;original
  parinfo[[8,12,16,20]].limits=[10,100]*!dtor

  parinfo[[6,10,14,18]].limited=[1,0] ;s0,s1,s2,...
  parinfo[[6,10,14,18]].limits=[0,0]
  parinfo[[7,11,15,19]].limited=[1,1] ;;theta
;  parinfo[[7,11,15,19]].limits=[[10,80],[100,170],[190,260],[280,350]]*!dtor ;;theta
;  parinfo[[7,11,15,19]].limits=[[0,90],[90,180],[180,270],[270,360]]*!dtor
;  parinfo[[7,11,15,19]].limits=[[0,100],[80,190],[170,280],[260,370]]*!dtor ;original
  parinfo[[7,11,15,19]].limits=[[0,110],[50,230],[140,310],[230,400]]*!dtor
  parinfo[23].tied = 'p[2]' ;fixing xsrc,ysrc=x0,y0 (?) ;original
  parinfo[24].tied = 'p[3]' ;original
;  parinfo[[23,24]].limited = [1,1] ;kaf
;  parinfo[[23,24]].limits = [-30,30] ;kaf
  parinfo[25].fixed = 1 ;wsrc ;original
;  parinfo[25].limited = [1,1] ;wsrc ;kaf
;  parinfo[25].limits = [2,7] ;wsrc ;kaf
  parinfo[26].limited = [1,0] ;fsrc ;original
;  parinfo[26].limited = [1,1] ;fsrc
  parinfo[26].limits = [0,8900] ;original
;  parinfo[26].limited = [1,1] ;kaf
;  parinfo[26].value=0.5
;  parinfo[26].limits = [0,1.0] ;kaf
endif else begin
   ;; fix all parameters except point source at best fit values from
;; lobes only fit
  ddof=4
;  parinfo[0:25].limited=[0,0]
  parinfo.parname=parname
  a0=1.18567744908752E-07
  a1=6604960.71917491
  x0=-4.23091015444864
  y0=-4.39669496656016
  sigr=13.4195548420504
  r0=56.180251776426
  s0=1.31879173578655
  theta0=0.64101527825816
  sigt0=0.857850289099333
  sigr0=2
  s1=1.26293450226898
  theta1=2.65779339626932
  sigt1=0.426513988466946
  sigr1=2
  s2=0.583054794074617
  theta2=4.06173926957818
  sigt2=0.349065840244293
  sigr2=2
  s3=1.52782883522185
  theta3=4.97574074363122
  sigt3=0.349065840244293
  sigr3=2
  counts=5628.7640065433
  xsrc=-4.23091015444864
  ysrc=-4.39669496656016
  wsrc=5.0
  fsrc=0.25
  a=[a0,a1,x0,y0,sigr,r0,s0,theta0,sigt0,sigr0,s1,theta1,sigt1,sigr1,$
     s2,theta2,sigt2,sigr2,s3,theta3,sigt3,sigr3,counts,$
     xsrc, ysrc, wsrc, fsrc]*1d
  parinfo.value=a
  parinfo[0:25].fixed=1
  parinfo[26].limited=[1,1]
  parinfo[26].limits=[0,500.]
endelse

;  weights=1./im^2
  weights=1./im ;;Poisson weights
;  weights=sqrt(im)
  err=sqrt(im) ;Poisson error
  wbad=where(im EQ 0,nwbad) ;setting errors for zero-value pixels
  if nwbad gt 0 then err[wbad]=0.01
;  IF nwbad GT 0 THEN weights[wbad]=0.
  
;  xpos=lindgen(200,200L)
;  xpos=reform(xpos,200L*200L)
;  image = reform(im,200L*200L)
  xr=indgen(200)
  yc=indgen(200)
  xpos=xr#(yc*0+1)
  ypos=(xr*0+1)#yc
  
  if keyword_set(lobes) then model='sn1987a_model_lobes'
  if keyword_set(ring) then model='sn1987a_model_ring'
  if keyword_set(ptsrc) then model='sn1987a_model_ptsrc'
  if keyword_set(bilat) then model='sn1987a_model_2lobes'
  if keyword_set(lobsrc) then model='sn1987a_model_lobes_ptsrc'

  if keyword_set(bilat) then parinfo[[7,11]].limits=[[-90,90],[90,270]]*!dtor
  
  if not keyword_set(manual) then begin
  case model of
     'sn1987a_model_lobes': begin
        parinfo[[9,13,17,21]].fixed=1 ; fix sigr0..sigr3 (lobe radial widths)
        ddof=18
     end 
     'sn1987a_model_ring': begin 
        parinfo[6:21].fixed=1 ; fix all lobe parameters
        parinfo[23:26].fixed=1 ; fix all ptsrc parameters
        ddof=6
     end 
     'sn1987a_model_ptsrc': begin 
        parinfo[[8,12,16,20,23,24,25,26]].fixed=1  ; fix all 4 sigt and all ptsrc (?) parameters
        ddof=18
     end 
     'sn1987a_model_2lobes': begin 
        parinfo[[9,13,14,15,16,17,18,19,20,21]].fixed=1 ; fix all 4 lobe sigr,
                                ; fix s,theta,sigt for lobes 2 and 3
         parinfo[23:26].fixed=1 ; fix all ptsrc parameters
        ddof=13
     end 
     'sn1987a_model_lobes_ptsrc': begin 
        parinfo[[9,13,17,21]].fixed=1  ; fix sigr0..sigr3 (lobe radial widths)
        ddof=19
     end 
     
  endcase
  endif

 ;;; Svet's test
;  parinfo[4].fixed=1
  ;stop
  result=mpfit2dfun(model,xpos,ypos,im,err,a,$;sigma,$
                    /noder,niter=iter,nprint=10,$
                    bestnorm=chisq,ftol=0.0001,$
                    perror=sigma,parinfo=parinfo,$
                    dof=dof,yfit=fitim,maxiter=500);,weights=weights)

;  a=result
;  result=yfit

;  result=mpcurvefit(xpos,image,weights,a,sigma,$
;                    FUNCTION_name=model,$
;                    /noderivative,iter=iter,nprint=10,$
;                    chisq=chisq,ftol=0.0001,$
;                    sigma=sigma,parinfo=parinfo)
  
;  fitim=reform(result,200L,200L)
;  wpar=where(a-result ne 0.,npar)
;  wim=where(im gt 0,nim)
;  dof2=nim-npar
;  sigma=sigma;*sqrt(chisq/dof2)
  print, 'r0sigma = ',sigma[5]
;  print, 'fsrcsigma = ',sigma[26]
  a=result
  parinfo.value=result
;  parinfo[5].limits=[0,100] ;r0
;  sigma[5]=0.5
  w=where(im ne 0 and fitim ne 0 and err ne 0,nw)
  dof=nw-ddof
  chisq=total((im[w]-fitim[w])^2/err[w]^2)

;  parinfo[26].limits=[0.0,5.0]
;  sigma[26]=parinfo[26].value/2.0
;  sigma[26]=parinfo[26].limits[1]
;  sigma[26]=sigma[26]*SQRT(chisq/dof)

  parmins=parinfo.value
  wlimited=where(parinfo.limited[0] eq 1,nwlimited) ;check for limited lower ranges
  if nwlimited gt 0 then parmins[wlimited]=parinfo[wlimited].limits[0]
;  parmins[26]=0.0
;  print, 'parmins = ',parmins

;  stop

;; calculate confidence interval (upper and lower erros, 90% default)
;; on specified parameters (calculate for more params by including in
;; wpar argument)
  
;  conf_error2d,xpos,ypos,im,err,a,sigma,model,perror,bestfit,pvarunit,bestchisq,yfit,psave,/d2,parinfo=parinfo,wpar=[4,5],weights=weights,doplot=doplot,delchi0=1.
errorpars=[5]
if fixedsigr eq 0 then errorpars=[errorpars,4]
if keyword_set(lobsrc) then errorpars=[errorpars,26]
conf_error2d,xpos,ypos,im,err,a,sigma,model,perror,bestfit,pvarunit,bestchisq,yfit,psave,/d2,parinfo=parinfo,wpar=errorpars,weights=weights,doplot=doplot,delchi0=1.,pmin0=parmins

  
  return
end 
