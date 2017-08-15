;from E. A. Helder

;reads model parameters from fitting sn1987a image to lobe+torus model
;  with new_sn1987a_image_analysis.pro

;modified by K. A. Frank to accept input parameter file as argument,
;as well as model

@modelpileup.pro
@tops.pro
@conf_error2d.pro

pro pileupimage, infile,outfile=outfile,model=model
   ;procedure to turn the modeled parameters of the torus 
   ;into an image with 'infinite' resolution. 



IF N_ELEMENTS(infile) EQ 0 THEN MESSAGE, "ERROR: must provide input parameter file."
IF NOT KEYWORD_SET(model) THEN model='lobes'

file = model+'/'+infile
;   file = 'lobes/sn1987a27_300-8000_cr_lobes_fpars.fits'
;IF N_ELEMENTS(outfile) EQ 0 THEN outfile = 'lobes/model_image.fits' ELSE outfile='lobes/'+outfile
IF N_ELEMENTS(outfile) EQ 0 THEN outfile = model+'/model_image.fits' ELSE outfile=model+'/'+outfile
;   outfile ='model_pileup_13238.fits' 

   tab = mrdfits(file, 1)

  
 a=[tab.a0,tab.a1,0d,0d,tab.sigr,tab.r0,tab.s0,tab.theta0,tab.sigt0,tab.sigr0,tab.s1,tab.theta1,tab.sigt1,tab.sigr1,tab.s2,tab.theta2,tab.sigt2,tab.sigr2,tab.s3,tab.theta3,tab.sigt3,tab.sigr3,tab.counts,tab.xsrc, tab.ysrc, tab.wsrc, tab.fsrc]*1d 

   parname=['a0','a1','x0','y0','sigr','r0',$
           's0','theta0','sigt0','sigr0',$
           's1','theta1','sigt1','sigr1',$
           's2','theta2','sigt2','sigr2',$
           's3','theta3','sigt3','sigr3','counts',$
           'xsrc', 'ysrc', 'wsrc', 'fsrc']
  
   parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0], tied:'',$
                       parname:''}, n_elements(a))
  
   parinfo.parname=parname
   parinfo.value=a
   parinfo[22].fixed=1 ;;counts
   parinfo[1].fixed=1 ;a1, rest schaalt mee..
   parinfo[[0,1]].limited=[1,0]
   parinfo[[0,1]].limits=[0,23] ;a0, a1
   parinfo[[2,3,5]].limited=[1,1]

   parinfo[[2,3]].limits=[-30,30] ;x0&y0
   parinfo[4].limited=[1,1]
   parinfo[4].limits=[5,30] ;width ring
   parinfo[4].fixed = 0
   parinfo[5].limits=[30,80]     ; r0
   parinfo[[8,12,16,20]].limited=[1,1] ;;sigt
   parinfo[[8,12,16,20]].limits=[20,90]*!dtor
   
   parinfo[[6,10,14,18]].limited=[1,0] ;s0,s1,s2,...
   parinfo[[6,10,14,18]].limits=[0,0]
   parinfo[[7,11,15,19]].limited=[1,1] ;;theta
;   parinfo[[7,11,15,19]].limits=[[10,80],[100,170],[190,260],[280,350]]*!dtor ;;theta
;   parinfo[[7,11,15,19]].limits=[[0,90],[90,180],[180,270],[270,360]]*!dtor
   parinfo[[7,11,15,19]].limits=[[0,100],[80,190],[170,280],[260,370]]*!dtor
   parinfo[23].tied = 'p[2]'
   parinfo[24].tied = 'p[3]'
   parinfo[25].fixed = 1
   parinfo[26].limited = [1,0]
   parinfo[26].limits = [0,1200]
 
   dum = (float(indgen(47)) - 23.0) - 0.5
   xpos = rebin(dum,47,47,/sample)
   ypos = transpose(xpos)

 
   makegrid
  
   fitmod = modelpileup(xpos, ypos, a,model=model)
   fitmod = fitmod/total(fitmod)
   mwrfits, fitmod, outfile, /create

stop


end
