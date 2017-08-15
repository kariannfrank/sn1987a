@conf_error2d

pro new_sn1987a_image_analysis_single,file,count,lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat, lobsrc=lobsrc ,all=all,nomin=nomin,fixedsigr=fixedsigr

;; runs image fitting on only a single image file (e.g. to be called
;; by new_sn1987a_image_analysis_parallel.pro)

;  for k=0,nf-1 do begin 
;     file=files[k]
;     count=counts[k]

if (count ge 300) or (keyword_set(nomin)) then begin ; KAF -> add check to ensure there are enough counts to do the fitting (unless nomin is set, then always fit)

     if keyword_set(all) then n=4 else n=0
     for i=0,n do begin 
        if keyword_set(all) then begin
           ring = 0 & lobes =0 & ptsrc=0 & bilat = 0 & lobsrc=0
           case i of
              0: ring=1
              1: lobes=1
              2: ptsrc=1
              3: bilat=1
              4: lobsrc=1
           endcase
        endif

        if keyword_set(ring) then mo='ring'
        if keyword_set(lobes) then mo='lobes'
        if keyword_set(ptsrc) then mo='ptsrc'
        if keyword_set(bilat) then mo='bilat'
        if keyword_set(lobsrc) then mo='lobsrc' 
        add='_'+mo
        if not keyword_set(fixedsigr) then sdir=mo+'/' else sdir=mo+'_fixedwidth/'

        if not exist(sdir) then spawn,'mkdir '+sdir
        nn=1
        fpars=create_struct('a0',0d,'a0err',0d,'a1',0d,'a1err',0d,$
                            'x0',0d,'x0err',0d,'y0',0d,'y0err',0d,$
                            'sigr',0d,'sigrerr',dblarr(2),'r0',0d,'r0err',dblarr(2),$
                            's0',0d,'s0err',0d,'theta0',0d,'theta0err',0d,$
                            'sigt0',0d,'sigt0err',0d,'sigr0',0d,'sigr0err',0d,$
                            's1',0d,'s1err',0d,'theta1',0d,'theta1err',0d,$
                            'sigt1',0d,'sigt1err',0d,'sigr1',0d,'sigr1err',0d,$
                            's2',0d,'s2err',0d,'theta2',0d,'theta2err',0d,$
                            'sigt2',0d,'sigt2err',0d,'sigr2',0d,'sigr2err',0d,$
                            's3',0d,'s3err',0d,'theta3',0d,'theta3err',0d,$
                            'sigt3',0d,'sigt3err',0d,'sigr3',0d,'sigr3err',0d,$
                            'counts',0d,$
                            'xsrc', 0d, 'xsrcerr', 0d, 'ysrc',0d, 'ysrcerr', 0d,$
                            'wsrc', 0d, 'wsrcerr', 0d, 'fsrc', 0d,$
                            'fsrcerr', dblarr(2), 'iter',0,'chisq',0d,'dof',0d)
        
        npars=n_tags(fpars)
        print, 'npars = ',npars
;        !p.multi=[0,1,3]
        !p.multi=[0,3,1]        

        im=mrdfits(file,0,inhead)
        im=im-min(im)
        im=im/max(im) ;;;LOSE FLUX INFO
        im=im/total(im)*count
;        PRINT, 'im = ',im
        deproject_image,im,newim
        
        print,total(im),total(newim)
        
        new_sn1987a_fit,newim,a,sigma,chisq,dof,iter,fitim,perror,lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat,lobsrc=lobsrc,fixedsigr=fixedsigr;,/manual
        rdis,im
        rdis,newim
        rdis,fitim

;;     for j=0,n_elements(a)-1 do fpars[i].(j)=a[j]
;        for j=0,npars-3,2 do fpars[i].(j)=a[j/2]
        
;;     dof=200.^2-22.
;;     csigma=sigma*sqrt(chisq/dof)

;        for j=1,npars-2,2 do fpars[i].(j)=sigma[(j-1)/2]

        fnames = tag_names(fpars) ;KAF added loop to account for COUNTS having no errors
        j = 0
        for i = 0,n_elements(a)-1 do begin
;           print, 'j = ',j
;           print, 'fname[j] = ',fnames[j]
           fpars.(j) = a[i]
           if fnames[j] ne 'COUNTS' then begin
              fpars.(j+1) = sigma[i]
              j = j+2
           endif else begin
              j = j+1
           endelse 
        endfor

        r0=a[5]
        r0err=perror[*,5]
        fpars.r0err=r0err
        fpars.sigrerr=perror[*,4]
        fpars.fsrcerr=perror[*,26]
;        print, 'fsrcerr = ',perror[*,26]
;        print, 'fsrcsigma = ',sigma[26]
        fpars.iter=iter
        fpars.chisq=chisq
        fpars.dof=dof

;  ofile='new_sn1987a'+add+'_fitim.fits'
;        root=str_sep(file,'.fits')
;        root=str_sep(root,'/') ;remove path to get base name
        root=FILE_BASENAME(file,'.fits')
        ofile=root+add+'_fitim.fits'
        mwrfits,newim,sdir+ofile,/create
        mwrfits,fitim,sdir+ofile
        diff=newim-fitim
        mwrfits,diff,sdir+ofile

        ;-write fit results-
        fparfile=root+add+'_fpars.fits'
        mwrfits,fpars,sdir+fparfile,/create
        print,r0,r0err

       ;-add header with keywords for input number of counts and input file (KAF)-
        head = headfits(sdir+fparfile)  ;read in header
        sxaddpar,head,'INCOUNTS',count  ;add keyword containing number of input counts
;        sxaddpar,head,'INFILE',file ;add keyword containing the input
;        image file
        sxaddpar,head,'INFILE',FILE_BASENAME(file) ;add keyword containing the input image file
        obsid=sxpar(inhead,'OBS_ID') ;get obs_id from input image file
        sxaddpar,head,'OBS_ID',obsid ;add obsid to fpars file
        modfits,sdir+fparfile,0,head    ;modify the primary header to include the updates

        !p.multi=0
     endfor 
;  endfor 
  endif else begin ;if there were not enough counts to fit
     print, "WARNING: Less than 300 counts in "+file+" -- image fitting skipped."
  endelse 

  return
end

