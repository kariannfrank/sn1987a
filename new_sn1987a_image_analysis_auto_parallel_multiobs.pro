@conf_error2d

pro new_sn1987a_image_analysis_auto_parallel,lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat, lobsrc=lobsrc ,all=all,band=band

  ;read file containing obsids and counts (created with get_counts.py)
  READCOL,'sn1987a_decon_counts.txt',/SILENT,obsids,counts,FORMAT='A,I'

  ;construct array of image file names (consistent with get_images.csh)
  IF NOT keyword_set(band) THEN band = '300-8000'
  files = STRARR(N_ELEMENTS(obsids))
  FOR f=0,N_ELEMENTS(files)-1 DO BEGIN
     print, obsids[f]
     CASE STRLEN(obsids[f]) OF
        3: obs = '00'+obsids[f]
        4: obs = '0'+obsids[f]
        5: obs = obsids[f]
     ENDCASE
     files[f] = 'obs'+obs+'_'+band+'_smoothed.img'
  ENDFOR

  IF NOT keyword_set(lobes) THEN lobes=0
  IF NOT keyword_set(ring) THEN ring=0
  IF NOT keyword_set(ptsrc) THEN ptsrc=0
  IF NOT keyword_set(bilat) THEN bilat=0
  IF NOT keyword_set(lobsrc) THEN lobsrc=0
  IF NOT keyword_set(all) THEN all=0

  nf=n_elements(files)

  ;regular for loop
;  for k=0,nf-1 do begin 
;     file=files[k]
;     count=counts[k]

     ;call image fitting
;     new_sn1987a_image_analysis_single,file,count,lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat,lobsrc=lobsrc,all=all

;  endfor 
  ;return

  ; parallel loop
  split_for, 0,nf-1, commands=[$
     'file=files[i]',$
     'count=counts[i]',$

     ;call image fitting
     'new_sn1987a_image_analysis_single,file,count,lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat,lobsrc=lobsrc,all=all'],$
             
     varnames=['files','counts','lobes','ring','ptsrc','bilat','lobsrc','all'],$
     nsplit=7
      
end
