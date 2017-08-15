;Name: survey_analysis.pro
;
;Date: February 28, 2012
;Author: Kari Frank
;
;Purpose: Read in CSV data file from survey and analyze results
;         
;Calling Sequence: 
;         survey_analysis
;      
;
;Input:
;      requires the file 'Physics_Graduate_Student_Survey.csv' to be
;located in ~/Dropbox/ folder
;      
;
;Output: 
;      
;      
; 
;Usage Notes:
;      
;Example:

PRO survey_analysis

;--------check for required arguments and set defaults--------

;--------set up paths--------
general_data_file = '~/Dropbox/Grad_Student_Campaign/Physics_Graduate_Student_Survey_General.csv'
ta_data_file = '~/Dropbox/Grad_Student_Campaign/Physics_Graduate_Student_Survey_TA.csv'
ra_data_file = '~/Dropbox/Grad_Student_Campaign/Physics_Graduate_Student_Survey_RA.csv'
plot_file = '~/Dropbox/Grad_Student_Campaign/Physics_Graduate_Student_Survey_Plots.ps'
stats_file = '~/Dropbox/Grad_Student_Campaign/Physics_Graduate_Student_Survey_Stats.ps'

;--------set up question strings--------
general_question_labels = ["How many years have you been in graduate school (round up)?", $
"I feel the Department has my best interests in mind.", $
"I would recommend the Purdue Physics graduate program to potential graduate students.", $
"I have been provided with opportunities to provide input into decisions that affect my experience as a graduate student.", $
"I have taken advantage of said opportunities to provide input.", $
"I have never experienced disparaging remarks or behavior from physics faculty members."]

ta_question_labels = ["How many years have you held a position as a teaching assistant (round up)?", $
"Have you been an RA and then gone back to being a TA?", $
"What is the largest number of teaching assignments you have held at one time? Anything that requires...", $	
"My supervisors treat me in a courteous, respectful, and professional manner.", $
"My supervisor manages the course effectively, including supervision of TAs (e.g. is well organized,...", $	
"I regularly spend more than 20 hours (10 hours if quarter-time) per week on my teaching duties.", $	
"My duties as a TA are clearly defined.", $	
"How many times have you written exams or exam questions for a course in which you were not the super...", $
"How many times have you written lecture notes (for the main lecture, not recitation or lab) for a co...", $	
"I am provided with adequate resources and equipment to do my job as a TA effectively."]

ra_question_labels = ["How many years have you been doing research (registered for 590 or 699, round up)?", $
"How many research advisors have you worked for?", $
"My research advisor has my best interests in mind.", $
"In what year of graduate school do you feel you really started your Ph.D. research?", $
"My research advisor cares about me as a person.", $
"My research advisor treats me in a courteous, respectful, and professional manner.", $
"My research advisor gives me useful feedback on my performance.", $
"My research advisor is a good communicator (e.g. responds promptly to emails, is available to speak...", $
"My research advisor is an effective group manager.", $
"How many times has your research advisor helped you fill out the Thesis Research Feedback form?", $
"The Thesis Research Feedback form successfully facilitated communication between me and my research advisor.", $
"My advisor has helped me write a CV or resume.", $
"How many conferences or large meetings have you attended?", $
"How many presentations have you given on your research (either oral or poster, including both on and off campus?", $
"My research advisor has introduced me to other physicists in my field.", $
"My research advisor encourages me to interact and/or collaborate with other physicists in my field.", $
"I am confident I will be able to obtain three positive letters of recommendation.", $
"I have been included in the grant or proposal application process.", $
"I have been provided with adequate resources and equipment to do my research effectively.", $
"I feel I am making satisfactory progress toward my degree."]

question_labels = [general_question_labels,ta_question_labels,ra_question_labels]
PRINT, where((question_labels EQ "How many years have you held a position as a teaching assistant (round up)?") OR (question_labels EQ "How many years have you been doing research (registered for 590 or 699, round up)?"))

;--------set up answer strings--------
agree_strs = ['Strongly Agree','Agree','Neutral','Disagree','Strongly Disagree']
numyears_zero_six = ['0 years','1 year','2 years','3 years','4 years','5 years','6+ years']
numyears_one_six = ['1 year','2 years','3 years','4 years','5 years','6+ years']
numyears_one_eight = ['1 year','2 years','3 years','4 years','5 years','6 years','7 years','8+ years']
na_agree_strs = [agree_strs,'Not Applicable']
nogroup_agree_strs = [agree_strs,'No Group']
yesno = ['Yes','No']

;--------read in data--------
READCOL,general_data_file,DELIMITER=';',SKIPLINE=1,/SILENT,Q1,Q2,Q3,Q4,Q5,Q6
READCOL,ta_data_file,DELIMITER=';',SKIPLINE=1,/SILENT,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16
READCOL,ra_data_file,DELIMITER=';',SKIPLINE=1,/SILENT,Q17,Q18,Q19,Q20,Q21,Q22,Q23,Q24,Q25,Q26,Q27,Q28,Q29,Q30,Q31,Q32,Q33,Q34,Q35,Q36


data = INTARR(N_ELEMENTS(Q1),36) ;row=person,column=question
FOR i=0,N_ELEMENTS(Q1)-1 DO BEGIN
   data[i,0]=Q1[i]
   data[i,1]=Q2[i]
   data[i,2]=Q3[i]
   data[i,3]=Q4[i]
   data[i,4]=Q5[i]
   data[i,5]=Q6[i]
   data[i,6]=Q7[i]
   data[i,7]=Q8[i]
   data[i,8]=Q9[i]
   data[i,9]=Q10[i]
   data[i,10]=Q11[i]  
   data[i,11]=Q12[i]
   data[i,12]=Q13[i]
   data[i,13]=Q14[i]
   data[i,14]=Q15[i]
   data[i,15]=Q16[i]
   data[i,16]=Q17[i]
   data[i,17]=Q18[i]
   data[i,18]=Q19[i]
   data[i,19]=Q20[i]
   data[i,20]=Q21[i]
   data[i,21]=Q22[i]
   data[i,22]=Q23[i]
   data[i,23]=Q24[i]
   data[i,24]=Q25[i]
   data[i,25]=Q26[i]
   data[i,26]=Q27[i]
   data[i,27]=Q28[i]
   data[i,28]=Q29[i]
   data[i,29]=Q30[i]
   data[i,30]=Q31[i]
   data[i,31]=Q32[i]
   data[i,32]=Q33[i]
   data[i,33]=Q34[i]
   data[i,34]=Q35[i]
   data[i,35]=Q36[i]
ENDFOR

total_years = data[*,0]
ra_years = data[*,WHERE(question_labels EQ ra_question_labels[0])]
ta_years = data[*,WHERE(question_labels EQ ta_question_labels[0])]

;--------set up plot--------
SET_PLOT, 'ps'
DEVICE, /COLOR
DEVICE, FILENAME=plot_file
LOADCT, 38
colors = INDGEN(9)*256/9
leg_colors_one_six = [colors[0:5]]
leg_colors_one_eight = [colors[0:7]]

;plot histograms by year
FOR p=0L,N_ELEMENTS(question_labels)-1 DO BEGIN ;loop through and plot each question
   ;--plot all years--
   array = data[*,p]
   good_subs = WHERE(array NE -99)
   data_size = N_ELEMENTS(good_subs)

   IF p LT 6 THEN avg_by_year = FLTARR(8) ELSE avg_by_year = FLTARR(6)
   IF p LT 6 THEN median_by_year = FLTARR(8) ELSE median_by_year = FLTARR(6)
   IF p LT 6 THEN error = FLTARR(8) ELSE error = FLTARR(6)

   ;set up answer labels
   CASE question_labels[p] of 
      "How many years have you been in graduate school (round up)?": answers = numyears_one_eight
      "I have taken advantage of said opportunities to provide input.": answers = na_agree_strs
      "How many years have you been doing research (registered for 590 or 699, round up)?": answers = numyears_zero_six
      "How many research advisors have you worked for?": answers = numyears_zero_six
      "In what year of graduate school do you feel you really started your Ph.D. research?": answers = numyears_one_six
      "My research advisor is an effective group manager.": answers = nogroup_agree_strs
      "How many times has your research advisor helped you fill out the Thesis Research Feedback form?": answers = numyears_zero_six
      "My advisor has helped me write a CV or resume.": answers = yesno
      "How many conferences or large meetings have you attended?": answers = numyears_zero_six
      "How many presentations have you given on your research (either oral or poster, including both on and off campus?":  answers = numyears_zero_six
      "How many years have you held a position as a teaching assistant (round up)?": answers = numyears_zero_six
      "Have you been an RA and then gone back to being a TA?": answers = yesno
      "What is the largest number of teaching assignments you have held at one time? Anything that requires...": answers = numyears_zero_six
      "How many times have you written exams or exam questions for a course in which you were not the super...": answers = numyears_zero_six
"How many times have you written lecture notes (for the main lecture, not recitation or lab) for a co...": answers = numyears_zero_six
      else: answers = agree_strs
   ENDCASE 

   IF N_ELEMENTS(good_subs) GT 1 THEN PLOTHIST,array[good_subs],BIN=1,TITLE=question_labels[p];,/BOXPLOT;,XTICKNAME=['',answers,''za]
   legend_labels = STRARR(N_ELEMENTS(answers))
   FOR l=0,N_ELEMENTS(answers)-1 DO BEGIN
      legend_labels[l] = STRING(l+1,FORMAT='(I0)')+'='+answers[l]
   ENDFOR
   LEGEND,legend_labels,/TOP,/RIGHT

   ;loop through each year
   IF (p NE 6) AND (p NE 16) AND (p NE 0) THEN BEGIN ;skip the 'how many years' questions
      IF p LT 6 THEN years = total_years
      IF p GT 6 THEN years = ta_years - 1
      IF p GT 16 THEN years = ra_years - 1
      FOR y = MIN(years[WHERE(years GT 0)]),MAX(years[WHERE(years GE 0)]) DO BEGIN
         year_subs = WHERE((array NE -99) AND (years EQ y))
         IF year_subs[0] NE -1 THEN BEGIN
            IF (MAX(array[year_subs]) NE MIN(array[year_subs])) AND (N_ELEMENTS(year_subs) GT 1) THEN BEGIN
               hist = HISTOGRAM(array[year_subs],BINSIZE=1,LOCATIONS=x)
               OPLOT,x,hist,PSYM=-4,COLOR=colors[y],THICK=3
               ;PLOTHIST, array[year_subs],BIN=1,/OVERPLOT,COLOR=colors[y]
            ENDIF 
         ENDIF 
      ENDFOR ;end for y
      IF p LT 6 THEN LEGEND,numyears_one_eight,TEXTCOLORS=leg_colors_one_eight,/TOP,/LEFT ELSE LEGEND,numyears_one_six,TEXTCOLORS=leg_colors_one_six,/TOP,/LEFT
      
      ;calculate and plot averages and median by year
      IF p LT 6 THEN maxyear=8 ELSE maxyear=6
      IF p LT 6 THEN x = [1,2,3,4,5,6,7,8] ELSE x = [1,2,3,4,5,6]
      FOR y = 1,maxyear DO BEGIN
         year_subs = WHERE((years EQ y) AND (array GE 0))
         IF year_subs[0] NE -1 THEN BEGIN
            avg_by_year[y-1] = MEAN(array[year_subs])
            median_by_year[y-1] = MEDIAN(array[year_subs])
            error[y-1] = (N_ELEMENTS(year_subs))^(-0.5)
         ENDIF ELSE BEGIN
            avg_by_year[y-1] = 0
            median_by_year[y-1] = 0
            error[y-1] = 0
         ENDELSE
      ENDFOR 
      PLOTERROR, x,avg_by_year,error,TITLE=question_labels[p],XTITLE='years experience',YTITLE='Average Answers',THICK=3,PSYM=-4,YRANGE=[1,6],ERRTHICK=3
      ;OPLOTERROR,x,median_by_year,error,THICK=3,PSYM=-5,COLOR=colors[3],ERRCOLOR=colors[3]
      LEGEND,legend_labels,/TOP,/RIGHT
   ENDIF 
ENDFOR

DEVICE, /CLOSE
SET_PLOT,'X'

END
