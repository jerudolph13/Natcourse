# Natcourse

This github includes the R programs used for the examples in our paper on the natural course. 

  **natcourse_1_data** : reads in the EAGeR data and manipulates as needed for analysis
  
  **natcourse_2_descr** : generates descriptive summaries of the EAGeR data

  **natcourse_3_noboot** : estimates the observed natural course and g-computation estimated natural course (without CI)
  
  **natcourse_4_boot**  :  estimates the natural course effect and ATE (risk difference & 95% CI)
  
  **natcourse_5_results** : summarizes and visualizes results
  
----------------------------------------------------------------------------------------------------------------------------  
  
There are additionally the following depracated programs in the tfix ("time-fixed") directory that are no longer used in the paper:

  **natcourse_tfix_simple** : the simulated example with time-fixed exposure and few confounders 

  **natcourse_tfix_eager** : the time-fixed EAGeR example (BMI exposure and conception outcome) 

  **natcourse_tfix_effect** : based on the above, uses g-computation to compare risk of conception setting all EAGeR women to 	have a healthy weight at baseline against the observed risk and uses bootstrap to get 95% confidence intervals

 
