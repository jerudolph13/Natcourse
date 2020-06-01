# Natcourse

This github includes the R programs used for the examples in our paper on the natural course. 

  **natcourse_tvar_eager** : program to estimate the observed natural course and g-computation estimated natural course in our  time-varying EAGeR example (compliance exposure and hCG confirmed pregnancy outcome)
  
  **natcourse_tvar_boot**  : program to estimate the natural course effect and ATE (RD plus 95% CI) for the above example 
  
----------------------------------------------------------------------------------------------------------------------------  
  
There are additionally the following depracated programs no longer used in the paper:

  **natcourse_tfix_simple** : the simulated example with time-fixed exposure and few confounders 

  **natcourse_tfix_eager** : the time-fixed EAGeR example (BMI exposure and conception outcome) 

  **natcourse_tfix_effect** : based on the above, uses g-computation to compare risk of conception setting all EAGeR women to 	have a healthy weight at baseline against the observed risk and uses bootstrap to get 95% confidence intervals

 
