For the TRF analysis (strf folder)
 0. cz_strf_pipeline is the main script to run TRF analysis.
Before that, data was prepared with cz_produce_data which prepares data in the format required for the TRF function and also prepares predictor matrices using cz_make_pred_mat1. Any change in predictors should be included in the cz_make_pred_mat1.  
 
 2. All models are run at the cluster using qsubcell fun, which is a Fieldtrip function
 see https://www.fieldtriptoolbox.org/faq/how_to_get_started_with_distributed_computing_using_qsub/
 and https://www.fieldtriptoolbox.org/tutorial/distributedcomputing_qsub/
 
 3. Model 6 is the largest model and takes longer to run, so it has it's own scripts so that it can be run in parallel sessions with the 1 to 5 models
