@camilazuga, 02.26.25




How predictors are selected in the TRF model:
- the predictor matrix for a specific model is created in line 96 (onwards) of cz_strf_main_server_v3
- to define which predictors to include, it matches the names defined in cz_modelDefinitions with the names of the predictors saved in a dictionary. The dictionary has key-value pairs of predictor name and predictor row position in the predictors matrix.
- the predictors matrix and dictionary are created (once) in the cz_make_pred_mat1 function, which has textgrids as inputs
- the cz_make_pred_mat1 is called within the cz_produce_cdata_v2 script
- cz_produce_cdata_v2 creates a fieldtrip struct with the data, the predictors and additional info, it is also only run once, the output servers as input to the strf main function

How to add a new predictor to a dataset:
- Modify cz_make_pred_mat1.m. to add the predictor to the predictor matrix
 -- Add the predictor to the pred variable (line 134)
 -- Add the predictor name to the dictionary keys in the k variable (line 179)
 -- Add the predictor position to the dictionary values in the v variable (line 191)
 
- Run cz_produce_cdata_v2, 
This scripts calls cz_make_pred_mat1 and saves the output in '\data\P023_SeussAdult\function_outputs', in a directory defined in line 58 (modify as desired).
The .mat files are named cdata_A*.mat. and serve as input to cz_strf_main_server_v3.

- Go to cz_strf_pipeline_v2, section 1, define the models you want to run in the corresponding input 'modelos' (line 17) and additional input parameters. Run.

Done! (hopefully)

