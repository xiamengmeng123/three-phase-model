**********Three-phase-model
######EVIRONMENT
The main toolkits unimol_tools、rdkit、sklearn、pytorch are required to run this model,
see requirement.txt file for details.

######DATASET
1.  mol_data --Used to store mol files downloaded from public databases. 
             
2.  csv/merged_interaction_energy.csv --Contains approximately 60,000 interaction energy  data sets.
    csv/mol_interaction_energy.csv   --Extracted a training dataset of 60,000 interaction energies in SMILES format.
    csv/mol_free_energy.csv    --Extracted a training dataset of 3,000 binding free energies in SMILES format.
3. Remove similar molecular structures
    python delete_duplicat_molecule.py

######RUN  
Here is some code involving the machine learning model 3phasemodel in the following steps:

1.Convert the counted data into data columns containing SMILES features.
    python ste1_gen_mol_2smitarget.py
    
2 Model training and predicting for about 60,000 interaction energy data,model is stored in clf_interaction.
    python interaction_energy_train.py
    python step2_regression_prediction.py

3.Model training and predicting for about 3000 free energy data,model is stored in clf_fep.
    python free_energy_train.py
    python step3_regression_prediction.py

There are also comparisons between models and cluster analysis.
    python unimol_repr_regression_compare.py   && python clustering.py


**********download_mol2 file is method used to download small molecules data from public datasets.

**********example-fep file is an example for FEP calculation
1. ./run_fep_res_complex.sh $molecule_id $partition
2. ./GMX-FEP-com-lambda.sh $partition {0..27} $run_path








