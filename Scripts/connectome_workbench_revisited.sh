# Modules to Load

# FSL
module load bio/FSL/6.0.5.1-centos7_64

# Connectome Workbench
module load bio/ConnectomeWorkbench/1.5.0-rh_linux64

# --------------------------------------------------------------------------------------------------------------------------------------

# Use this command:

# wb_command -cifti-parcellate

# What does it do? --> Parcellate a CIFTI file.

# Syntax: wb_command -cifti-parcellate <cifti-in> <cifti-label> COLUMN <cifti-out> 
   
#       <cifti-in> - the cifti file to parcellate
#       <cifti-label> - a cifti label file to use for the parcellation
#       <direction> - which mapping to parcellate (integer, ROW, or COLUMN)
#       <cifti-out> - output - output cifti file
      
# Example:


wb_command -cifti-parcellate /external/rprshnas01/public_datasets/HCP/HCP_S900/100610/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii /external/rprshnas01/netdata_kcni/jglab/Code/libraries_of_others/github/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/Schaefer2018_200Parcels_7Networks_order.dlabel.nii COLUMN /nethome/kcni/hharita/Code/whobpyt/scratch/testing_ConnectomeWB/100610_rfMRI_200Schaefer_7Ntwx_cifti_parcellated.ptseries.nii

# Use this command:

# wb_command -cifti-correlation

# What does it do? --> Generate correlation of rows in a CIFTI file. 

# Syntax: wb_command -cifti-correlation <cifti> <cifti-out>

      <cifti> - input cifti file
      <cifti-out> - output - output cifti file

wb_command -cifti-correlation /nethome/kcni/hharita/Code/whobpyt/scratch/testing_ConnectomeWB/100610_rfMRI_200Schaefer_7Ntwx_cifti_parcellated.ptseries.nii /nethome/kcni/hharita/Code/whobpyt/scratch/testing_ConnectomeWB/100610_rfMRI_200Schaefer_7Ntwx_cifti_parcellated.pconn.nii