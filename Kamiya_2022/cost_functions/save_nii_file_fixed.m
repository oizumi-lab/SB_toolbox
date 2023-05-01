function save_nii_file_fixed(savename, data_array, num_parcels)
% allocate ratio to all voxels and get nii file of map_all
addpath(genpath('/Users/shunsukekamiya/Desktop/1_Lab/6_toolboxes/HCP_Toolbox')); % We use Schaefer's code
% convert ratio data to allocate to all voxels
labels = ft_read_cifti(['/Users/shunsukekamiya/Desktop/1_Lab/6_toolboxes/HCP_Toolbox/clustering/Schaefer2018_100Parcels_7Networks_order.dlabel.nii'],'mapname','array');
dlabel = labels.dlabel;
allocated_data = zeros(32492*2,1);
for i = 1:num_parcels          
   indices = find(dlabel==i);
   allocated_data(indices) = data_array(i);
end 
labeltmp = ft_read_cifti(['/home/kamiya/toolboxes/HCP_Toolbox/clustering/Schaefer2018_cifti/Schaefer2018_100Parcels_7Networks_order.dlabel.nii'],'mapname','array');
labeltmp.data = allocated_data;
ft_write_cifti(savename, labeltmp,'parameter','data');
end
