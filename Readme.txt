The code is related to the paper: XXXXX (will be updated soon)

The code is used to run three-dimensional (3D) and two-dimensional (2D) order index analysis on images from cytoskeletal structures.

1. Software
The code runs in the MATLAB programming language. To install, download MATLAB from: http://www.mathworks.com/products/matlab/
2. Data format
Theoretically, any format that can be read or loaded by MATLAB is suitable to be tested by the code. We include the example data in the .zip file “Example data”.
3. Run the analysis
Here are totally 5 MATLAB files, where ‘orderindexmain.m’ is the main program, and the others are functions that will be called during the running of the main program. The explanations for the variables within the code, as well as the ideas in organizing each part of the code, have been detailed in the main program. The possible parameters that should be modified accordingly to your data sets have been highlighted as ‘modify x’ (x refers to the numbering).
Generally, the main program can be divided into 6 parts:
1)	Load images to create a 3D stack
Here the 2D images are stacked up and form a 3D stack, which is then used for the 3D orientation determination and 3D order index calculation.
2)	Create the binary mask selecting the fibrous regions
Here a binary mask will be created mainly based on the signal intensity. Only the fibrous regions of the 3D stack will be identified by this mask, which will be used in acquiring the order index map. 
3)	Acquire the voxel-wise 3D orientation
Here the voxel-wise 3D orientation of the 3D stack is acquired. The method is described in our previous papers (Biomed. Opt. Express 6, 2294–2310 (2015); Biomaterials, 116, 34-47 (2017); Biomaterials, 179, 96-108 (2018)).
4)	Calculate the 3D order index 
The assessing window used for generating order index is typical of the same size as that used for the 3D orientation calculation, with its width normally two to three times the size of fiber diameter. Under each window, orientation information from about 5-10 fibers contributes to the order index, and it focuses on local fiber characteristics or architecture within microenvironment. Of course, the window size can be defined personally according to the research question.
5)	Perform post-processing
Here we generate ‘pretty’ images of orientation and order index. To acquire these images, the raw intensity image is used to provide the contrast of fiber features, and the orientation or order index maps are labeled by different colors to show the orientation or order index information. We then save the order index maps one by one which are generated in a 3D format.
6)	Calculate the 2D order index
Besides 3D order index, within the code we also provide the calculation method for the generation of order index map in a 2D context and save it.
