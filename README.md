# 3D_Identification
The code of  3D Identification

# Enviroment
Linux16.04/Windows10/MacOS10.15.5
mex Compiler : Clang
MatLab : R2020b

# Files
We will introduce several directories.
- core: the core function about alignment and split strategies.
- example: some examples of how to use our method.
- gaussian: the code of resisting gaussian noise.
- mex: some C++ dependencies which need to be compiled before using our method.
- model: some models provided for tests.
- retrieval: the code used to retrieve 4 dataset.
- segmentation: the code showing the effectiveness of segmentation strategies.
- thirdparty: the thirdparty code used in our method.
- util: some utilization function used in our method such as reading mesh.
- watermarking: the code used to test watermarking models.
- setup.m: used to prepare environment: compile dependencies.


# Run
To run the code in examples, you should call `setup.m` firstly.
Then just use the .m file in fold example.
For example, input `dragon_rh` command in matlab window.

## Retrieval
Each sub-directory in /retrieval has two files. One is used to prepare data (generate cvs), and the other is used to retrieve dataset.

##
We did some extra experiments on PointNetVLAD Datasets. And the code is in the file /util/retrieve_residual.m. You should download the benchmark dataset first when try to do it. And you should modify the dataset path in /util/retrieve_residual.m.