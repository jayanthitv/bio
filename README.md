
# Anchor Clustering Pipeline



## Introduction

The Anchor Clustering Pipeline is a Python script with integrated C++ codes designed to cluster junction sequences of immune repertoires. The pipeline processes and filters the junction sequences, performs Point Packing and BIRCH algorithms, and conducts pairwise comparisons for determining clonal groups by hierarchical clustering. The script can handle million scaled nucleotide junction sequences with several user-defined options to customize the analysis. 

<img width="613" alt="Picture1" src="https://user-images.githubusercontent.com/35077981/228913743-d0f84b6c-e707-456d-8288-d9532bd61f58.png">


## Requirements

- Python 3.8 or later
- Pandas
- NumPy
- SciPy
- Scikit-learn
- Fastparquet 
- Networkx 
- Prettytable 
- Pyarrow 


## Installation

Clone the repository to your local machine. (Highly recommend Anchor_Clustering_PreVJ_Code)

```{bash}
git clone https://github.com/skylerchang/Anchor_Clustering_Nt/Anchor_Clustering_PreVJ_Code
```
or
```{bash}
git clone https://github.com/skylerchang/Anchor_Clustering_Nt/Anchor_Clustering_PostVJ_Code
```

## Requirements of Python packages: 
```{bash}
pip install pandas numpy scipy scikit-learn fastparquet networkx prettytable pyarrow 
```

## Requirement of g++ Compiler: 
Anchor Clustering requires the `g++` compiler to build. Please follow installation instructions for different operating systems and verify the installation of `g++`. 


## Input data format
The input data file should be a tab-separated file with the following columns in the given order:

- sequence_id
- v_call
- j_call
- junction_length
- junction
- truth (If no ground truth provided, can be added some information else but with truth as column name) 


## Usage
The pipeline can be run using the following command:

```{bash}
python Anchor_Clustering_PreVJ.py -F input_file -d distance_type -t clustering_threshold -s vjgrouping -b birch_radius -f fraction_data -z size_threshold -m minimum_distance_ratio -p population_size -r random_material_rate -l linkage_method
```
or

```{bash}
python Anchor_Clustering_PostVJ.py -F input_file -s vjgrouping -m minimum_distance_ratio -p population_size -r random_material_rate -b birch_radius -f fraction_data -d distance_type -t clustering_threshold -l linkage_method
```

## Arguments
* -F, --input_file: Path to the input data file. Required.
* -s, --vjgrouping: Specifies whether to use V and J gene segments or not. Options are V, VJ, or None. Optional. Default is VJ.
* -m, --minimum_distance_ratio: The normalized minimum distance for point packing. Optional. Default is 0.6.
* -p, --population_size: The population size for determining the initial populations of anchor sets in Point Packing. Optional. Default is 1000.
* -r, --random_material_rate: The random material rate for introducing the number of mutations in Point Packing. Optional. Default is 50.
* -b, --birch_radius: The radius threshold for BIRCH algorithm. Optional. Default is 0.5.
* -f, --fraction_data: The fraction of data to be fitted in BIRCH algorithm. Optional. Default is 1.
* -z, --size_threshold: The size threshold for determining when to stop split up for subcluster. Optional. Default is 1000.
* -d, --distance_type: Specifies how to calculate distance. "hd" refers to Hamming distance threshold, "norm_hd" refers to normalized Hamming distance threshold. Optional. Default is norm_hd.
* -t, --clustering_threshold: The distance threshold for clonal grouping. Optional. Default is 0.12.
* -l, --linkage_method: Specifies the linkage method in hierarchical clustering. Options are single, average, complete, ward, centroid, median, or weighted. Optional. Default is single.


## Usage example
To run the scripts, please make sure put the dataset and all the codes in the same folder, then use the following command in your terminal path to the codes:

You can try the following commands with the default parameters: 

```{bash}
python Anchor_Clustering_PreVJ.py -F input_file.tsv
```

```{bash}
python Anchor_Clustering_Post.py -F input_file.tsv
```

For saving memory consumption of million scale data, please modify data fration (-f) to be trained in BIRCH model (example with fraction 0.1): 

```{bash}
python Anchor_Clustering_PreVJ.py -F input_file.tsv -f 0.1
```

```{bash}
python Anchor_Clustering_Post.py -F input_file.tsv -f 0.1
```

If needs, you can specify the parameter settings as following: 
 
```{bash}
python Anchor_Clustering_PreVJ.py -F input_file.tsv -s V -m 0.6 -p 1000 -r 50 -b 0.5 -f 1 -z 1000 -d norm_hd -l single
```
```{bash}
python Anchor_Clustering_PostVJ.py -F input_file.tsv -s None -m 0.6 -p 1000 -r 100 -b 0.5 -f 1 -z 1000 -d norm_hd -l single
```


## Output
The output will be a CSV file with the prefix AC_Pre_ or AC_Post_ followed by various parameters, such as <br> AC_Pre_VJ_0.6_1000_50_0.5_0.1_1000_0.12_input_file.csv (with the default settings) or <br>
AC_Post_VJ_0.6_1000_50_0.5_0.1_1000_0.12_input_file.csv (with the default settings). <br>

The file will contain the original input file columns along with their cluster assignments.
Temporary files generated during the process will be deleted after the pipeline is completed.


## General Tips:
Anchor Clustering parameters are categorized into two types: those related to Point Packing and those related to BIRCH clustering.
 
The default settings are designed to accommodate a wide range of immune repertoire datasets. In specific cases, it may be necessary to further optimize these parameters as described below.
 
### Specific Cases:
#### Case 1: Large Datasets (Million-scale)
Memory Management: Adjust the -f, --fraction_data parameter to match your system's memory capacity. Lowering the value to, for instance, 0.01, reduces the portion of data used by the BIRCH model during training, conserving memory during multi-processing.
 
#### Case 2: Medium Datasets (10K to 1M)
Stringent Clustering: If precise clustering is required, adjust the -t, --clustering_threshold parameter. Set a lower threshold (e.g., 0.1 for normalized Hamming distance and 1 for Hamming distance) to achieve more granular clustering.
 
#### Case 3: Small Datasets (less than 10K)
* Consider Clonality: For datasets with high similarity, assess whether the data already exhibits clonal relatedness.
* Distance Ratio Adjustment: The -m, --minimum_distance_ratio can be reduced to enhance the selection of genetically similar anchor sequences. For instance, setting it below 0.5 might slow the point packing process but increases similarity sensitivity. A balance is achieved at 0.6, which is suitable for most naturally occurring immune repertoires.
* Size Threshold Tuning: With the -z, --size_threshold, you can opt for a higher threshold or direct pairwise comparisons (setting it to 10K) for a more detailed analysis in smaller datasets.
 
We hope this guide empowers you to effectively utilize Anchor Clustering for your research. For further assistance or to provide feedback, please reach out through our GitHub issues page. 


## License
This project is licensed under the terms of the MIT License.

## Contact
For help or issues using Anchor Clustering, please submit a GitHub issue.
For other communications, please contact hchang02@uoguelph.ca
