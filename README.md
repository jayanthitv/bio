
# Anchor Clustering Pipeline



## Introduction

The Anchor Clustering Pipeline is a Python script designed to cluster junction sequences of immune repertoire. The pipeline processes and filters the junction sequences, performs Point Packing and BIRCH algorithms, and conducts pairwise comparisons for determining clonal groups by hierarchical clustering. The script can handle nucleotide junction sequences with several user-defined options to customize the analysis. 


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
git clone https://github.com/skylerchang/AnchorClusteringNt/Anchor_Clustering_PreVJ_Code
```
or
```{bash}
git clone https://github.com/skylerchang/AnchorClusteringNt/Anchor_Clustering_PostVJ_Code
```

Requirments of Python packages: 
```{bash}
pip install pandas numpy scipy scikit-learn fastparquet networkx prettytable pyarrow 
```

## Usage
To run the scripts, please make sure put the dataset and all the codes in the same folder, then use the following command in your terminal path to the codes:

```{bash}
python anchor_clustering_pipeline.py -F input_file.tsv -f data_format -d distance_type -t clustering_threshold -s subcluster -b birch_radius -p partial_data -z size_threshold -m normalized_minimum_distance -l linkage_method
```

### Arguments
* -F, --input_file: Input single file of junction sequences for grouping (tab separated format).
* -d, --distance_type: Specifies how to calculate distance (hd for Hamming distance threshold, norm_hd for normalized Hamming distance threshold).
* -t, --clustering_threshold: The distance threshold for clonal grouping (float).
* -s, --vjgrouping: Choice of grouping clusters based on V and J gene segments (v, vj, none).
* -b, --birch_radius: The radius threshold for the BIRCH algorithm (float).
* -f, --fraction_data: The fraction of data to be fitted in the BIRCH algorithm (float).
* -z, --size_threshold: The size threshold for determining when to stop split up (integer).
* -m, --minimum_distance_ratio: The normalized minimum distance for point packing (float).
* -l, --linkage_method: Specifies the linkage method for hierarchical clustering (single, average, complete, ward, centroid, median, weighted).

## Example
```{bash}
python anchor_clustering_pipeline.py -F example_input.tsv -f nt -d norm_hd -t 0.15 -s vj -b 0.5 -p 1 -z 1000 -m 0.6 -l single
```

## Output
The output will be a CSV file with the prefix HDC_ followed by various parameters, such as HDC_0.6_1000_1.2_2.3_0.5_0.15_1_example_input.csv. The file will contain the clustered CDR3 sequences along with their cluster assignments.

Temporary files generated during the process will be deleted after the pipeline is completed.

## License

This project is licensed under the terms of the MIT License.

## Contact
For help or issues using Anchor Clustering, please submit a GitHub issue.
For other communications, please contact hchang02@uoguelph.ca
