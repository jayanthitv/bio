import sys
import os
import time
import warnings
import resource
import argparse
import subprocess
import concurrent.futures
import glob
from tkinter import N
import pandas as pd
import numpy as np
from prettytable import HEADER
from sklearn.cluster import Birch
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.exceptions import ConvergenceWarning
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import linkage, fcluster
import tracemalloc
from sqlalchemy import true
import networkx as nx
import pyarrow
import fastparquet


def filter_and_save(small_df,length,name,size_threshold,file_suffix):
    """
    Filter clusters based on size threshold and save the data to either a feather file or a FASTA file.
    
    Args:
        small_df (pandas.DataFrame): A DataFrame containing a subset of junction sequence data.
        length (int): The junction length to filter the data.
        name (str): The name of the gene cluster.
        size_threshold (int): The size threshold for filtering subclusters.
        file_suffix (str): A suffix to add to the output file names.
    """

    small_df_junction = small_df.junction.reset_index(drop=True).unique()
    # Filter clusters based on size threshold
    if len(small_df_junction) > size_threshold:
        small_df.reset_index(drop=True).to_feather(f'seq_{length}_{name}_cluster.feather')
        # Save subclusters to FASTA file format
        with open(f"seq_{length}_{name}.fasta", "w") as fa:
            for line in small_df_junction:
                fa.write(f'>sequence\n{line.strip()}\n')
    else:
        # Save smaller clusters to feather file format
        small_df.reset_index(drop=True).to_feather(f'norm_seq_{length}_{name}_cluster.feather')


def process_length(length, junction_data, size_threshold, subcluster):
    """
    Process the data by filtering and clustering based on junction lengths.
    
    Args:
        length (int): The junction length to filter the data.
        junction_data (pandas.DataFrame): A DataFrame containing junction sequence data.
        size_threshold (int): The size threshold for filtering subclusters.
        subcluster (str): The type of subclustering to perform. Accepts 'V', 'VJ', or any other value as default.
    """

    junction_sub_length = junction_data[junction_data.junction_length == length].reset_index(drop=True)
    
    if subcluster in ['V', 'VJ']:
        # Cluster the data based on V gene segment
        junction_subv = gene_clustering(junction_sub_length, type='v_call')
        groupedv = junction_subv.groupby('cluster')

        for name, group in groupedv:
            small_df = pd.DataFrame(group.iloc[:, :6])
            
            if subcluster == 'VJ':
                # Iterate through the J gene clusters
                junction_subj = gene_clustering(small_df, type='j_call')
                groupedj = junction_subj.groupby('cluster')
                
                for name1, group1 in groupedj:
                    small_df2 = pd.DataFrame(group1.iloc[:, :6])
                    filter_and_save(small_df2, length, f'{name}-{name1}', size_threshold, '')
            else:
                filter_and_save(small_df, length, name, size_threshold, '')
    else:
        filter_and_save(junction_sub_length, length, 'NoneVJ', size_threshold, '')


def parallel_process(junction_data,size_threshold,subcluster):
    """
    Process the data in parallel for each junction length.

    Args:
        junctions_data (pandas.DataFrame): A DataFrame containing junction data.
        size_threshold (int): The size threshold for filtering clusters.
    """

    # Process data in parallel for each junction length
    with concurrent.futures.ProcessPoolExecutor() as executor:
        process_list = [executor.submit(process_length, length, junction_data, size_threshold,subcluster) for length in junction_data.junction_length.unique()]
        concurrent.futures.wait(process_list)



def parallel_process(junction_data,size_threshold,subcluster):
    """
    Process the data in parallel for each junction length.

    Args:
        junctions_data (pandas.DataFrame): A DataFrame containing junction data.
        size_threshold (int): The size threshold for filtering clusters.
    """

    # Process data in parallel for each junction length
    with concurrent.futures.ProcessPoolExecutor() as executor:
        process_list = [executor.submit(process_length, length, junction_data, size_threshold,subcluster) for length in junction_data.junction_length.unique()]
        concurrent.futures.wait(process_list)




def data_filter(junction_data, size_threshold,subcluster):

    """
    Analyzes the junction length in a table of junction sequences.
    
    Args:
        junction_data (pandas.DataFrame): A table of junction sequences.
        size_threshold (int): The minimum size of a subcluster to write to a file.
        
    Returns:
        seq{length}.feather: feather files of varied junction lengths
        seq{length}.fasta: FASTA files of varied junction lengths
        norm_seq_{length}_{name}-{name1}_cluster.feather: Subluster file meets the size threshold and waits to run pairwise comparison.
    """

    # Remove ambiguous sequences (including N)
    print("###Initial junction nucleotide data with total sequence:",len(junction_data))
    junction_data = junction_data[~junction_data['junction'].str.contains('N')] 
    print("###After removing ambiguous sequences(including N):",len(junction_data))
 
    # Write data to feather files for each junction length
    parallel_process(junction_data,size_threshold,subcluster)



def create_coordinate(junction_data,length,name,population_size,min_anchor_distance,random_material_rate,clustername,size_threshold,birch_radius,fraction):
    """Creates a coordinate table of evaluating Hamming Distances between junction sequences and selected anchors

    Args:
        junction_data (pandas.DataFrame): A table of junction sequences.
        length (int): The junction length to analyze.
        population_size (int): The population size for the point packing algorithm.
        min_anchor_distance (int): The minimum distance between anchors.
        random_material_rate (float): The random material rate for the point packing algorithm.
        clustername (str): The name of the cluster.
        size_threshold (int): The minimum size of a cluster to write to a file.
        birch_radius (float): The radius parameter of BIRCH algorithm. 
        fraction (float): The fraction of data to be fitted into birch model.

    Returns:
        norm_seq_{length}_{name}_clusternt_{clustername}_{i}.feather: Subluster file meets the size threshold and waits to run pairwise comparison.
        x_seq_{length}_{name}_clusternt_{clustername}_{i}.feather: Subluster file is above the size threshold and needs to futher split.
    """

    # Generate anchor sequence set with point packing algorithm
    start_pp = time.time()
    subprocess.run(["./a.out",f'seq_{length}_{name}.fasta',f'{length}',f'{name}',f'{population_size}',f'{min_anchor_distance}',f'{random_material_rate}'], stdout=subprocess.PIPE)
    first_line = open(f'best.pcld_{length}_{name}').readline().rstrip()
    num = int(first_line.split(" ")[0])
    if num >= 2:
        end_pp = time.time()
        spend_pp = round((end_pp-start_pp),5)

        # Extract the anchor sequences
        anchors=[]
        with open(f'best.pcld_{length}_{name}') as f:
            for line in f:
                if line.endswith("-fitness\n"):
                    continue
                else:
                    anchors.append([list(x) for x in line.rstrip().split(",")])

        # Prepare junction sequence list 
        seqs= junction_data['junction']
        sequences = []
        for line in seqs:
            sequences.append([list(x) for x in line.rstrip().split(",")])

        # Calculate Hamming distances between sequences to each anchor sequence 
        a=np.array(anchors).squeeze()           
        b=np.array(sequences).squeeze()
        arr = pd.DataFrame((a!=b[:, None]).sum(axis=2),dtype=np.uint8).add_prefix('dis') 
   
        # Normalize Hamming distance coordinate data into scaled numeric data 
        normalized_df = normalize(StandardScaler().fit_transform(arr)) 
        normalized_df = pd.DataFrame(normalized_df,dtype=np.float16)  
        

        tracemalloc.start()
        warnings.filterwarnings("ignore", category=ConvergenceWarning, module="sklearn")

        anc_model = Birch(n_clusters=num, threshold=birch_radius, branching_factor=5000)

        # Fit the model to your data and get the cluster labels
        labels = anc_model.fit_predict(normalized_df.sample(int(len(normalized_df) * fraction), random_state=91291640))

        # Get the number of unique labels (clusters)
        n_clusters = len(np.unique(labels))

        if (n_clusters < num):
            print(f"######!!! BIRCH can only get cluster number is {n_clusters}", f"but Anchor clustering specify cluster number is {num} !!!######")
            # Create a new Birch model with the number of clusters equal to the result of the fit_predict method
            birch_model = Birch(n_clusters=n_clusters,threshold=birch_radius, branching_factor=5000)

            # Fit the model to your data again
            birch_model.fit(normalized_df.sample(int(len(normalized_df) * fraction), random_state=91291640))
  
            yhat = birch_model.predict(normalized_df)
            num=n_clusters

        else: 
             yhat = anc_model.predict(normalized_df)

        yhat = pd.DataFrame(yhat,columns=['cluster'],dtype="category")
        

        # Calculate the memory usage of peak 
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")


        sum_table = pd.concat([junction_data, yhat], axis=1, ignore_index=True)
        grouped = sum_table.groupby(sum_table.columns[-1])
        cluster = globals()


        # Extract the clusters and determine if it needs to be further split
        for i in range(0,int(num)):
            cluster['double_%d' % i] = grouped.get_group(i)
            clustern = cluster['double_%d' % i][cluster['double_%d' % i].columns[0:6]]

            if len(clustern) > size_threshold:
                clustern.columns = ['sequence_id','v_call','j_call','junction_length','junction','truth']
                clustern.reset_index(drop=True).to_feather(f'x_seq_{length}_{name}_clusternt_{clustername}_{i}.feather')
            else: 
                clustern.columns = ['sequence_id','v_call','j_call','junction_length','junction','truth']
                clustern.reset_index(drop=True).to_feather(f'norm_seq_{length}_{name}_clusternt_{clustername}_{i}.feather')
                
    else:
        junction_data.reset_index(drop=True).to_feather(f'norm_seq_{length}_{name}_clusternt_{clustername}.feather')





def point_packing(junction_data_name,size_threshold,min_distance_ratio,population_size,random_material_rate,birch_radius,fraction):
    """Determine if junction data needs to be further processed by point packing algorithm

    Args:
        length (int): The junction length to analyze.
        size_threshold (int): The minimum size of a cluster to write to a file.
        min_distance_ratio (float): a distance ratio equals the minimum distance divided by junction length.
        birch_radius (float): The radius parameter of BIRCH algorithm. 
        fraction (float): The fraction of data to be fitted into birch model.

    Returns:
        norm_seq_{length}_clusternt.txt: Subcluster file meets the size threshold and waits to run pairwise comparison.
    """ 

    junction_data = pd.read_feather(junction_data_name)
    sequences = junction_data['junction']

    length = int(junction_data_name.split("_")[1])
    name = junction_data_name.split("_")[2]
    min_anchor_distance = int(length * min_distance_ratio)
  

    # Creates a coordinate table of evaluating Hamming Distances between junction sequences and selected anchors 
    clustername = ''
    create_coordinate(junction_data,length,name,population_size,min_anchor_distance,random_material_rate,clustername,size_threshold,birch_radius,fraction)
    
    # For subcluster file needs to be further separated:
    files_processed = True
    while files_processed:
        files_processed = False

        for xfile in glob.glob(f"x_seq_{length}*.feather"):
            files_processed = True
            length = int(xfile.split('_')[2]) 
            names = xfile.split('_')[3] 
            clusternamesub =xfile.split('_')[-1]
            clusternamesub = clusternamesub.replace('.feather','')
            sub_data = pd.read_feather(xfile)
            namenew=f"{names}-{clusternamesub}"
            
            # Extract the unique sequences 
            sub_sequences = sub_data['junction']
            unique_sub_sequences = set(sub_sequences)

            # If it's greater than the size threshold, process the data and run point packing and generate coordinate table
            if len(unique_sub_sequences) > size_threshold: 
                with open(f"seq_{length}_{namenew}.fasta", "w") as fasub:
                    for line in unique_sub_sequences:
                        fasub.write(f'>seq\n{line.strip()}\n')
                min_anchor_distance = int(length * min_distance_ratio)
                os.remove(xfile)
                create_coordinate(sub_data,length,namenew,population_size,min_anchor_distance,random_material_rate,clusternamesub,size_threshold,birch_radius,fraction)
            # If it's smaller than the size threshold, the subcluster file will be in queue for pairwise comparison
            else:
                os.rename(xfile,f'norm_{xfile}')

        
        ###### Sleep 5 seconds for sometimes Point Packing has longer time for selecting anchors and causing x_seq_{length}*.feather files left ####### 
        ###### If this case still happened, you can either adjust the size_threshold (greater than 1000) or increase the sleep time ############        
        time.sleep(5)  
        if not files_processed:
            break


def hier_clustering(postdata,max_d,distance,link):
    """
    Hierarchical clustering to determine the clusters from subcluster files

    Args:
        postdata (pandas.DataFrame): A table of junction sequences
        max_d (int or float): A distance as the cut-off to determine the clusters
        distance (str): A method of evaluating distance, normalized Hamming distance or Hamming distance
        link (str): Linkage method

    Returns:
        numpy.ndarray: Cluster results of each subcluster file.
    """ 
    # Get the number of sequences
    N = len(postdata)

    # Initialize the distance matrix
    dist_matrix = np.zeros((N, N))

    # Loop through the sequences
    for i in range(N):
        for j in range(N):
            # Skip if i and j are the same sequence
            if i == j:
                continue
            # Calculate the distance between sequences i and j
            elif j > i:
                a = postdata[i]
                b = postdata[j]
                if distance == 'norm_hd':
                    dist_matrix[i, j] = sum(i != j for i, j in zip(a, b))/len(a) 
                else:
                    dist_matrix[i, j] = sum(i != j for i, j in zip(a, b))
            # Set the distance in the distance matrix and the other direction
            else:
                dist_matrix[i, j] = dist_matrix[j, i]
    # Convert the distance matrix to a distance array
    dist_array = ssd.squareform(dist_matrix)

    # Perform specified linkage method
    Z = linkage(dist_array, link)

    # Get the clusters
    clusters = fcluster(Z, max_d, criterion='distance')

    # Return the clusters
    return(clusters)


def gene_clustering(table,type):
    """
    Further split clusters based on if sequences in one cluster have the same gene segments.
    E.G, A complicated case of inferring seq1, seq2 and seq3 as one cluster (seq1 and seq2 have a shared IGHV3 V gene segment, and seq2 and seq3 have a shared IGHV2)
       
       sequences    v_segment       
       seq1         IGHV3,IGHV4
       seq2         IGHV3,IGHV2
       seq3         IGHV2 

    Args:
        table (pandas.DataFrame): A table of junction sequences with their associated gene families.
        type (str): The type of gene segments to be considered for clustering, e.g., 'v_call' or 'j_call'.

    Returns:
        pandas.DataFrame: A table with updated cluster assignments based on gene segment.
    """


    # Clean up the gene_type column and remove duplicates from the new_segment column
    table['new_segment'] = table[f'{type}'].str.replace('\\*..','',regex=True)
    table.new_segment = table.new_segment.str.replace('D','',regex=True)
    table.new_segment = table.new_segment.apply(lambda x: pd.Series(str(x).split(','))).replace(regex=r'(.*?)-(.*?)-(.*)',value='\\1-\\2').apply(lambda x: ','.join(x.dropna()),axis=1)
    table.new_segment= table.new_segment.apply(lambda x: x.split(',')).apply(lambda x: list(set(x))).str.join(',')

    # Explode the new_segment column
    table_exp=table.assign(new_segment=table['new_segment'].str.split(',')).explode('new_segment')
    temp=table_exp[['junction', 'new_segment']]
    g = nx.from_pandas_edgelist(temp, source='junction', target='new_segment')

    # get clusters
    clusters = list(nx.connected_components(g))

    # get membership
    membership = {}
    for i, c in enumerate(clusters):
        for node in c:
            membership[node] = i

    # convert membership to dataframe
    membership_df = pd.DataFrame(list(membership.items()), columns=['junction', 'cluster'])

    # join with original dataframe
    table = pd.merge(table, membership_df, on='junction')
    table.drop(columns=['new_segment'],inplace=True)

    return(table)



def pairwise_process(file,cut_off,distance,link):
    """
    Perform pairwise comparisons on subclusters with hierarchical clustering and further split clusters.
    
    Args:
        file (str): The filename of a subcluster table of junction sequences (in .feather format).
        cut_off (int or float): A distance as the cut-off to determine the clusters.
        distance (str): A method of evaluating distance, normalized Hamming distance or Hamming distance.
        subcluster (str): A choice of further separate clusters based on V and J gene segments.
        link (str): Linkage method for hierarchical clustering.
    Returns:
        clu_{file}.txt: cluster results of each subcluster file.

    """ 


    # Read in the file as a pandas dataframe
    tablesub = pd.read_feather(file)
    # Get the number of unique junctions
    lines = len(tablesub.junction.unique()) 

    if lines!=1:
        start_pw = time.time()
        clusters = hier_clustering(tablesub.junction.unique(),cut_off,distance,link)
        end_pw = time.time()
        spend_pw = round((end_pw-start_pw),5)

        # Create a dataframe from the clusters
        clone = pd.DataFrame(clusters,columns=['cluster']).reset_index()
        # Add an index column to the tablesub dataframe
        tablesub['index'] = tablesub.groupby('junction',sort=False).ngroup()
        # Merge the tablesub and clone dataframes on the index column
        tablesub = pd.merge(tablesub, clone, how="outer", on=["index"]).drop(columns=['index'])
        tablesub=tablesub.sort_values(by=['cluster'])
        tablesub.reset_index(drop=True).to_feather(f'clu_{file}.feather')
    else:
    # Add a cluster column to the tablesub dataframe
        tablesub['cluster'] = '0'
        tablesub.reset_index(drop=True).to_feather(f'clu_{file}.feather')



def ArgParserCommands():
    print (''' Running Anchor Clustering PreVJ Pipeline:''')
    parser = argparse.ArgumentParser(description='Running Anchor Clustering pipeline for junction nucleotide data. Please provide the data file to be clustered.' 
                                    '''
                                    The dataset should be separated by tab delimeter and the column names and orders as following: 
                                    sequence_id v_call  j_call  junction_length    junction    truth
                                    '''
                                    )

    parser.add_argument('-F', '--input_file', dest='infile', default='', required=True,
                        help='Input single file of junction sequences for clustering.') 
    parser.add_argument('-d', '--distance_type', dest='distance', required=False, choices=('hd', 'norm_hd'), default='norm_hd',
                        help='Specifies how to calculate distance. "hd" refers to Hamming distance threshold, "norm_hd" refers to normalized Hamming distance threshold.')
    parser.add_argument('-t', '--clustering_threshold', dest='thres', type=float, required=False, default=0.12,
                        help='The distance threshold for clonal grouping.')
    parser.add_argument('-s', '--vjgrouping', dest='sub', type=str, required=False, choices=('V', 'VJ', 'None'), default='VJ',
                        help='Default is using V and J gene segments. If "none", this program will only use junction sequences and omit the V and J gene information.')
    parser.add_argument('-b', '--birch_radius', dest='radius', type=float, required=False, default=0.5,
                        help='The radius threshold for BIRCH algorithm.')
    parser.add_argument('-f', '--fraction_data', dest='fraction', type=float, required=False, default=1,
                        help='The fraction of data to be fitted in BIRCH algorithm.')
    parser.add_argument('-z', '--size_threshold', dest='size', type=int, required=False, default=1000,
                        help='The size threshold for determining when to stop split up for subcluster.')
    parser.add_argument('-m', '--minimum_distance_ratio', dest='mini', type=float, required=False, default=0.6,
                        help='The normalized minimum distance for point packing')
    parser.add_argument('-p', '--population_size', dest='population', type=int, required=False, default=1000,
                        help='The population size for determining the initial populations of anchor sets in Point Packing.')
    parser.add_argument('-r', '--random_material_rate', dest='random', type=int, required=False, default=50,
                        help='The random material rate for introducing the number of mutations in Point Packing.')
    parser.add_argument('-l', '--linkage_method', dest='link', required=False, choices=('single', 'average', 'complete', 'ward', 'centroid', 'median', 'weighted'), default='single',
                        help='Specifies the linkage method in hierarchical clustering.')

    return parser.parse_args()


def main():
    parsed_args = ArgParserCommands()
    distance = parsed_args.distance
    cut_off = parsed_args.thres
    birch_radius = parsed_args.radius
    fraction = parsed_args.fraction
    size_threshold = parsed_args.size
    min_distance_ratio = parsed_args.mini
    population_size = parsed_args.population
    random_material_rate = parsed_args.random
    subcluster = parsed_args.sub
    link = parsed_args.link



    print("Processing and filtering junction sequences")
    File = parsed_args.infile
    junction_data = pd.read_table(File, delimiter='\t')
    # len_list = list(junction_data.junction_length.unique())
    tstart_pr = time.time()
    data_filter(junction_data, size_threshold, subcluster)
    tend_pr = time.time()
    tspend_pr = round((tend_pr-tstart_pr),3)
    print(f"total processing VJ is {tspend_pr}")


    # Perform Point Packing and BIRCH algorithms
    print("Performing Point Packing and BIRCH algorithms")
    subprocess.run(["g++","-lm","-O3","SelectDNAII_PreVJ.cpp","stat.cpp"])
    tstart_pp = time.time()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures= [executor.submit(point_packing,file,size_threshold,min_distance_ratio,population_size,random_material_rate,birch_radius,fraction) for file in glob.glob("seq*.feather")]
    tend_pp = time.time()
    tspend_pp = round((tend_pp-tstart_pp),3)
    print(f"Total Point Packing and BIRCH algotithms time is {tspend_pp}")

    # Pairwise comparisons of determining hierarchical clustering
    tstart_pw = time.time()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(pairwise_process, file, cut_off, distance, link) for file in glob.glob("norm*.feather")]
    tend_pw = time.time()
    tspend_pw = round((tend_pw-tstart_pw),3)
    print(f"Total pairwise time is {tspend_pw}")


    df = pd.concat(map(pd.read_feather, glob.glob('./clu*.feather')))
    df['cluster'] = df.groupby([(df.index == 0).cumsum(), 'cluster']).ngroup().add(1)

    print("###Current total length is:",len(df))
    df.to_csv(f"AC_Pre_{subcluster}_{min_distance_ratio}_{population_size}_{random_material_rate}_{birch_radius}_{fraction}_{size_threshold}_{cut_off}_{File}.csv",index=False,sep='\t')
    for f in glob.glob("best.pcld*"):
        os.remove(f) 
    for f in glob.glob("clu*.feather"): 
        os.remove(f)
    for f in glob.glob("norm*.feather"): 
        os.remove(f)
    for f in glob.glob("seq*.feather"): 
        os.remove(f)
    for f in glob.glob("seq*.fasta"): 
        os.remove(f)

    print("Anchor clustering PreVJ pipeline has finished")
    





if __name__ == "__main__":
    start = time.perf_counter()
    main()
    end = time.perf_counter()
    print(f"Total time {round((end-start)/60,2)} mins")




















