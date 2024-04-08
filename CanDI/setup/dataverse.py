import os
import requests
from tqdm import tqdm
import sys

"""This file contains all metadata of datasets deposited for CanDI

The datasets are deposited in the Harvard Dataverse repository:
- https://doi.org/10.7910/DVN/JIAT0H

Codes are adapted from:
- https://github.com/mims-harvard/TDC/blob/main/tdc/utils/load.py OR
- https://github.com/snap-stanford/GEARS/blob/master/gears/utils.py
"""

CANDI_DATAVERSE_DOI = 'doi:10.7910/DVN/JIAT0H'

depmap_dataset_names = [
    'CCLE_expression',
    'CCLE_fusions',
    'CCLE_gene_cn',
    'CCLE_mutations',
    'CCLE_RNAseq_reads',
    'CRISPR_gene_dependency',
    'CRISPR_gene_effect',
    'sample_info',
    'README',
]

name2type = {
    'CCLE_expression': 'csv',
    'CCLE_fusions': 'csv',
    'CCLE_gene_cn': 'csv',
    'CCLE_mutations': 'csv',
    'CCLE_RNAseq_reads': 'csv',
    'CRISPR_gene_dependency': 'csv',
    'CRISPR_gene_effect': 'csv',
    'sample_info': 'csv',
    'README': 'txt',
}

name2id = {
    'CCLE_expression': 8076862,
    'CCLE_fusions': 10085763,
    'CCLE_gene_cn': 8076861,
    'CCLE_mutations': 8076857,
    'CCLE_RNAseq_reads': 8076859,
    'CRISPR_gene_dependency': 8076863,
    'CRISPR_gene_effect': 8076860,
    'sample_info': 10085764,
    'README': 8151459,
}


def print_sys(s):
    """system print

    Args:
        s (str): the string to print
    """
    print(s, flush = True, file = sys.stderr)


def dataverse_download(url, path, name, types):
    """dataverse download helper with progress bar

    Args:
        url (str): the url of the dataset
        path (str): the path to save the dataset
        name (str): the dataset name
        types (dict): a dictionary mapping from the dataset name to the file format
    """
    save_path = os.path.join(path, f"{name}.{types[name]}")
    response = requests.get(url, stream=True)
    total_size_in_bytes = int(response.headers.get("content-length", 0))
    block_size = 1024
    progress_bar = tqdm(total=total_size_in_bytes, unit="iB", unit_scale=True)
    with open(save_path, "wb") as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()


def download_wrapper(name, path, return_type=None):
    """wrapper for downloading a dataset given the name and path, for csv,pkl,tsv or similar files

    Args:
        name (str): the rough dataset query name
        path (str): the path to save the dataset
        return_type (str, optional): the return type. Defaults to None. Can be "url", "name", or ["url", "name"]

    Returns:
        str: the exact dataset query name
    """
    server_path = "https://dataverse.harvard.edu/api/access/datafile/"

    url = server_path + str(name2id[name])

    if not os.path.exists(path):
        os.mkdir(path)

    file_name = f"{name}.{name2type[name]}"

    if os.path.exists(os.path.join(path, file_name)):
        print_sys("Found local copy...")
        os.path.join(path, file_name)
    else:
        print_sys("Downloading...")
        dataverse_download(url, path, name, name2type)
    
    if return_type == "url":
        return url
    elif return_type == "name":
        return file_name
    elif return_type == ["url", "name"]:
        return url, file_name


def depmap_dataverse_download(path, return_type=None):
    """download all datasets to the path

    Args:
        path (str): the path to save the datasets
        return_type (str, optional): the return type. Defaults to None. Can be "url", "name", or ["url", "name"]
    """
    url_list = []
    file_names = []

    for name in depmap_dataset_names:
        url, file_name = download_wrapper(name, path, return_type=["url", "name"])
        url_list.append(url)
        file_names.append(file_name)
    
    if return_type == "url":
        return url_list
    elif return_type == "name":
        return file_names
    elif return_type == ["url", "name"]:
        return url_list, file_names


# from pyDataverse.api import NativeApi, DataAccessApi
# from pyDataverse.models import Dataverse

# BASE_URL = 'https://dataverse.harvard.edu'
# api = NativeApi(BASE_URL)
# data_api = DataAccessApi(BASE_URL)
# dataset = api.get_dataset(doi)
# files_list = dataset.json()['data']['latestVersion']['files']
# for file in files_list:
#     filename = file["dataFile"]["filename"]
#     file_id = file["dataFile"]["id"]
#     print(filename, file_id)
