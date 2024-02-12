#!/usr/bin/env python
# coding: utf-8

# # Gen3/Terra Data Table Utility Functions for [BioData Catalyst](https://biodatacatalyst.nhlbi.nih.gov/) <a class="tocSkip">

# **Version:** 11/20/2020  
# This Notebook is occasionally updated with enhancements or fixes.  
# The latest version of this Notebook is available in the Terra [BioData Catalyst Collection](https://terra.biodatacatalyst.nhlbi.nih.gov/#workspaces/biodata-catalyst/BioData%20Catalyst%20Collection) Workspace [here](https://terra.biodatacatalyst.nhlbi.nih.gov/#workspaces/biodata-catalyst/BioData%20Catalyst%20Collection/notebooks/launch/terra_data_table_util.ipynb).

# # Purpose
# 
# This Notebook was developed for users of [BioData Catalyst](https://biodatacatalyst.nhlbi.nih.gov/) data on Terra. It was designed to work with the current (November 2020) Gen3 BioData Catalyst data model and the harmonzied clinical data provided for many BioData Catalyst studies.
# 
# The primary purpose of this Notebook is combining multiple Gen3 graph-structured data tables in Terra to create a single consolidated table that is easier to use. The content of the combined table produced is configurable.
# 
# The default behavior is to produce a table keyed by subject id, with one row per subject, for all subjects in a Terra Workspace. This table may include the genomic data, harmonized clinical metadata, or both, along with the associated administrative information.
# 
# This has been successfully tested with the largest single BDCat project, `parent-WHI_HMB-IRB_` with 117,675 subjects.
# 
# Additionally, convenience functions are provided for working with Terra data tables, including uploading, downloading, modification and deletion.

# # Requirements and Assumptions
# 
# **Run in Terra**  
# This Notebook is intended to be used within the Terra Notebook environment using a Python 3 Jupyter kernel. 
# 
# **Workspace Data**   
# The consolidation is performed for all Gen3 data for the BioData Catalyst program in a Terra workspace. The data may be for subjects from one or more projects/cohorts.
# 
# **Libraries**   
# The following libraries are expected to be available in the Terra Notebook environment, either by being preinstalled the 'Terra Notebook Runtime, Container Image, or explicit installation by the user:
# * `fiss` (version 0.16.23 or later)
# * `numpy` (version 1.15.2 or later)
# * `pandas` (0.25.3 or laster)
# * `tenacity` (6.1.0 or later)

# # How to Use
# There are two ways this Notebook may be used.
# 
# ## Simple Single Notebook Use
# This provides a simple/easy way to get familiar with this Notebook, and may be all that is needed for some users.  
# To use this approach:
# 
# **First set the following variable to `True`, then continue to the "Simple Single Notebook Use" instructions at the end of this Notebook.**

# In[2]:


single_notebook_use = False


# ## Multiple Notebook Use
# The functions in this Notebook may instead be called from other Notebooks.
# In this case, the functions of this Notebook are "imported" into the other Notebook using the `%run` command. 
# For example, adding the the following steps added to antoher Notebook is all that is needed to create a consolidated data table. In this example, a consolidated table containing both the genomic data and harmonized metadata.  
# To use this approach:
# 
# **Step 1. Ensure the `single_notebook_use` variable above is set to `False`.**
# 
# **Step 2.** Add the following code to the *other* Notebook:
# ```
# %run terra_data_table_util.ipynb  
# 
# PROJECT = os.environ['GOOGLE_PROJECT']  
# WORKSPACE = os.environ['WORKSPACE_NAME']
# 
# consolidate_gen3_geno_pheno_tables(PROJECT, WORKSPACE, "consolidated_metadata")
# ```
# 
# **Step 3.** Run the *other* Notebook and observe the output showing the creation of the consolidated data table.
# 
# **Step 4.** In the Terra UI "Data" tab, examine the new data table (named `"consolidated_metadata"`) that was created.
# 

# # How it Works
# 
# This Notebook uses the Broad FireCloud API to read each Terra data table identified in the merge specification into a Pandas DataFrame and performs SQL-style joins on the tables using the Pandas `merge` operation to produce a single, consolidated table.
# References in the Gen3 data model are only the direction of the graph leaf/bottom nodes to the root/top node.
# 
# During this consolidation process, the name of each column in a table is prefixed with the name of the table it is from. Additionally, the columns containing entity ids have the `_entity_id` suffix appended.

# ## Functions Provided
# 
# This Notebook provides the following pre-defined functions for creating consolidated tables:
# * `consolidate_gen3_geno_pheno_tables(project: str, workspace: str, new_table_name: str) -> None`
# * `consolidate_gen3_geno_tables(project: str, workspace: str, new_table_name: str) -> None`
# * `consolidate_gen3_pheno_tables(project: str, workspace: str, new_table_name: str) -> None`
# 
# These functions is available for processing custom merge specifications:
# * `consolidate_to_terra_table(project: str, workspace: str, merge_spec: dict, table_name: str) -> pd.DataFrame`
# * `consolidate_to_tsv(project: str, workspace: str, merge_spec: dict) -> str`
# * `consolidate_to_df(project: str, workspace: str, merge_spec: dict) -> pd.DataFrame`
# 
# The Terra data tables that are included in the consolidated table, and how they are combined, is defined by a merge specification defined as a Python dictionary.
# The merge specification supports standard SQL-style join operations and can be customized as desired.
# 
# The following functions explicitly merge individual tables.  
# These may be useful, for example, for merging user-provided data with a consolidated table produced by the functions above.
# * `merge_terra_tables_to_table(project: str, workspace: str, left_table_name: str, left_table_previously_consolidated: bool, right_table_name: str, right_table_previously_consolidated: bool, join_type: str, join_column: str, final_index_source_column: str, result_table_name: str, **kwargs) -> None`
# * `merge_terra_tables_to_df(project: str, workspace: str, left_table_name: str, left_table_previously_consolidated: bool, right_table_name: str, right_table_previously_consolidated: bool, join_type: str, join_column: str, **kwargs) -> pd.DataFrame`
# * `merge_df_and_terra_table_to_df(project: str, workspace: str, left_df: pd.DataFrame, right_table_name: str, right_table_previously_consolidated: bool, join_type: str, join_column: str, **kwargs) -> pd.DataFrame`
# 
# These convenience functions are also available:
# * `set_entity_attribute_value(project: str, workspace: str, entity_type: str, entity_name: str, attribute_name: str, value: Union[list, bool, int, float, str, None]) -> None:`
# 
# * `delete_terra_table(project: str, workspace: str, table_name: str) -> None`
# * `delete_all_gen3_tables(project: str, workspace: str) -> None`
# * `get_terra_table_to_df(project: str, workspace: str, table_name: str, attributeNames=None, model="flexible") -> pd.DataFrame:`
# * `rename_df_column(df: pd.DataFrame, current_column_name: str, new_column_name: str) -> None`
# * `write_df_to_tsv_file(df: pd.DataFrame, filename: str) -> None:`
# * `upload_entities_df(project: str, workspace: str, df: pd.DataFrame, chunk_size=500) -> None`

# # Dependencies and Imports

# Python packages only need to be installed once per Terra Cloud Environment.  
# The Jupyter kernel may need to be restarted after installing a new package for the first time.  
# Check for the packages that are required and install any that are missing.
#   
# Ensure that a recent version of firecloud is installed.
# The version must be 0.16.23 or later for flexible entity model support.

# In[3]:


import sys
import pkg_resources

installed_packages = {pkg.key for pkg in pkg_resources.working_set}

if "tenancity" not in installed_packages:
    get_ipython().run_line_magic('pip', 'install tenacity')

if "firecloud" not in installed_packages:
    get_ipython().run_line_magic('pip', 'install --upgrade firecloud')

# For debugging:
# ! pip install pysnooper


# In[4]:


import csv
from datetime import datetime
import io
import json
import math
import os
import re
import resource
import sys
from typing import Union

from firecloud import fiss
from firecloud.errors import FireCloudServerError
import firecloud.api as fapi
import numpy as np
import pandas as pd
import tenacity
import time
from tenacity import retry, after_log, before_sleep_log, retry_if_exception_type, stop_after_attempt, wait_exponential


# In[5]:


import logging
from logging import INFO, DEBUG
logger = logging.getLogger()
logger.setLevel(INFO)


# # Commonly Used Merge Specifications and Convenience Functions

# ## Create a consolidated data table containing both genomic and phenotypic data

# In[6]:


GEN3_GENO_PHENO_MERGE_SPEC = {
  "default_join_type": "outer",
  "merge_sequence": [
    {
      "join_column": "simple_germline_variation",
      "table_names": ["simple_germline_variation", "germline_variation_index"]
    },
    {
      "join_column": "submitted_aligned_reads",
      "table_names": ["submitted_aligned_reads", "aligned_reads_index"]
    },
    {
      "join_column": "read_group",
      "table_names": ["read_group", "submitted_unaligned_reads", "read_group_qc"]
    },
    {
      "join_column": "aliquot",
      "table_names": ["aliquot", "submitted_cnv_array", "submitted_snp_array"]
    },
    {
      "join_column": "sample",
      "table_names": ["sample"]
    },
    {
      "join_column": "subject",
      "table_names": ["subject", "blood_pressure_test", "cardiac_mri", "demographic", "electrocardiogram_test", "exposure", "lab_result", "medical_history", "medication"]
    },
    {
      "join_column": "study",
      "table_names": ["study"]
    },
    {
      "join_column": "project",
      "table_names": ["project"]
    },
    {
      "join_column": "program",
      "table_names": ["program"]
    }
  ],
  "final_index_source_column": "subject_submitter_id"
}


# In[7]:


def consolidate_gen3_geno_pheno_tables(project: str, workspace: str, new_table_name: str, **kwargs) -> None:
    consolidate_to_terra_table(project, workspace, GEN3_GENO_PHENO_MERGE_SPEC, new_table_name, **kwargs)


# ## Create a consolidated data table containing only genomic (no phenotypic) data

# In[8]:


GEN3_GENO_MERGE_SPEC =  {
  "default_join_type": "outer",
  "merge_sequence": [
    {
      "join_column": "simple_germline_variation",
      "table_names": ["simple_germline_variation", "germline_variation_index"]
    },
    {
      "join_column": "submitted_aligned_reads",
      "table_names": ["submitted_aligned_reads", "aligned_reads_index"]
    },
    {
      "join_column": "read_group",
      "table_names": ["read_group", "submitted_unaligned_reads", "read_group_qc"]
    },
    {
      "join_column": "aliquot",
      "table_names": ["aliquot", "submitted_cnv_array", "submitted_snp_array"]
    },
    {
      "join_column": "sample",
      "table_names": ["sample"]
    },
    {
      "join_column": "subject",
      "table_names": ["subject"]
    },
    {
      "join_column": "study",
      "table_names": ["study"]
    },
    {
      "join_column": "project",
      "table_names": ["project"]
    },
    {
      "join_column": "program",
      "table_names": ["program"]
    }
  ],
  "final_index_source_column": "subject_submitter_id"
}


# In[9]:


def consolidate_gen3_geno_tables(project: str, workspace: str, new_table_name: str, **kwargs) -> None:
    consolidate_to_terra_table(project, workspace, GEN3_GENO_MERGE_SPEC, new_table_name, **kwargs)


# ## Create a consolidated data table containing only phenotypic (no genomic) data

# Note: Here the "sample" table is being included in the phenotypic data because it contains useful identifier information (e.g., the "NWD" identifier).

# In[10]:


GEN3_PHENO_MERGE_SPEC =  {
  "default_join_type": "outer",
  "merge_sequence": [
    {
      "join_column": "subject",
      "table_names": ["subject", "sample", "blood_pressure_test", "cardiac_mri", "demographic", "electrocardiogram_test", "exposure", "lab_result", "medical_history", "medication"]
    },
    {
      "join_column": "study",
      "table_names": ["study"]
    },
    {
      "join_column": "project",
      "table_names": ["project"]
    },
    {
      "join_column": "program",
      "table_names": ["program"]
    }
  ],
  "final_index_source_column": "subject_submitter_id"
}


# In[11]:


def consolidate_gen3_pheno_tables(project: str, workspace: str, new_table_name: str, **kwargs) -> None:
    consolidate_to_terra_table(project, workspace, GEN3_PHENO_MERGE_SPEC, new_table_name, **kwargs)


# # Custom Merge Specification and Use

# In[12]:


GEN3_USER_CUSTOM_MERGE_SPEC =  {
  "default_join_type": "inner",
  "merge_sequence": [
    {
      "join_column": "simple_germline_variation",
      "table_names": ["simple_germline_variation", "germline_variation_index"]
    },
    {
      "join_column": "submitted_aligned_reads",
      "table_names": ["submitted_aligned_reads", "aligned_reads_index"]
    },
    {
      "join_column": "read_group",
      "table_names": ["read_group", "submitted_unaligned_reads", "read_group_qc"]
    },
    {
      "join_column": "aliquot",
      "table_names": ["aliquot", "submitted_cnv_array", "submitted_snp_array"]
    },
    {
      "join_column": "sample",
      "table_names": ["sample"]
    },
    {
      "join_column": "subject",
      "join_type": "left",
      "table_names": ["subject", "blood_pressure_test", "cardiac_mri", "demographic", "electrocardiogram_test", "exposure", "lab_result", "medical_history", "medication"]
    },
    {
      "join_column": "study",
      "table_names": ["study"]
    },
    {
      "join_column": "project",
      "table_names": ["project"]
    },
    {
      "join_column": "program",
      "table_names": ["program"]
    }
  ],
  "final_index_source_column": "subject_submitter_id"
}


# In[13]:


def consolidate_gen3_custom_tables(project: str, workspace: str, new_table_name: str, **kwargs) -> None:
    consolidate_to_terra_table(project, workspace, GEN3_USER_CUSTOM_MERGE_SPEC, new_table_name, **kwargs)


# # Related Convenience Functions

# Perform the merges defined in the merge specification and upload the resulting table to Terra the given name.

# In[14]:


def consolidate_to_terra_table(project: str, workspace: str, merge_spec: dict, table_name: str, **kwargs) -> None:
    
    # dest_project and dest_workspace can be used to write the consolidated data table
    # to a different billing project and workspace than the one contining the Gen3 data tables.
    # This can be useful for working around some of the current Terra scale limitations, 
    # and for other researcher use cases.
    
    dest_project = kwargs.get("dest_project", project)
    dest_workspace = kwargs.get("dest_workspace", workspace)
    
    if 'final_index_source_column' in merge_spec and len(merge_spec['final_index_source_column']):
        entity_id_column = merge_spec['final_index_source_column']
    else:
        logger.error("The merge specification field \"final_index_source_column\" is missing or has an empty value.")
        return
    
    logger.info("Starting data consolidation to table \"{}\".".format(table_name)); _flush_log();
    
    _check_and_log_existing_table(dest_project, dest_workspace, table_name)
    
    consolidated_df = consolidate_to_df(project, workspace, merge_spec)
 
    # Add "entity:{table_name}_id" as the first column, as required by Terra.
    consolidated_df.insert(0, f"entity:{table_name}_id", consolidated_df[entity_id_column])
    
    upload_entities_df_and_verify(dest_project, dest_workspace, consolidated_df) 
        
    logger.info(_get_python_resource_usage())


# Perform the merges defined in the merge specification and return the resulting table in a TSV format string.

# In[15]:


def consolidate_to_tsv(project: str, workspace: str, merge_spec: dict)  -> str:
    return consolidate_to_df(project, workspace, merge_spec).to_csv(sep="\t")


# Perform the merges defined in the merge specification and return the resulting table as a Pandas DataFrame.

# In[16]:


# @pysnooper.snoop
def consolidate_to_df(project: str, workspace: str, merge_spec: dict)  -> pd.DataFrame:
    default_merge_parameters = merge_spec['default_merge_parameters'] if 'default_merge_parameters' in merge_spec else dict(how="outer")
    if "default_join_type" in merge_spec:
        default_merge_parameters['how'] = merge_spec['default_join_type']
        
    merged_df = None
    for merge_info in merge_spec['merge_sequence']:
        merge_parameters = _create_combined_merge_parameters(default_merge_parameters, merge_info)
        _substitute_entity_id_column_name(merge_parameters)
        merged_df = _consolidate_tables_to_df(project, workspace, merge_info['table_names'], merge_parameters, merged_df)
    return merged_df


# Merge two Terra tables to a Terra table using the specified merge parameters

# In[17]:


def merge_terra_tables_to_table(project: str, workspace: str,
                                left_table_name: str, left_table_previously_consolidated: bool,
                                right_table_name: str, right_table_previously_consolidated: bool,
                                join_type: str, join_column: str,
                                final_index_source_column: str, result_table_name: str,
                                **kwargs) -> None:
    
    dest_project = kwargs.get("dest_project", project)
    dest_workspace = kwargs.get("dest_workspace", workspace)
    
    _check_and_log_existing_table(dest_project, dest_workspace, result_table_name)
    
    merged_df = merge_terra_tables_to_df(project, workspace,
                                left_table_name, left_table_previously_consolidated,
                                right_table_name, right_table_previously_consolidated,
                                join_type, join_column, **kwargs)
    
    # Add "entity:{entity_name}_id" as the first column, as required by Terra.
    merged_df.insert(0, f"entity:{result_table_name}_id", merged_df[final_index_source_column])
    
    upload_entities_df_and_verify(dest_project, dest_workspace, merged_df)


# Merge two Terra tables to a Pandas DataFrame using the specified merge parameters

# In[18]:


def merge_terra_tables_to_df(project: str, workspace: str,
                            left_table_name: str, left_table_previously_consolidated: bool,
                            right_table_name: str, right_table_previously_consolidated: bool,
                            join_type: str, join_column: str, **kwargs) -> pd.DataFrame: 
    if left_table_previously_consolidated:
        left_table_df = get_terra_table_to_df(project, workspace, left_table_name)
    else:
        left_df = get_gen3style_terra_table_to_df(project, workspace, left_table_name)
        
    merged_df = merge_df_and_terra_table_to_df(project, workspace,
                                               left_df,
                                               right_table_name, right_table_previously_consolidated,
                                               join_type, join_column, **kwargs)
        
    return merged_df


# Merge a DataFrame and a Terra table to a Panda DataFrame using the specified merge parameters

# In[19]:


def merge_df_and_terra_table_to_df(project: str, workspace: str,
                            left_df: pd.DataFrame,
                            right_table_name: str, right_table_previously_consolidated: bool,
                            join_type: str, join_column: str, **kwargs) -> pd.DataFrame:
        
    if right_table_previously_consolidated:
        right_table_df = get_terra_table_to_df(project, workspace, right_table_name)
    else:
        right_df = get_gen3style_terra_table_to_df(project, workspace, right_table_name)
        
    merge_parameters = dict(how=join_type, on=join_column, **kwargs)
    merged_df = left_df.merge(right_df, **merge_parameters)
    
    logger.info("Merged table \"{}\" on column \"{}\" with join type: \"{}\". New merged table dimmensions: ({}x{})".format(
    right_table_name, merge_parameters['on'], merge_parameters['how'], merged_df.shape[0], merged_df.shape[1]))
    
    if logger.isEnabledFor(DEBUG):
        logger.debug("The in-memory merged data frame size is: {} rows x {} columns".format(merged_df.shape[0], merged_df.shape[1]))
        _debug_write_df_to_tsv_file(merged_df, "merged_df")
        
    return merged_df


# Create a Pandas DataFrame containing the contents of the given Terra data table.

# In[20]:


def get_terra_table_to_df(project: str, workspace: str, table_name: str, attributeNames=None, model="flexible") -> pd.DataFrame:
    data_table_info = DataTableInfo(project, workspace)
    row_count, _, _ = data_table_info.get_table_info(table_name)
    single_read_max_size = 5000
    if row_count <= single_read_max_size:
        # Process as a single read operation
        response = _fapi_get_entities_tsv(project, workspace, table_name, attributeNames, model=model)
        table_df = pd.read_csv(io.StringIO(response.text), sep='\t')
    else:
        table_df = _get_large_terra_table_to_df(project, workspace, table_name, attributeNames)
    
    # Change the dataframe index from the default numeric index to the the entity id column.
    # TODO - Resetting the index below had the unexpected effect of causing the subsequent merge
    #        operation to fail due to a key error, even though the intended key was present
    #        in both tables. Omit the following until it can be investigated and resolved.
    # table_df.set_index(f"entity:{table_name}_id", inplace=True)
    
    return table_df


# Create a Pandas DataFrame containing the contents of the given Terra data table, with columns renamed to facilitate merging:
# * In general, column names are prefixed with the name of the table to address conflicts that would otherwise occur due to fields having the same name in multiple different tables.  
# * Columns representiong relationships between tables are suffixed with `_entity_id`.

# In[21]:


def get_gen3style_terra_table_to_df(project: str, workspace: str, table_name: str, model="flexible") -> pd.DataFrame:
    table_df = get_terra_table_to_df(project, workspace, table_name)
    columns = table_df.columns
    rename_df_column(table_df, f"entity:{table_name}_id", f"{table_name}_entity_id") # Column 0
    for column in columns[1:]:
        canonical_column_name = _get_attribute_canonical_name(column)
        if canonical_column_name in _GEN3_TABLE_NAMES:
            rename_df_column(table_df, column, f"{canonical_column_name}_entity_id")
        else:
            rename_df_column(table_df, column, f"{table_name}_{canonical_column_name}")
    # Deduplicate "*_entity_id" columns
    table_df = table_df.loc[:,~table_df.columns.duplicated()]
    return table_df


# Rename a column in the given Pandas DataFrame.

# In[22]:


def rename_df_column(df: pd.DataFrame, current_column_name: str, new_column_name: str) -> None:
    df.rename(columns={current_column_name : new_column_name}, inplace=True)


# Write the contents of the Pandas DataFrame to given filename on the file system.

# In[23]:


def write_df_to_tsv_file(df: pd.DataFrame, filename: str) -> None:
    with open(filename, mode="w") as tsv_file:
        tsv_string = df.to_csv(sep="\t", index=False)
        tsv_file.write(tsv_string)


# Upload the contents of the Pandas DataFrame to a Terra data table.  
# This includes support for "chunking" large tables into smaller sections that can be successfully uploaded individually.
# 
# Note: The format of the table within the Pandas DataFrame must comform to the format described here: https://support.terra.bio/hc/en-us/articles/360025758392-Managing-data-with-tables-

# In[24]:


def upload_entities_df(project: str, workspace: str, df: pd.DataFrame, chunk_size=500) -> None:
    logger.info("Starting upload of data table to Terra."); _flush_log()
    chunk_start = chunk_end = 0
    row_count = df.shape[0]
    first_iteration = True
    start_time = time.time()
    while chunk_start < row_count:
        chunk_end = min(chunk_start + chunk_size, row_count)
        chunk_df = df.iloc[chunk_start:chunk_end]
        chunk_tsv = chunk_df.to_csv(sep="\t", index=False)
        _fapi_upload_entities(project, workspace, chunk_tsv, "flexible")
        chunk_start = chunk_end
        if first_iteration:
            first_iteration = False
            iteration_duration = time.time() - start_time
            estimated_duration = _estimate_total_duration(row_count, chunk_size, iteration_duration)
            logger.info("Estimated time to upload table with {} rows: {}".format(row_count, estimated_duration))
            _output_now("Uploading ")
        _output_now(".")
    _output_now("\n")
    total_duration = time.time() - start_time
    logger.info("\nFinished uploading data table in {} minutes.".format(str(round((total_duration / 60), 1))))    


# In[25]:


def upload_entities_df_and_verify(project: str, workspace: str, df: pd.DataFrame, chunk_size=500) -> None:
    
    entity_id_column_name = df.columns[0]
    if not re.match("entity:.+_id", entity_id_column_name):
        raise TerraDataUtilException("The first column name does not match the pattern required to upload: {}".formt(entity_id_column_name))
    table_name = entity_id_column_name.replace("entity:", "").replace("_id", "")  
    
    df_rows, df_columns = df.shape
    if logger.isEnabledFor(DEBUG):
        logger.info("The in-memory data table size is: {} rows x {} columns".format(df.shape[0], df.shape[1]))
        _debug_write_df_to_tsv_file(df, "upload_df")
    
    upload_entities_df(project, workspace, df, chunk_size)
    
    # Compare the in-memory and actual uploaded data table sizes and output the results.
    data_table_info = DataTableInfo(project, workspace)
    actual_rows, actual_columns, _ = data_table_info.get_table_info(table_name)
    if (df_rows == actual_rows and df_columns == actual_columns):
        logger.info("The data table \"{}\" size is: {} rows x {} columns".format(
            table_name, actual_rows, actual_columns))
    else:
        if (df_rows > actual_rows or df_columns > actual_columns):
            logger.error("Data table truncation error."
                         " The in-memory data table has more rows or columns ({}x{}) than the data table \"{}\" uploaded to Terra ({}x{})".format(
                           df_rows, df_columns, table_name, actual_rows, actual_columns))
        else:
            logger.warning("Data table size mismatch warning."
                           " The in-memory data table has fewer rows or columns ({}x{}) than the data table \"{}\" uploaded to Terra ({}x{})".format(
                           df_rows, df_columns, table_name, actual_rows, actual_columns)) 


# Add or update the given attribute and value into a Terra data table or set.  
# This may be used, for example, to set an array of strings into a cell, for subsequent use as a workflow input parameter.

# In[26]:


def set_entity_attribute_value(project: str, workspace: str, entity_type: str, entity_name: str, attribute_name: str,
                               value: Union[list, bool, int, float, str, None]) -> None:
    # See: https://software.broadinstitute.org/firecloud/documentation/article?id=10892
    set_attribute_json = [fapi._attr_set(attribute_name, value)]
    _fapi_update_entity(project, workspace, entity_type, entity_name, set_attribute_json)


# Delete the Terra data table with the given billing project, workspace and name.

# In[27]:


def delete_terra_table(project: str, workspace: str, table_name: str) -> None:
    logger.info("Preparing to delete data table \"{}\" ...".format(table_name)); _flush_log()
    if table_name not in DataTableInfo(project, workspace).get_table_names():
        logger.warning("Data table \"{}\" not found.".format(table_name))
        return
    
    logger.info("Starting deletion of data table \"{}\"".format(table_name)); _flush_log()
    entity_id_series = get_table_entity_ids_to_series(project, workspace, table_name)
    chunk_size = 100
    num_chunks = math.ceil(entity_id_series.size / chunk_size)
    first_iteration = True
    start_time = time.time()
    for chunk in  np.array_split(entity_id_series, num_chunks):
        # The FireCloud API requires entity ids to be strings, not a numeric type.
        # chunk_as_strings = [str(id) for id in chunk]
        response = _fapi_delete_entity_type(project, workspace, table_name, chunk)
        if first_iteration:
            first_iteration = False
            iteration_duration = time.time() - start_time
            estimated_duration = _estimate_total_duration(entity_id_series.size, chunk_size, iteration_duration)
            logger.info("Estimated time to delete table with {} rows: {}".format(entity_id_series.size, estimated_duration))
            _output_now("Deleting ")
        _output_now(".")
    _output_now("\n")
    total_duration = time.time() - start_time
    logger.info("\nFinished deleting data table \"{}\" in {} minutes.".format(table_name, str(round((total_duration / 60), 1))))    


# Delete all Gen3 data tables in the given billing project and workspace.

# In[28]:


def delete_all_gen3_tables(project: str, workspace: str):
    logger.info("Deleting all Gen3 tables in workspace \"{}\". This may require a very long time depending on the number and size of the Gen3 tables.".format(workspace))
    data_table_info = DataTableInfo(project, workspace)
    for gen3_table_name in _GEN3_TABLE_NAMES:
        if gen3_table_name in data_table_info.get_table_names():
            delete_terra_table(project, workspace, gen3_table_name)
    logger.info("Finished deleting all Gen3 tables in workspace \"{}\".".format(workspace))


# Return the full list of entity ids for a Terra data table as a Series of strings

# In[29]:


def get_table_entity_ids_to_series(project: str, workspace: str, table_name: str) -> pd.Series:
    entity_id_column_name = f"entity:{table_name}_id"
    response = _fapi_get_entities_tsv(project, workspace, table_name, attributeNames=[entity_id_column_name], model="flexible")
    table_df = pd.read_csv(io.StringIO(response.text), sep='\t', usecols=[entity_id_column_name], dtype="str")
    entity_id_series = table_df[entity_id_column_name]
    return entity_id_series


# # Internals

# Data and functions used internally and not intended for user modification.  
# *The code in the rest of this document will likely be moved to a new Python library at some point.*

# This is the list of tables defined in the Gen3 data model for BioData Catalyst, for use internal to this Notebook.  
# All of the tables used in merge specications must exist in this list, yet this list may contain additional tables names are not used in the merge specifications and do not exist in the current workspace.
# 
# For use when deleting all Gen3 tables, this list must be a partial ordering based on the Gen3 dependencies between tables, from the leaves of the Gen3 graph data model to the root.

# In[30]:


_GEN3_TABLE_NAMES = [
    "germline_variation_index",
    "simple_germline_variation",
    "aligned_reads_index",
    "submitted_unaligned_reads",
    "submitted_aligned_reads",
    "read_group_qc",
    "read_group",
    "submitted_cnv_array",
    "submitted_snp_array",
    "aliquot",
    "sample",
    "blood_pressure_test",
    "cardiac_mri",
    "demographic",
    "electrocardiogram_test",
    "exposure",
    "lab_result",
    "medical_history",
    "medication",
    "subject",
    "study",
    "reference_file_index",
    "reference_file",
    "project",
    "program"
]


# In[31]:


if logger.isEnabledFor(DEBUG):
    get_ipython().run_line_magic('xmode', 'Verbose')
    import pysnooper


# In[32]:


# @pysnooper.snoop()
def _consolidate_tables_to_df(project: str, workspace: str, table_names: list, merge_parameters: dict, initial_df = None) -> pd.DataFrame:
    if initial_df is None:
        assert len(table_names) >= 2, "At least two table names are required." 
        table_name = table_names[0]
        table_names = table_names[1:]
        first_df = get_gen3style_terra_table_to_df(project, workspace, table_name)
        if table_name == "sample":
            _deduplicate_merge_data(None, first_df, "sample", _get_entity_id_column_name("subject"))
        merged_df = first_df
    else:
        assert len(table_names) >= 1, "At least one table names is required to merge with previous data."
        merged_df = initial_df
        
    data_table_info = DataTableInfo(project, workspace)
    for table_name in table_names:
        if table_name not in data_table_info.get_table_names():
            logger.info("The table \"{}\" was not found in this workspace and will be ignored.".format(table_name))
            continue            
        current_df = get_gen3style_terra_table_to_df(project, workspace, table_name)
        if table_name == "sample":
            _deduplicate_merge_data(merged_df, current_df, "sample", _get_entity_id_column_name("subject"))
            
        if logger.isEnabledFor(DEBUG):
            _debug_write_df_to_tsv_file(merged_df, "merged_df")
            _debug_write_df_to_tsv_file(current_df, "current_df")
            
        logger.debug("Merging table \"{}\" using column \"{}\" with join type: {}".format(
            table_name, merge_parameters['on'], merge_parameters['how']))
        logger.debug("Full merge parameters: {}".format(merge_parameters))
        
        merged_df = _merge_dataframes(merged_df, current_df, **merge_parameters)
        
        logger.info("Merged table \"{}\" on column \"{}\" with join type: \"{}\". New merged table dimmensions: ({}x{})".format(
            table_name, merge_parameters['on'], merge_parameters['how'], merged_df.shape[0], merged_df.shape[1]))
        
    return merged_df


# In[33]:


def _get_entity_id_column_name(entity_type: str):
    return f"{entity_type}_entity_id"


# In[34]:


def _check_and_log_existing_table(project: str, workspace:str , table_name:str) -> None:
    data_table_info = DataTableInfo(project, workspace)
    if (table_name in data_table_info.get_table_names()):
        existing_rows, existing_columns, _ = data_table_info.get_table_info(table_name)
        logger.info("A data table with the name \"{}\" already exists with dimmesions ({}x{}). Corresponding data will be updated and any existing additional data will be left unchanged.".format(
        table_name, existing_rows, existing_columns))


# In[35]:


def _get_large_terra_table_to_df(project: str, workspace: str, table_name: str, attributeNames=None) -> pd.DataFrame:
    total_row_count, _, _ = DataTableInfo(project, workspace).get_table_info(table_name)
    page_size = 5000
    num_pages = int(math.ceil(float(total_row_count) / page_size))
    entity_results_list = []
    for i in range(1, num_pages + 1):
        entity_results_list.append(_fapi_get_entities_query(project, workspace, table_name, i, page_size).json())

    row_jsons = []
    field_names = set()
    for results in entity_results_list:
        for result_json in results['results']:
            row_json = _format_row_json(result_json)
            field_names = field_names.union(row_json.keys())
            row_jsons.append(row_json)

    tsv_data = io.StringIO()
    try:
        field_name_list = sorted(list(field_names))
        dict_writer = csv.DictWriter(tsv_data, field_name_list, dialect=csv.excel_tab)
        dict_writer.writeheader()
        dict_writer.writerows(row_jsons)
        table_df = pd.read_csv(io.StringIO(tsv_data.getvalue()), sep='\t')

        entity_id_column_name = f"entity:{table_name}_id"
        if attributeNames is not None:
            columns = sorted(attributeNames)
        else:
            columns = list(table_df.columns)
            columns.remove(entity_id_column_name)
        columns.insert(0, entity_id_column_name)
        table_df = table_df[columns]
    finally:
        tsv_data.close()

    return table_df

def _format_row_json(result_json: dict) -> dict:
    row_json = result_json['attributes']
    for key, value in row_json.items():
        # Process a reference
        if type(value) == dict and 'entityType' in value:
            assert _get_attribute_canonical_name(key) == value['entityType']
            row_json[key] = value['entityName']
    entity_id_name = f"entity:{result_json['entityType']}_id"
    entity_id_value = result_json['name']
    row_json[entity_id_name] = entity_id_value
    return row_json


# In[36]:


def _create_combined_merge_parameters(default_merge_parameters: dict, merge_info: dict) -> dict:
    standard_pandas_default_parameters = dict(how="inner", on=None, left_on=None, right_on=None, left_index=False, right_index=False, sort=False, suffixes=("_x", "_y"), copy=True, indicator=False, validate=None)
    combined_parameters = standard_pandas_default_parameters.copy()
    combined_parameters.update(default_merge_parameters)
    if 'merge_parameters' in merge_info:
        combined_parameters.update(merge_info['merge_parameters'])
    if 'join_column' in merge_info:
        combined_parameters['on'] = merge_info['join_column']
    if 'join_type' in merge_info:
        combined_parameters['how'] = merge_info['join_type']
    return combined_parameters


# In[37]:


def _substitute_entity_id_column_name(merge_parameters: dict) -> dict:
    for key in 'on', 'left_on', 'right_on':
        if key in merge_parameters and merge_parameters[key]:
            merge_parameters[key] = _get_entity_id_column_name(merge_parameters[key])
            # TODO - Add support for case where value is a list/array - requires careful testing


# In[38]:


def _get_attribute_canonical_name(attribute_name: str) -> str:
    return attribute_name.split(":")[-1]


# In[39]:


def _deduplicate_merge_data(merged_df: pd.DataFrame, current_df: pd.DataFrame,
                           current_table_name: str, current_dedup_key: str) -> None:
    # Some TOPMed projects (COPDGene, MESA, maybe others) are known to have multiple sample
    # entries for the same subject. According to BioData Catalyst data experts,
    # the duplicates should be equivalent, so just keep the first entry found in each case.

    # Identify duplicates in the given column of the current table and obtain
    # a list of entity ids for the rows containing duplicates.
    # Then remove the duplicate rows from the current table.
    current_dups = current_df[current_dedup_key].duplicated(keep="first")
    current_dups_values = current_df[current_dups][current_dedup_key].tolist()
    if len(current_dups_values) == 0:
        logger.debug("No duplicates found in table {} for key {}".format(current_table_name, current_dedup_key))
        return
    current_table_entity_id = _get_entity_id_column_name(current_table_name)
    common_key_values_for_dupes = current_df[current_dups][current_table_entity_id].tolist()
    current_df.drop(current_df[current_dups].index, inplace=True)
    logger.warning("Removed {} duplicate entries from table \"{}\" in column \"{}\". Retained the first entry found. Deleted rows with ids: {}".format(
        len(current_dups_values), current_table_name, current_dedup_key, current_dups_values))

    # From the results that have been merged thus far, remove the rows that would have been joined
    # to the rows that were deleted as duplicates from the current table. This will prevent "orphan"
    # rows from being created in the consolidated dataframe, which would otherwise happen with
    # some join types (e.g. "outer").
    if merged_df is not None and current_table_entity_id in merged_df.columns:
        mask = merged_df[current_table_entity_id].isin(common_key_values_for_dupes)
        merged_df.drop(merged_df[mask].index, inplace=True)


# In[40]:


def _merge_dataframes(left_df: pd.DataFrame, right_df: pd.DataFrame, **merge_parameters) -> pd.DataFrame:
    merged_df = left_df.merge(right_df, **merge_parameters)
        
    # Deduplicate "*_entity_id" columns
    merged_df = merged_df.loc[:,~merged_df.columns.duplicated()]
        
    return merged_df


# In[41]:


def _debug_write_df_to_tsv_file(df: pd.DataFrame, filename: str) -> None:
    filename += "_" + datetime.now().strftime("%Y%m%d_%H%M%S%f") + ".tsv"
    write_df_to_tsv_file(df, filename)


# List the table names and attributes in the given workspace.  
#   
# Note: The FireCloud Orchestration layer accessed by FISS
# imposes a 1 minute timeout on the operation, which is not
# long enough for workspacs with many large tables. Use `_rawls_list_entity_types` instead, which has a longer, 2 minute timeout.

# In[42]:


@retry(reraise=True,
       retry=retry_if_exception_type(FireCloudServerError), 
       stop=stop_after_attempt(5),
       wait=wait_exponential(multiplier=4, min=10, max=60),
       after=after_log(logger, logging.DEBUG),
       before_sleep=before_sleep_log(logger, logging.INFO))
def _fapi_list_entity_types(project: str, workspace: str):
    response = fapi.list_entity_types(project, workspace)
    fapi._check_response_code(response, 200)
    return response


# List the table names and attributes in the given workspace.  
# 
# This produces the same result as FISS `list_entity_types`,
# yet accesses the Terra Rawls service directly to avoid the 1 minute timeout
# imposted by the FireCloud Orchestration layer. Rawls instead has a 2 minute timeout for this operation, which makes the difference between failure and success for some workspaces with very large data tables.

# In[43]:


@retry(reraise=True,
       retry=retry_if_exception_type(FireCloudServerError), 
       stop=stop_after_attempt(5),
       wait=wait_exponential(multiplier=4, min=10, max=60),
       after=after_log(logger, logging.DEBUG),
       before_sleep=before_sleep_log(logger, logging.INFO))
def _rawls_list_entity_types(namespace: str, workspace: str):
    RAWLS_PRODUCTION_HOSTNAME = "rawls.dsde-prod.broadinstitute.org"
    
    # For convenience and consistency, use the FISS session and auth infrastructure.
    if fapi.__SESSION is None:
        fapi._set_session()
        
    headers = {"Content-type":  "application/json"}
    url = f"https://{RAWLS_PRODUCTION_HOSTNAME}/api/workspaces/{namespace}/{workspace}/entities"
    response = fapi.__SESSION.get(url, headers=headers)
    fapi._check_response_code(response, 200)
    return response


# In[44]:


@retry(reraise=True,
       retry=retry_if_exception_type(FireCloudServerError), 
       stop=stop_after_attempt(5),
       wait=wait_exponential(multiplier=4, min=10, max=60),
       after=after_log(logger, logging.DEBUG),
       before_sleep=before_sleep_log(logger, logging.INFO))
def _fapi_get_entities_tsv(project: str, workspace: str, table_name: str, attributeNames=None, model="flexible"):
    response = fapi.get_entities_tsv(project, workspace, table_name, attributeNames, model=model)
    fapi._check_response_code(response, 200)
    return response


# In[45]:


@retry(reraise=True,
       retry=retry_if_exception_type(FireCloudServerError),
       stop=stop_after_attempt(5),
       wait=wait_exponential(multiplier=4, min=10, max=60),
       after=after_log(logger, logging.DEBUG),
       before_sleep=before_sleep_log(logger, logging.INFO))
def _fapi_get_entities_query(project: str, workspace: str, table_name: str,
                             page=1, page_size=100, sort_direction="asc", filter_terms=None):
    response = fapi.get_entities_query(project, workspace, table_name, page, page_size, sort_direction, filter_terms)
    fapi._check_response_code(response, 200)
    return response


# In[46]:


@retry(reraise=True,
       retry=retry_if_exception_type(FireCloudServerError), 
       stop=stop_after_attempt(5),
       wait=wait_exponential(multiplier=4, min=10, max=60),
       after=after_log(logger, logging.DEBUG),
       before_sleep=before_sleep_log(logger, logging.INFO))
def _fapi_upload_entities(project: str, workspace: str, entity_tsv: str, model: str):
    response = fapi.upload_entities(project, workspace, entity_tsv, model)
    fapi._check_response_code(response, 200)


# In[47]:


@retry(reraise=True,
       retry=retry_if_exception_type(FireCloudServerError), 
       stop=stop_after_attempt(3),
       wait=wait_exponential(multiplier=4, min=5, max=30),
       after=after_log(logger, logging.DEBUG),
       before_sleep=before_sleep_log(logger, logging.INFO))
def _fapi_update_entity(project: str, workspace: str, entity_type: str, entity_name: str, updates: list):
    # See: https://software.broadinstitute.org/firecloud/documentation/article?id=10892
    response = fapi.update_entity(project, workspace, entity_type, entity_name, updates)
    fapi._check_response_code(response, 200)


# In[48]:


@retry(reraise=True,
       retry=retry_if_exception_type(FireCloudServerError), 
       stop=stop_after_attempt(3),
       wait=wait_exponential(multiplier=4, min=5, max=30),
       after=after_log(logger, logging.DEBUG),
       before_sleep=before_sleep_log(logger, logging.INFO))
def _fapi_delete_entity_type(namespace: str, workspace: str, etype: str, ename) -> dict:
    response = fapi.delete_entity_type(namespace, workspace, etype, ename)
    if response.status_code == 409:
        message = f"Please remove existing references to entities in {etype} and try again. References: {response.content}"
        logger.warning(message)
        raise TerraDataUtilException(message)
    fapi._check_response_code(response, 204)
    return response


# In[49]:


def _get_python_resource_usage() -> str:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    memory_use_mb = math.ceil(usage.ru_maxrss / 1024)
    display_string = f"Python Memory Use: {memory_use_mb} mb"
    return display_string


# In[50]:


def _estimate_total_duration(total_count: int, batch_size: int, batch_duration_seconds: float) -> str:
    batches = math.ceil(total_count / batch_size)
    total_duration_seconds = batches * batch_duration_seconds
    return f"{math.ceil(total_duration_seconds/60)} minutes"


# In[51]:


def _flush_log() -> None:
    sys.stderr.flush()


# In[52]:


def _output_now(message: str) -> None:
    sys.stderr.write(message)
    sys.stderr.flush()


# In[53]:


class DataTableInfo:
    
    def __init__(self, project: str, workspace: str):
        self.project = project
        self.workspace = workspace
        self._data_table_info = None
        self._data_table_names = None

    def refresh(self):
        response = _rawls_list_entity_types(self.project, self.workspace)
        assert response.status_code == 200
        self._data_table_info = json.loads(response.text)
        self._data_table_names = list(self._data_table_info.keys())

    def get_data_table_info(self, refresh=False):
        if not self._data_table_info or refresh:
            self.refresh()
        return self._data_table_info.copy()

    def get_table_names(self, refresh=False):
        if not self._data_table_names or refresh:
            self.refresh()
        return self._data_table_names.copy()

    def get_table_info(self, table_name, refresh=False):
        if not self._data_table_info or refresh:
            self.refresh()
        row_count = None
        column_count = None
        attributes = None
        if table_name in self._data_table_names:
            row_count = self._data_table_info[table_name]['count']
            attributes = self._data_table_info[table_name]['attributeNames'].copy()
            column_count = len(attributes) + 1  # Add one for the entity id column
        return row_count, column_count, attributes


# Patch the `FireCloudServerError` exception string to include the HTTP response code.

# In[54]:


def _code_and_message(self):
    return f"code: {self.code} {self.message}"
FireCloudServerError.__str__ = _code_and_message


# In[55]:


class TerraDataUtilException(Exception):
    def __init__(self, message):
        self.message = message
        Exception.__init__(self, message)


# ## Simple Single Notebook Use
# 
# This section describes running this Notebook directly to consolidate all the subject-level data in this same workspace into a single table. The code below will include both the genomic and phenotypic data in a new table named "consolidated metadata". Both the data that is included and the name of the resulting table may be customized.
# 
# **Step 1.** Ensure the `single_notebook_use` variable is set to `True` (above, near top of this Notebook)
# 
# **Step 2.** Run this entire Notebook, for example using the menu item:
#   "Kernel": "Restart & Run All"
#   
# **Step 3.** Observe the output of the cell below to see what the program is doing.
# 
# **Step 4.** In the Terra UI "Data" tab, examine the new data table (named "consolidated_metadata") that was created.

# In[64]:


if single_notebook_use:
    PROJECT = os.environ['GOOGLE_PROJECT']  
    WORKSPACE = os.environ['WORKSPACE_NAME']
    consolidate_gen3_geno_pheno_tables(PROJECT, WORKSPACE, "consolidated_metadata")


# (Optional) If desired, delete the data table created above.  
# Uncommnet the following lines, then run the cell below.

# In[57]:


# if single_notebook_use:
#     PROJECT = os.environ['GOOGLE_PROJECT']  
#     WORKSPACE = os.environ['WORKSPACE_NAME']
#     delete_terra_table(PROJECT, WORKSPACE, "consolidated_metadata")


# **End of Notebook**
