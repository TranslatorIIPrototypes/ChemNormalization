import os
import json
import hashlib
import requests

from datetime import datetime

import pandas as pd
import redis
from neo4j import GraphDatabase, Driver

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdmolops import RemoveStereochemistry
import rdkit.RDLogger as Rkl


##############
# Class: ChemNormalization
#
# By: Phil Owen
# Date: 1/10/2020
# Desc: Class that gets all chemical substances from a neo4j graph database and:
#   gets a simplified SMILES value for each
#   groups the simplified SMILES and inserts them into a redis database
##############
class ChemNormalization:
    # Storage for the configuration params
    _config: json = None

    # get a reference to the rdkit logger
    logger: Rkl.logger = Rkl.logger()

    # storage for a debug limit on the number of chemical substance records processed
    _debug_record_limit: str = ''

    # flags to turn on/off KGX and/or redis output
    def __init__(self):
        """ class constructor """
        self._config = self.get_config()
        self._debug_record_limit = self._config['debug_record_limit']
        self._do_node_norm = self._config['do_node_normalization']
        self._node_norm_chunk_size = self._config['node_norm_chunk_size']
        self._node_norm_endpoint = self._config['node_norm_endpoint']
        self._do_KGX = self._config['do_KGX']
        self._do_redis = self._config['do_redis']
        self._do_curie_update = self._config['do_curie_update']

        # adjust the rdkit logging level
        self.rdkit_logging(Rkl.ERROR)

        pass

    def load(self) -> bool:
        """ Loads two redis databases with simplified chemical substance SMILES data. """
        # init the return
        rv: bool = True

        # init local variables
        id_to_simple_smiles_pipeline = None
        simple_smiles_to_similar_smiles_pipeline = None
        out_node_f = None
        out_edge_f = None

        # did we get a valid output type
        if self._do_KGX != 1 and self._do_redis != 1:
            self.print_debug_msg(f'Test/Debug mode enabled.', True)

        if self._do_KGX == 1:
            self.print_debug_msg(f'KGX output enabled.', True)

        if self._do_redis == 1:
            self.print_debug_msg(f'Redis output enabled.', True)

        if self._do_node_norm == 1:
            self.print_debug_msg(f'Node normalization enabled.', True)

        if self._do_curie_update == 1:
            self.print_debug_msg(f'CURIE prefix update enabled.', True)

        self.print_debug_msg(f'Start of load.', True)

        try:
            self.print_debug_msg(f'Collecting simplified SMILES.', True)

            # get the grouped and simplified SMILES
            df: pd.DataFrame = self.get_simplified_smiles_for_chemicals()

            # are we normalizing the chemical substance node data
            if self._do_node_norm == 1:
                self.print_debug_msg(f'Normalizing chemical substance node data.', True)

                # normalize the Node id/name data
                df = self.normalize_node_data(df)

            # get the simplified SMILES in groups
            df_gb = df.set_index('simplified_SMILES').groupby('simplified_SMILES')

            # are we doing KGX file output
            if self._do_KGX == 1:
                # get the KGX output file handles
                out_node_f = open(self._config['output_node_file'], 'w', encoding="utf-8")
                out_edge_f = open(self._config['output_edge_file'], 'w', encoding="utf-8")

                # write out the headers
                out_node_f.write(f'id,name,simple_smiles,category\n')
                out_edge_f.write(f'id,subject,edge_label,object\n')

            self.print_debug_msg(f'Putting away chemical substances.', True)

            # are we doing redis output
            if self._do_redis == 1:
                # get the pipelines for redis loading
                id_to_simple_smiles_pipeline = self.get_redis(self._config, 0).pipeline()
                simple_smiles_to_similar_smiles_pipeline = self.get_redis(self._config, 1).pipeline()

            # loop through each record and create the proper dict for redis insertion
            for simplified_SMILES, similar_SMILES_group in df_gb:
                # init the storage for the data members (ID, original SMILES, simplified SMILES
                members = []

                try:
                    # for each row in the group
                    for index, row in similar_SMILES_group.iterrows():
                        # save the chemical substance ID and simplified SMILES lookup record to the redis cache
                        self.print_debug_msg(f'ID to simplified SMILES -> Chem ID: {row["chem_id"]}, Simplified SMILES: {simplified_SMILES}, <Original SMILES: {row["original_SMILES"]}>')

                        # are we doing redis output
                        if self._do_redis == 1:
                            id_to_simple_smiles_pipeline.set(row["chem_id"], f'{simplified_SMILES}')

                        # save each element of the group to generate a list of similar SMILES
                        members.append({'id': row['chem_id'], 'ORIGINAL_SMILES': row['original_SMILES']})

                        # are we doing KGX file output
                        if self._do_KGX == 1:
                            # write out the node data to the file
                            out_node_f.write(f"{row['chem_id']},\"{row['name']}\",\"{simplified_SMILES}\",chemical_substance\n")

                    # create an object for all the member elements
                    similar_smiles = {'members': [member for member in members], 'simplified_smiles': simplified_SMILES}

                    # are we doing KGX file output
                    if self._do_KGX == 1:
                        # write out the edges
                        for pass1 in members:
                            for pass2 in members:
                                # insure that we dont have a loopback
                                if pass1['id'] != pass2['id']:
                                    # get the MD5 value for the edge ID
                                    edge_id = hashlib.md5(f"{pass1['id']}{pass2['id']}similar_to".encode('utf-8')).hexdigest()

                                    # write out the data
                                    out_edge_f.write(f"{edge_id},{pass1['id']},similar_to,{pass2['id']}\n")

                    # convert the data object into json format
                    final = json.dumps(similar_smiles)

                    # are we doing redis output
                    if self._do_redis == 1:
                        self.print_debug_msg(f'Simplified SMILES to similar SMILES list -> Simplified SMILES: {simplified_SMILES}, Similar SMILES list: [{final}]\n')
                        simple_smiles_to_similar_smiles_pipeline.set(simplified_SMILES, final)
                except Exception as e:
                    self.print_debug_msg(f'Exception thrown: {e}', True)
                    continue

            # are we doing redis output
            if self._do_redis == 1:
                self.print_debug_msg(f'Dumping to chemical substance id to lookup db.')
                id_to_simple_smiles_pipeline.execute()

                self.print_debug_msg(f'Dumping to simple SMILES to array of similar SMILES db.')
                simple_smiles_to_similar_smiles_pipeline.execute()

            # are we doing KGX file output
            if self._do_KGX == 1:
                out_node_f.close()
                out_edge_f.close()

        except Exception as e:
            self.print_debug_msg(f'Exception thrown: {e}', True)
            rv = False

        # return to the caller
        return rv

    def normalize_node_data(self, df) -> pd.DataFrame:
        """ This method calls the NodeNormalization web service to get the normalized identifier and name of the chemical substance node. """

        # init the indexs
        start_index: int = 0

        # get the last index of the list
        last_index: int = len(df)

        # declare the id to be the index
        df.set_index('chem_id')

        # grab chunks of the dataframe
        while True:
            if start_index < last_index:
                # define the end index of the slice
                end_index: int = start_index + self._node_norm_chunk_size

                # collect a group of records from the dataframe
                data_chunk = df[start_index: end_index]

                # get the data
                resp = requests.get(self._node_norm_endpoint + '?curie=' + '&curie='.join(data_chunk['chem_id'].tolist()))

                # did we get a good status code
                if resp.status_code == 200:
                    # convert to json
                    rvs = resp.json()

                    # for each row in the slice add the new id and name
                    for rv in rvs:
                        # did we find a normalized value
                        if rvs[rv] is not None:
                            # find the name and replace it with label
                            if 'label' in rvs[rv]['id']:
                                df.loc[df['chem_id'] == rv, 'name'] = rvs[rv]['id']['label']

                            # find the id and replace it
                            df.loc[df['chem_id'] == rv, 'chem_id'] = rvs[rv]['id']['identifier']
                        else:
                            self.print_debug_msg(f'{rv} has no normalized value.', False)

                # move on down the list
                start_index += self._node_norm_chunk_size
            else:
                break

        # return to the caller
        return df

    def get_simplified_smiles_for_chemicals(self) -> pd.DataFrame:
        """ This method gets SMILES for every chemical substance in the robokop neo4j graph database and creates a simplified SMILES from each.
            The simplified SMILES values will be used as a grouping mechanism and saved in the redis database.
        """
        # Create a target data frame for the processed data
        df: pd.DataFrame = pd.DataFrame(columns=['chem_id', 'original_SMILES', 'simplified_SMILES', 'name'])

        try:
            # Create the query. This is of course robokop specific
            # Query modified to exclude all chemical substances that have wildcard definitions
            c_query: str = f'match (c:chemical_substance) where c.smiles is not NULL and c.smiles <> "" and NOT c.smiles CONTAINS "*" RETURN c.id, c.smiles, c.name order by c.smiles {self._debug_record_limit}'

            self.print_debug_msg(f"Querying target database for chemical substances.", True)

            # check to see if we are in test mode
            if self._do_KGX != 0 or self._do_redis != 0:
                # execute the query
                records: list = self.run_neo4j_query(c_query)

                # to create a test data file
                # d = pd.DataFrame(records, columns=['c.id', 'c.smiles', 'c.name'])
                # d.to_json('datafile.json.test', orient='records')
            else:
                # open the test data file and use that instead of the database
                with open('./tests/datafile.json') as json_file:
                    records = json.load(json_file)

            self.print_debug_msg(f"Target database queried, {len(records)} chemical substance records returned.", True)

            # did we get some records
            if len(records) > 0:
                self.print_debug_msg(f"Processing SMILES.", True)

                # init a counter
                rec_count: int = 0

                # loop through the records
                for r in records:
                    # increment the record counter
                    rec_count = rec_count + 1

                    try:
                        # Construct a molecule from a SMILES string
                        molecule: Mol = Chem.MolFromSmiles(r['c.smiles'])
                    except Exception as e:
                        # alert the user there was an issue and continue
                        self.print_debug_msg(f"Error - Exception trying to get a molecule for record {rec_count}, chem id: {r['c.id']} with original SMILES: {r['c.smiles']}, Exception {e}. Proceeding.", True)
                        continue

                    # did we get the molecule
                    if molecule is None:
                        # Couldn't parse the molecule
                        self.print_debug_msg(f"Error - Got an empty molecule for record {rec_count}, chem id: {r['c.id']} with smiles: {r['c.smiles']}. Proceeding.", True)
                        continue
                    try:
                        # get the uncharged version of the largest fragment
                        molecule_uncharged: Mol = rdMolStandardize.ChargeParent(molecule)

                        # Remove all stereo-chemistry info from the molecule
                        RemoveStereochemistry(molecule_uncharged)

                        # get the simplified SMILES value
                        simplified_smiles: str = Chem.MolToSmiles(molecule_uncharged)

                        # convert the curie prefix to the new standard
                        if self._do_curie_update == 1:
                            chem_id = r['c.id'].replace("KEGG:", "KEGG.COMPOUND:").replace("CHEMBL:", "CHEMBL.COMPOUND:")
                        else:
                            chem_id = r['c.id']

                        # check to see if there is a name
                        if r['c.name'] is None or r['c.name'] == '' or r['c.name'] == 'NULL':
                            name_fixed = chem_id
                        else:
                            # insure there are no dbl quotes in the name, it throws off the CSV file
                            name_fixed = r['c.name'].replace('"', "'")

                        # save the new record
                        record = {'chem_id': chem_id, 'original_SMILES': r['c.smiles'], 'simplified_SMILES': simplified_smiles, 'name': name_fixed}

                        # append the new record to the data frame
                        df = df.append(record, ignore_index=True)
                    except Exception as e:
                        # alert the user that something was discovered in the original graph record
                        self.print_debug_msg(f"Error - Could not get a simplified SMILES for record {rec_count}, chem id: {r['c.id']}, Original SMILES: {r['c.smiles']}, Exception: {e}")
            else:
                self.print_debug_msg(f"No records to process.", True)

        except Exception as e:
            raise e

        # return to the caller
        return df

    @staticmethod
    def get_config() -> json:
        """ Loads the configuration settings into a class object. """

        # get the location of the configuration file
        cname = os.path.join(os.path.dirname(__file__), '..', 'config.json')

        # open the config file
        with open(cname, 'r') as json_file:
            # parse the contents of the configuration settings file
            data = json.load(json_file)

        # return to the caller
        return data

    def print_debug_msg(self, msg: str, force: bool = False):
        """ Prints a debug message if enabled in the config file """
        if self._config['debug_messages'] == 1 or force:
            now = datetime.now()

            print(f'{now.strftime("%Y/%m/%d %H:%M:%S")} - {msg}')

    @staticmethod
    def get_redis(config: json, db_id: int):
        """ Returns a connection to a redis instance """
        return redis.StrictRedis(host=config['redis_host'],
                                 port=int(config['redis_port']),
                                 db=db_id,
                                 password=config['redis_password'])

    def rdkit_logging(self, level: Rkl):
        """ Adjusts RDKit logging level. """

        # set the new logging level
        self.logger.setLevel(level)

    def get_driver(self) -> Driver:
        # Gets a connection to the graph database
        driver: Driver = GraphDatabase.driver(self._config['neo4j_uri'], auth=(self._config['neo4j_user'], self._config['neo4j_password']))

        # return to the caller
        return driver

    def run_neo4j_query(self, cypherquery: str) -> list:
        """ Runs a cypher query """

        # get a connection to the graph db
        driver = self.get_driver()

        # start a session
        with driver.session() as session:
            # get the results of the query
            results = session.run(cypherquery)

        # return the data to the caller
        return list(results)
