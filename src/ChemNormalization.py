import os
import json
from datetime import datetime

import pandas as pd
import redis
from neo4j import GraphDatabase

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

    # init the output type flags
    _do_KGX: int = 0
    _do_Redis: int = 0

    # flags to turn on/off KGX and/or redis output
    def __init__(self):
        """ class constructor """
        self._config = self.get_config()
        self._debug_record_limit = self._config['debug_record_limit']
        self._do_KGX = self._config['do_KGX']
        self._do_Redis = self._config['do_Redis']

        # did we get a valid output type
        if self._do_KGX != 1 and self._do_Redis != 1:
            raise Exception('Config file settings do not specify at least one output format')

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

        self.print_debug_msg(f'Start of load.', True)

        try:
            # get the grouped and simplified SMILES
            df_gb: pd.DataFramGroupBy = self.get_simplified_smiles_for_chemicals()

            self.print_debug_msg(f'Working chemical substances...', True)

            # are we doing KGX file output
            if self._do_KGX == 1:
                # get the KGX output file handles
                out_node_f = open(self._config['KGX_output_node_file'], 'w', encoding="utf-8")
                out_edge_f = open(self._config['KGX_output_edge_file'], 'w', encoding="utf-8")
                self.print_debug_msg(f'KGX output enabled...', True)

                # write out the headers
                out_node_f.write(f'id,name,simple_smiles,category\n')
                out_edge_f.write(f'subject,edge_label,object\n')

            # are we doing redis output
            if self._do_Redis == 1:
                # get the pipelines for redis loading
                id_to_simple_smiles_pipeline = self.get_redis(self._config, 0).pipeline()
                simple_smiles_to_similar_smiles_pipeline = self.get_redis(self._config, 1).pipeline()
                self.print_debug_msg(f'Redis output enabled...', True)

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
                        if self._do_Redis == 1:
                            id_to_simple_smiles_pipeline.set(row["chem_id"], f'{simplified_SMILES}')

                        # save each element of the group to generate a list of similar SMILES
                        members.append({'id': row['chem_id'], 'ORIGINAL_SMILES': row['original_SMILES']})

                        # are we doing KGX file output
                        if self._do_KGX == 1:
                            # write out the node data to the file
                            out_node_f.write(f"{row['chem_id']},\"{row['name']}\",\"{row['original_SMILES']}\",chemical_substance\n")

                    # create an object for all the member elements
                    similar_smiles = {'members': [member for member in members], 'simplified_smiles': simplified_SMILES}

                    # are we doing KGX file output
                    if self._do_KGX == 1:
                        # write out the edges
                        for pass1 in members:
                            for pass2 in members:
                                # insure that we dont have a loopback
                                if pass1['id'] != pass2['id']:
                                    # write out the data
                                    out_edge_f.write(f"{pass1['id']},similar_to,{pass2['id']}\n")

                    # convert the data object into json format
                    final = json.dumps(similar_smiles)

                    # are we doing redis output
                    if self._do_Redis == 1:
                        self.print_debug_msg(f'Simplified SMILES to similar SMILES list -> Simplified SMILES: {simplified_SMILES}, Similar SMILES list: [{final}]\n')
                        simple_smiles_to_similar_smiles_pipeline.set(simplified_SMILES, final)
                except Exception as e:
                    self.print_debug_msg(f'Exception thrown: {e}', True)
                    continue

            # are we doing redis output
            if self._do_Redis == 1:
                self.print_debug_msg(f'Dumping to chemical substance id to lookup db ...')
                id_to_simple_smiles_pipeline.execute()

                self.print_debug_msg(f'Dumping to simple SMILES to array of similar SMILES db ...')
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

    def get_simplified_smiles_for_chemicals(self) -> pd.DataFrame:
        """ This method gets SMILES for every chemical substance in the robokop neo4j graph database and creates a simplified SMILES from each.
            The simplified SMILES values will be used as a grouping mechanism and saved in the redis database.
        """
        # Create a target data frame for the processed data
        df: pd.DataFrame = pd.DataFrame(columns=['chem_id', 'original_SMILES', 'simplified_SMILES', 'name'])

        try:
            # Create the query. This is of course robokop specific
            # DEBUG: to return 2 that have the same simplified SMILES use this in the where clause ->  and (c.id="CHEBI:140593" or c.id="CHEBI:140451")
            # DEBUG: to return 2 that have the exact same SMILES use this in the where clause ->   and (c.id="CHEBI:85764" or c.id="CHEBI:140773")
            c_query: str = f'match (c:chemical_substance) where c.smiles is not NULL and c.smiles <> "" and c.smiles <> "**" and c.smiles <> "*" and c.name is not null and c.name <> "" RETURN c.id, c.smiles, c.name order by c.smiles {self._debug_record_limit}'

            self.print_debug_msg(f"Querying target database...", True)

            # execute the query
            records: list = self.run_neo4j_query(c_query)

            self.print_debug_msg(f"Target database queried, {len(records)} chemical substance records returned.", True)

            # did we get some records
            if len(records) > 0:
                # init a counter
                rec_count: int = 0

                # loop through the records
                for r in records:
                    # increment the record counter
                    rec_count = rec_count + 1

                    # self.logger.error(f"{r['c.id']}")
                    try:
                        # Construct a molecule from a SMILES string
                        molecule: Mol = Chem.MolFromSmiles(r['c.smiles'])
                    except Exception as e:
                        # alert the user there was an issue and continue
                        self.print_debug_msg(f"Error - Exception trying to get a molecule for chem id: {r['c.id']} with original SMILES: {r['c.smiles']}, Execption {e}. Proceeding.")
                        continue

                    # did we get the molecule
                    if molecule is None:
                        # Couldn't parse the molecule
                        self.print_debug_msg(f"Error - Got an empty molecule for chem id: {r['c.id']} with smiles: {r['c.smiles']}. Proceeding.")
                        continue
                    try:
                        # get the uncharged version of the largest fragment
                        molecule_uncharged: Mol = rdMolStandardize.ChargeParent(molecule)

                        # Remove all stereo-chemistry info from the molecule
                        RemoveStereochemistry(molecule_uncharged)

                        # get the simplified SMILES value
                        simplified_smiles: str = Chem.MolToSmiles(molecule_uncharged)

                        # insure there are no dbl quotes in the name, it throws off the CSV file
                        if r['c.name'] is not None and r['c.name'] != '':
                            name_fixed = r['c.name'].replace('"', "'")
                        else:
                            name_fixed = 'No chemical name given'

                        # save the new record
                        record = {'chem_id': r['c.id'], 'original_SMILES': r['c.smiles'], 'simplified_SMILES': simplified_smiles, 'name': name_fixed}

                        # append the new record to the data frame
                        df = df.append(record, ignore_index=True)
                    except Exception as e:
                        # alert the user that something was discovered in the original graph record
                        self.print_debug_msg(f"Error - Could not get a simplified SMILES for record {rec_count}, chem id: {r['c.id']}, Original SMILES: {r['c.smiles']}, Exception: {e}")

                # get the simplified SMILES in groups
                df = df.set_index('simplified_SMILES').groupby('simplified_SMILES')
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

    def get_driver(self) -> GraphDatabase.driver:
        # Gets a connection to the graph database
        driver = GraphDatabase.driver(self._config['neo4j_uri'], auth=(self._config['neo4j_user'], self._config['neo4j_password']))

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
