import os
import json
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
#   gets a simplified SMILEs value for each
#   groups the simplified SMILES and inserts them into a redis database
##############
class ChemNormalization:
    # Storage for the configuration params
    _config: json = None

    # storage for a debug limit on the number of chemical substance records processed
    _debug_record_limit: str = ''

    def __init__(self):
        """ class constructor """
        self._config = self.get_config()
        self._debug_record_limit = self._config['debug_record_limit']

        # adjust the rdkit logging level
        self.rdkit_logging(Rkl.ERROR)

        pass

    def load_redis(self) -> bool:
        """ Loads two redis databases with simplified chemical substance SMILES data. """
        # init the return
        rv: bool = True

        try:
            # get a connection to the redis instance for the (ID -> Simple SMILES) and (Simple SMILES -> [Similar SMILES, ..]) data
            #id_to_simple_smiles_redis = self.get_redis(self._config, 0)
            #simple_smiles_to_similar_smiles_redis = self.get_redis(self._config, 1)

            # get the pipelines for redis loading
            #id_to_simple_smiles_pipeline = id_to_simple_smiles_redis.pipeline()
            #simple_smiles_to_similar_smiles_pipeline = simple_smiles_to_similar_smiles_redis.pipeline()

            # get the grouped and simplified SMILES
            df_gb: pd.DataFramGroupBy = self.get_simplified_smiles_for_chemicals()

            # loop through each record and create the proper dict for redis insertion
            for simplified_SMILES, similar_SMILES_group in df_gb:
                # init the storage for the data members (ID, original SMILES, simplified SMILES
                members = []

                # for each row in the group
                for index, row in similar_SMILES_group.iterrows():
                    # save the chemical substance ID and simplified SMILES lookup record to the redis cache
                    #self.print_debug_msg(f'ID to simplified SMILES -> Chem ID: {row["chem_id"]}, Simplified SMILES: {simplified_SMILES}, <Original SMILES: {row["original_SMILES"]}>')
                    #id_to_simple_smiles_pipeline.set(row["chem_id"], f'{simplified_SMILES}')

                    # save each element of the group to generate a list of similar SMILES
                    members.append({'id': row['chem_id'], 'ORIGINAL_SMILES': row['original_SMILES']})

                # create an object for all the member elements
                similar_smiles = {'members': [member for member in members], 'simplified_smiles': simplified_SMILES}

                # convert the data object into json format
                final = json.dumps(similar_smiles)

                #self.print_debug_msg(f'Simplified SMILES to similar SMILES list -> Simplified SMILES: {simplified_SMILES}, Similar SMILES list: [{final}]\n')
                #simple_smiles_to_similar_smiles_pipeline.set(simplified_SMILES, final)

            #self.print_debug_msg(f'Dumping to chemical substance id to lookup db ...')
            #id_to_simple_smiles_pipeline.execute()

            #self.print_debug_msg(f'Dumping to simple SMILES to array of similar SMILES db ...')
            #simple_smiles_to_similar_smiles_pipeline.execute()
        except Exception as e:
            self.print_debug_msg(f'Exception thrown: {e}')
            rv = False

        # return to the caller
        return rv

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

    def print_debug_msg(self, msg: str):
        """ Prints a debug message if enabled in the config file """
        if self._config['debug_messages'] == 1:
            print(msg)

    @staticmethod
    def get_redis(config: json, db_id: int):
        """ Returns a connection to a redis instance """
        return redis.StrictRedis(host=config['redis_host'],
                                 port=int(config['redis_port']),
                                 db=db_id,
                                 password=config['redis_password'])

    def get_simplified_smiles_for_chemicals(self) -> pd.DataFrame:
        """ This method gets SMILES for every chemical substance in the robokop neo4j graph database and creates a simplified SMILES from each.
            The simplified SMILES values will be used as a grouping mechanism and saved in the redis database.
        """

        # Create a target data frame for the processed data
        df: pd.DataFrame = pd.DataFrame(columns=['chem_id', 'original_SMILES', 'simplified_SMILES'])

        try:
            # Create the query. This is of course robokop specific
            # DEBUG: to return 2 that have the same simplified SMILES use this in the where clause ->  and (c.id="CHEBI:140593" or c.id="CHEBI:140451")
            # DEBUG: to return 2 that have the exact same SMILES use this in the where clause ->   and (c.id="CHEBI:85764" or c.id="CHEBI:140773")
            c_query: str = f'match (c:chemical_substance) where c.smiles is not NULL and c.smiles <> "" and c.smiles <> "**" and c.smiles <> "*" RETURN c.id, c.smiles order by c.smiles {self._debug_record_limit}'

            # execute the query
            records: list = self.run_neo4j_query(c_query)

            # did we get some records
            if len(records) > 0:
                # loop through the records
                for r in records:
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

                        record = {'simplified_SMILES': simplified_smiles, 'chem_id': r['c.id'], 'original_SMILES': r['c.smiles']}

                        # append the record to the data frame
                        #df = df.append(record, ignore_index=True)
                    except Exception as e:
                        # alert the user that something was discovered in the original graph record
                        self.print_debug_msg(f"Error - Could not get a simplified SMILES for chem id: {r['c.id']}, Original SMILES: {r['c.smiles']}, Exception: {e}")

                # get the simplified SMILES in groups
                df = df.set_index('simplified_SMILES').groupby('simplified_SMILES')
        except Exception as e:
            raise e

        # return to the caller
        return df

    @staticmethod
    def rdkit_logging(level: Rkl):
        """ Adjusts RDKit logging level. """

        # get a reference to the rdkit logger
        logger = Rkl.logger()

        # set the new loggin level
        logger.setLevel(level)

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
