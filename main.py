'''Script to analyze NIAS from plastic data

First provide a settings file and a dataset. 
Dataset files should be in mzML format if MZMine is used, otherwise MGFs are required.
Generate a workflow instance (a combination of dataset and settings) with the workflowrunner class.
Then run those instances with workflowrunner.run()
Repeat for all settings and datasets as required

'''

import json
import pickle
from time import ctime
from typing import Any,Iterable
import os.path
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

failed_imports = []
try:
    import MS2LDA
except ImportError:
    failed_imports.append("MS2LDA")


from src.myworkflow.integrating import *
from src.myworkflow.parsing import *
# from src.myworkflow.visualisation import *
# from src.myworkflow.unrefined import *
from src.myworkflow.running import * #run_mzmine, run_sirius, run_mzmine run_toxtree,run_classyfire, run_ms2lda_post_to_classyfire, check_classyfire_done, get_classyfire_results
from src.myworkflow.util import validate_dictkeys,json_from_string_or_file,get_filelist_from_folder,convert_missing,integrate_df_cols_to_df
from src.myworkflow.external import batched,filter_component,prune_component,get_edges_of_component

class WorkflowRunner():
    '''Holds the full workflow

    Used to  use the workflow for a specific combination of parameters and settings. 

    Parameters
    ----------
    settings : str | dict
        fn of settings file or dict in json format
        should contain settings defined in workflowsettings required settings class property     

    Attributes
    ----------
    available_tools : list
        Defines which tools can be run. runtool should also be set True in settings file to actually run the tool.
    available_integrations : list
        Defines which integrations can be selected. Integrate_tool should also be set True in settings file to actually integrate the tool.

    '''
    
    available_tools: list[str] = ["mzmine","ms2lda","sirius","toxtree","classyfire", "matchms"]
    available_integrations: list[str] = ["mzmine","ms2lda","sirius","toxtree","classyfire", "matchms","plastchemdb","sirius_db"]

    def __init__(self,settings_json: str | dict) -> None:
        self.settings: WorkflowSettings = WorkflowSettings(settings_json,available_tools=WorkflowRunner.available_tools,available_integrations=WorkflowRunner.available_integrations)
        self.name = self.settings.name
        self.output_folder = self.settings.output_folder
        self.graph: nx.Graph = nx.Graph
        self.network_df, self.network_edgelist = pd.DataFrame,nx.edgelist

    def get_internal_settings(self,settings_json: str | dict) -> dict:
        """Check if internal settings were provided and return them as a json dict
        
        Parameters
        ----------
        settings_json : str | dict
            fn of settings file or settings dict in json format

        Returns
        -------
        Settings dict in json format
        """
        try:
            fn = json_from_string_or_file(settings_json)["paths"].get("internal_settings")
        except KeyError:
            print("settings file must contain 'internal_setting' path in paths")
        return json.load(fn)
    
    def run_all(self) -> None:
        self.run_and_integrate_per_depth()
        print("Sucesfully ran and integrated all tools")

    def get_depths(self,tools: Iterable,goal:str) -> dict[Any,str]:
        """"""
        depths: dict = {}
        for tool in tools:
            config = self.settings.config.get(tool,{})
            goal_depth = config.get(goal,{}).get("depth")
            depth = goal_depth if goal_depth else config.get("depth","N/A")
            match depth:
                case  "1" | 1:
                    depths[tool] = 1
                case  "2" | 2:
                    depths[tool] = 2
                case  "3" | 3:
                    depths[tool] = 3
                case  "4" | 4:
                    depths[tool] = 4
                case  "N/A":
                    depths[tool] = "N/A"
        return depths

    def run_and_integrate_per_depth(self) -> None:
        """For every depth in ascending order, run and integrate tools that are of that depth, finally produce the integrated network graphml file.

        Every tool should have a depth for running and integration specified in internal settings. For every depth, first the tools will be ran, then integrated. 
        This is only attempted for tools that were selected in class property `available tools` or `available integrations`. 
        As the first thing on depth 3, the graph is made so any tools requiring the graph should be of depth 3+.
        As the first thing on depth 4, the network is made so any tools requiring the network should be of depth 4+.

        """

        running_depths: dict[Any,str] = self.get_depths(self.settings.to_run,goal="running")
        integration_depths: dict[Any,str] = self.get_depths(self.settings.to_integrate,goal="integration")
        for current_depth in [1,2,3,4,"N/A"]:
            running: set[Any] = {tool for tool,depth in running_depths.items() if depth==current_depth}
            integrating: set[Any] = {tool for tool,depth in integration_depths.items() if depth==current_depth}               
            if current_depth == 3:
                self.graph = nx.read_graphml(self.settings.paths.get("base_network"))
            if current_depth == 4:   
                self.network_df, self.network_edgelist = network_to_edgelist_and_nodes_df(self.graph)
                self.integration_col = "smiles"  
                self.network_df["smiles"] = self.network_df.apply(self.get_consensus_smiles,axis=1)
                self.network_df["Molecular formula"] = self.network_df.apply(self.get_consensus_formula,axis=1)
            for tool in running:
                self.run_tool(tool)
            for integration in integrating:
                self.integrate_tool(integration)

        for cf_class in ["subclass","class","superclass"]:
                    self.network_df[cf_class] = self.network_df.apply(lambda x: self.get_consensus_class(row=x,classtype=cf_class),axis=1)
       # target_df[target_df["Molecular formula"].notna()].merge(source_df[source_df["molecular_formula"].notna()],left_on="Molecular formula",right_on="molecular_formula")
        self.produce_integrated_graphml()
        
    def get_consensus_formula(self,row):
        """Get formula from database sirius formula if exists, otherwise normal formula"""
        # if row["library_smiles"]
        if not row.get("sirius:molecularFormula"): 
            row["sirius:molecularFormula"] = "N/A"
        if row.get("sirius_db:molecularFormula"):
            return row["sirius_db:molecularFormula"] if isinstance(row["sirius_db:molecularFormula"],str) else row["sirius:molecularFormula"]
        else:
            return row["sirius:molecularFormula"] 
        
    def get_consensus_smiles(self,row):
        """Get smiles from database sirius smiles if exists, otherwise normal smiles"""
        # if row["library_smiles"]
        if not row.get("csifingerid:smiles"): 
            row["csifingerid:smiles"] = "N/A"
        if row.get("csifingerid_db:smiles"):
            return row["csifingerid_db:smiles"] if isinstance(row["csifingerid_db:smiles"],str) else row["csifingerid:smiles"]
        else:
            return row["csifingerid:smiles"] 
    
    def get_consensus_class(self,row,classtype):
        if not row.get(f"canopus:CF_{classtype}"): 
            row[f"canopus:CF_{classtype}"] = "N/A"
        if row.get(f"cf:{classtype}"):
            return row[f"cf:{classtype}"] if isinstance(row[f"cf:{classtype}"],str) else row[f"canopus:CF_{classtype}"]
        else: 
            return row[f"canopus:CF_{classtype}"]
    
    def get_features_per_motif(self,motif_col: pd.DataFrame) -> dict[str,list[str]]:
        features_per_motif= {}
        for feature_id, motifs in motif_col.to_dict().items():
            for motif in motifs:
                if not motif in features_per_motif:
                    features_per_motif[motif] = []
                    features_per_motif[motif].append(feature_id)
                if not feature_id in features_per_motif[motif]:
                    features_per_motif[motif].append(feature_id)
        return features_per_motif

    def produce_integrated_graphml(self) -> None:
        """Cleans, filters and creates graphml from network DF and edges
        
        Topology filter is always used: cluster can be only of 100 nodes. 
        Converts all missing values to "N/A" and turns all properties into strings for easier reading in cytoscape!
        """
        integrated_df: pd.DataFrame = self.network_df
        integrated_network = nx.from_pandas_edgelist(self.network_edgelist,edge_attr=True)
        integrated_df_as_strings = integrated_df.astype(str)
        integrated_df_cleaned = integrated_df_as_strings.map(convert_missing)
        nodes_dict = integrated_df_cleaned.to_dict(orient='index')
        for node_ID in nodes_dict:
            if not node_ID in integrated_network:
                integrated_network.add_node(node_ID)
            #nodes_dict[node_ID] = {}
            for attribute in nodes_dict[node_ID]:
                integrated_network.nodes[node_ID][attribute] = nodes_dict[node_ID][attribute]
        filter_component(integrated_network,100)

        nx.write_graphml(integrated_network, f"{self.output_folder}/{self.name}_premotifs.graphml")

        df_motifs = integrated_df["MS2LDA:allmotifs"][integrated_df["MS2LDA:allmotifs"] != "N/A"]
        features_per_motif = self.get_features_per_motif(df_motifs)
        for motif in features_per_motif:
            combinations = list(itertools.combinations(features_per_motif[motif],2))
            # print(combinations)
            combinations = [tuple([combination[0],combination[1],{motif:True}]) for combination in combinations]
            integrated_network.add_edges_from(combinations,motif=motif,type="motif")
        nx.write_graphml(integrated_network, f"{self.output_folder}/{self.name}_postmotifs.graphml")

    def run_tool(self,tool: Any) -> None:
        """Get running parameters and run selected tool
        
        Here, all running behavirous for a specific tool is defined. First the required parameters are obtained from properties that were set & checked in `WorkflowSettings`.
        Simple format converions of input are done here as well.  
        No order is defined here: use depth in internal settings for that.
        Will only do something for tools defined in `available_tools` class property.

        """
        match tool:
            case "classyfire":
                #print('ran classyfire (not, debug)')
                input_csv = self.create_smiles_csv_from_df(output_name = f"{self.output_folder}/classyfire_input.csv",header=False)
                output_path = self.settings.paths.get("classyfire_output")
                print('running classyfire')
                run_classyfire(input_path=input_csv,output_path=output_path)
                # parse_classyfire_sdf("/lustre/BIF/nobackup/hendr218/temp/sdftest.txt")
            case "matchms":
                print("running matchms")
                try:
                    graph: nx.Graph = create_network_from_mgf(self.settings.paths.get("input_mgf"))
                except ValueError:
                    raise InputError("Input MGF was not valid")
                output_path = self.settings.paths.get("base_network")
                graph.export_to_graphml(output_path)
                #nx.write_graphml(graph, output_path)
            case "ms2lda":
                output_path: str = self.settings.paths.get("ms2lda_output")
                self.settings.ms2lda["dataset_parameters"]["output_folder"] = output_path
                ms2lda_params: dict = self.settings.get_ms2lda_params()
                input_mgf: str = self.settings.paths.get("input_mgf")
                #print("ran ms2lda (not, debug)")
                print('running ms2lda')
                run_ms2lda(dataset=input_mgf,params=ms2lda_params)
            case "mzmine":
                mzmine_location = self.settings.paths.get("mzmine_location")
                mzmine_userfile_location = self.settings.paths.get("mzmine_userfile_location")
                base_batchfile = self.settings.paths.get("mzmine_base_batchfile")
                output_loc = f'{self.output_folder}/mzmine'
                #get filelist if it was provided, else fall back to the folder
                self.data: list = self.settings.paths.get("file_list", get_filelist_from_folder(self.settings.paths.get("data_folder")))
                batchfile_writer = mzbatch_writer(base_mzbatch=base_batchfile,data=self.data,mzmine_params=self.settings.mzmine)
                batchfile = batchfile_writer.write_mzbatch(output_filename=f'{output_loc}_mzbatch.mzbatch')
                #print("ran mzmine, not, debug")
                print('running mzmine')
                run_mzmine(batchfn=batchfile,mzmine_output_loc=output_loc,mzmine_location=mzmine_location,mzmine_userfile_location=mzmine_userfile_location,temp_folder=f'{output_loc}/temp/')
            case "sirius":
                output_loc: str = self.settings.paths.get("sirius_output")
                sirius_path: str = self.settings.paths.get("sirius_path")
                input_mgf: str= self.settings.paths.get("input_mgf")
                instrument: str = self.settings.sirius.get("instrument")
                formula_db: str = self.settings.sirius.get("formula_db")
                #print("ran sirius (not, debug)")
                print('running sirius')
                run_sirius(input_path=input_mgf,output_path=output_loc,sirius_path=sirius_path,formula_db=formula_db,instrument=instrument)
            case "toxtree":
                toxtree_loc = self.settings.paths.get("toxtree_path")
                toxtree_module = self.settings.toxtree.get("module")
                config: dict = self.settings.config.get("toxtree",{}).get("running",{}).get("translations",{})
                if not toxtree_module in config.get("modules",{}):
                    raise Settingserror(f"toxtree module {toxtree_module} not supported, add to module translation in internal settings to support")
                toxtree_module = config.get("modules",{}).get(toxtree_module)
                input_csv = self.create_smiles_csv_from_df(output_name = f"{self.output_folder}/toxtree_input.csv",header=True)
                output_path = f"{self.output_folder}/toxtree_results.csv"
                #print("ran toxtree (not, debug)")
                print('running toxtree')
                run_toxtree(input_path=input_csv,output_path=output_path,toxtree_path=toxtree_loc,module_path=toxtree_module)
       
    def integrate_tool(self,tool: Any) -> None:
        """Get needed paramaters and integrate selected tool
        
        Here, all integration behaviour for a specific tool is defined. First the required parameters are obtained from properties that were set & checked in `WorkflowSettings`.
        Simple format converions of input are done here as well.  
        No order is defined here: use depth in internal settings for that.
        Will only do something for tools defined in `available_integrations` class property.

        """
        match tool:
            case "classyfire":
                print('integrating classyfire')
                output_folder: str = self.settings.paths.get("classyfire_output")
                classyfire_records: pd.DataFrame= parse_classyfire_sdf(output_folder)
                config: dict  = self.settings.config.get("classyfire", {}).get("integration", {})
                integration_settings: tuple[dict,str] = get_merging_settings(config,self.integration_col)
                self.network_df: pd.DataFrame = integrate_df_cols_to_df(self.network_df,classyfire_records, *integration_settings)
            case "ms2lda":
                print('integrating ms2lda')
                output_folder = self.settings.paths.get("ms2lda_output")
                dataset = self.settings.paths.get("input_mgf")
                output_database = self.settings.paths.get("ms2lda_results_db")
                self.network_df = integrate_ms2lda_to_df(self.network_df,dataset,output_folder,output_database)
            case "mzmine":
                #temp to get some quick results, will be improved later dataframe: pd.DataFrame,input_mgf: str,metadata_csv: str,quant_table: str
                print('integrating mzmine')
                input_mgf = self.settings.paths.get("sirius_mgf")
                quant_table = self.settings.paths.get("quant_table")
                metadata_csv = self.settings.paths.get("metadata_csv")
                self.network_df = temp_metadata_adder(self.network_df,input_mgf,metadata_csv,quant_table)
            case "sirius":
                print('integrating sirius')
                config: dict = self.settings.config.get("sirius",{}).get("integration",{}).get("translations",{})
                translators = config.get("sirius",{}),config.get("csi:fingerid",{}),config.get("canopus",{})
                outputs = self.settings.paths.get("sirius_tool_output"), self.settings.paths.get("csi:fingerid_output"), self.settings.paths.get("canopus_output")
                self.graph = integrate_sirius_to_graph(graph=self.graph,translators=translators,outputs=outputs)
            case "sirius_db":
                print('integrating sirius from db')
                config: dict = self.settings.config.get("sirius_db",{}).get("integration",{}).get("translations",{})
                translators = config.get("sirius",{}),config.get("csi:fingerid",{}),config.get("canopus",{})
                outputs = self.settings.paths.get("sirius_tool_db_output"), self.settings.paths.get("csi:fingerid_db_output"), self.settings.paths.get("canopus_db_output")
                self.graph = integrate_sirius_to_graph(graph=self.graph,translators=translators,outputs=outputs)
            case "toxtree":
                print('integrating toxtree')
                toxtree_output: str = self.settings.paths.get("toxtree_output")
                toxtree_df: pd.DataFrame = parse_cramer_classifications(toxtree_output)
                config: dict  = self.settings.config.get("toxtree", {}).get("integration", {})
                integration_settings = get_merging_settings(config,self.integration_col)
                self.network_df = integrate_df_cols_to_df(self.network_df,toxtree_df, *integration_settings)
            case "plastchemdb":
                print('integrating plastchem')
                config: dict  = self.settings.config.get("plastchemdb", {}).get("integration", {})
                plastchem_db: str = self.settings.paths.get("plastchem_path")
                plastchem_df =  pd.read_csv(plastchem_db, sep='\t', encoding='windows-1251', low_memory=False, header=1)
                integration_settings = get_merging_settings(config,self.integration_col) 
                self.network_df = integrate_df_cols_to_df(self.network_df,plastchem_df, *integration_settings)

    def create_smiles_csv_from_df(self,header: bool,output_name: str) -> bool:
        """Creates a smiles file from a dataframe column with rows on newlines"""
        network_df = self.network_df
        #df_smiles = network_df[network_df['csifingerid:smiles'].notna()].rename(columns={"csifingerid:smiles": "smiles"})      
        df_smiles = network_df[network_df['smiles'].notna()]  
        df_smiles.to_csv(output_name,columns=["smiles"],header=header, index=False,na_rep="N/A")
        return output_name 
    
class Settingserror(Exception):
    """Error for indicating invalid settings files"""
    pass

class InputError(Exception):
    """Error for indicating invalid input files"""
    pass

class WorkflowSettings:
    """Class to process and hold settings for the workflow

    First check if the always required fields from `required_settings' and `required_paths` class attributes are present. 
    If so, loads paths and adds new ones depending on the selected tools. (i.e. if a tool creates a file that might be needed downstream, it should add it to paths)
    Then, checks if all tools have the input paths that they need as defined in the internal_settings file
    Finally, tests if 

    Parameters
    ----------
    settings : str
        fn of settings file
    available_tools : list[str]
        tools supported by worklflow: others will be ignored
    available_integrations : list[str]
        integrations supported by worklflow: others will be ignored

    Attributes
    ----------
    required_settings : list[str]
        settings that must be provided in settings json 
    required_settings tools : list[str]
        paths that must be provided in settings settings json paths field

    """
    required_settings: list[str] = ["paths","run_tools","integrate_tools"]
    required_paths: list[str] = ["base_output_folder","internal_settings"]

    def __init__(self,settings_json:str | dict,available_tools: list[str],available_integrations: list[str]) -> None:
        self.input: dict = json_from_string_or_file(settings_json)   
        self.validate_required_settings(settings=self.input)
        self.config: dict = json_from_string_or_file(self.input["paths"]["internal_settings"])

        #Add properties for convenience in use and code
        self.paths: dict = self.input["paths"]
        self.name = self.paths.get("name",ctime().replace(" ","_"))
        self.output_folder = f'{self.paths.get("base_output_folder")}/{self.name}'
        Path(self.output_folder).mkdir(exist_ok=True)
        self.to_run: list[str] = self.select_used_tools(available_tools=available_tools,setting="run_tools")
        self.to_integrate: list[str] = self.select_used_tools(available_tools=available_integrations,setting="integrate_tools")

        #Let tools set additional settings and then check if all requirements are met
        self.set_tool_settings(set(self.to_run+self.to_integrate))
        self.check_tool_requirements(tools=self.to_run,goal="running")
        self.check_tool_requirements(tools=self.to_integrate,goal="integration")

        self.check_tool_requirements(tools=["all_tools"],goal="running")
        self.check_tool_requirements(tools=["all_tools"],goal="integration")


    def validate_required_settings(self,settings: dict) -> None:
        """Check if settings json contains required fields"""
        if not validate_dictkeys(settings,WorkflowSettings.required_settings):
            raise Settingserror(f"settings json must contain {', '.join(WorkflowSettings.required_settings)}")
        if not validate_dictkeys(settings["paths"],WorkflowSettings.required_paths):
            raise Settingserror(f"settings json paths must contain {', '.join(WorkflowSettings.required_paths)}")
    

    def check_tool_requirements(self,tools: Iterable[str],goal: str) -> None:
        for tool in tools:
            requirements = self.get_tool_requirements(tool,goal)
            self.validate_tool_requirements(tool,requirements)               
    
    def set_tool_settings(self,tools: Iterable[str]) -> None:
        """set tool settings as property and add any settings that depend on the tool""" 
        for tool in tools:
            tool_settings = self.input.get(tool,{})
            match tool:
                case "mzmine":
                    self.mzmine = tool_settings
                    self.paths["input_mgf"] = f'{self.output_folder}/mzmine/mzmine_iimn_fbmn.mgf'
                    self.paths["sirius_mgf"] = f'{self.output_folder}/mzmine/mzmine_sirius.mgf'
                    self.paths["quant_table"] = f'{self.output_folder}/mzmine/mzmine_iimn_fbmn_quant_full.csv'
                case "sirius":
                    self.sirius = tool_settings
                    sirius_loc = f'{self.output_folder}/sirius/'

                    Path(sirius_loc).mkdir(exist_ok=True)
                    
                    if "sirius" in self.to_run:
                        self.paths["csi:fingerid_output"] = f"{sirius_loc}structure_identifications.tsv"
                        self.paths["sirius_tool_output"] = f"{sirius_loc}formula_identifications.tsv"
                        self.paths["canopus_output"] = f"{sirius_loc}canopus_structure_summary.tsv"

                    self.paths["sirius_output"] = sirius_loc
                case "sirius_db":
                    self.siriusdb = tool_settings
                    siriusdb_loc = f'{self.output_folder}/sirius_db/'
                    self.paths["canopus_db_output"] = f"{siriusdb_loc}canopus_structure_summary.tsv"
                    self.paths["csi:fingerid_db_output"] = f"{siriusdb_loc}structure_identifications.tsv"
                    self.paths["sirius_tool_db_output"] = f"{siriusdb_loc}formula_identifications.tsv"
                    self.paths["sirius_db_output"] = siriusdb_loc
                case "matchms":
                    self.matchms = tool_settings
                    self.paths["base_network"] = f'{self.output_folder}/base_network.graphml'
                case "toxtree":
                    self.toxtree = tool_settings
                    self.paths["toxtree_output"] = f"{self.output_folder}/toxtree_results.csv"
                    print("warning: classyfire requires SMILES to be in the dataframe, currently checks for that are not implemented")
                case "classyfire":
                    self.classyfire = tool_settings
                    self.paths["classyfire_output"] = f"{self.output_folder}/classyfire_results.sdf"
                    print("warning: classyfire requires SMILES to be in the dataframe, currently checks for that are not implemented")
                case "ms2lda":
                    self.ms2lda = tool_settings
                    ms2lda_loc = f"{self.output_folder}/ms2lda"
                    self.paths["ms2lda_output"] = f"{ms2lda_loc}"
                    self.paths["ms2lda_results_db"] = f"{self.output_folder}/ms2lda/motifDB_optimized.xlsx"

    def get_tool_requirements(self,tool: str, goal:str) -> tuple[list,list,list,list]:
        """Get requirements if available, otherwise return emoty list"""
        config = self.config.get(tool,{}).get(goal,{}).get("requirements",{})
        required_paths: Iterable[Any] = config.get("paths",[])
        required_settings: Iterable[Any] = config.get("settings",[])
        required_optional_paths: Iterable[Iterable[Any]] = config.get("optional_paths",[])
        required_optional_settings: Iterable[Iterable[Any]] = config.get("optional_settings",[])
        required_modules: Iterable[Iterable[Any]] = config.get("modules",[])
        return required_paths,required_settings,required_optional_paths,required_optional_settings,required_modules
    
    def validate_tool_requirements(self,tool, requirements):
        """Raise errors if required settings are not present in settings"""
        required_paths,required_settings,required_optional_paths,required_optional_settings,required_modules = requirements
            # if not validate_dictkeys(self.input["paths"],required_paths):
        if not validate_dictkeys(self.paths,required_paths):
            raise Settingserror(f"all of {required_paths} required to run {tool}, but one is missing in paths settings")
        if not validate_dictkeys(self.input.get(tool),required_settings):
            raise Settingserror(f"One of {required_settings} missing in {tool} settings")
        for combination in required_optional_paths:
            if not validate_dictkeys(self.paths,combination):
                    raise Settingserror(f"Any one of {combination} required to run {tool}, but one is missing in paths settings")
        for combination in required_optional_settings:
            if not validate_dictkeys(self.input.get(tool),combination):
                    raise Settingserror(f"Any of {combination} is required in {tool} settings")
        for module in required_modules:
            if not module in sys.modules:
                raise Settingserror(f"all of {required_modules} required to run {tool}, but one is missing in imported modules")

    def select_used_tools(self,available_tools: list,setting: str) -> list:
        """Check which available tools were selected in settings file"""
        return [tool for tool in available_tools if self.input[setting].get(tool) == "True"]
    
    def get_ms2lda_params(self):
        """Get the parameters ms2lda needs froms settings"""
        preprocessing_parameters = self.ms2lda["preprocessing_parameters"]
        convergence_parameters = self.ms2lda["convergence_parameters"]
        annotation_parameters = self.ms2lda["annotation_parameters"]
        model_parameters = self.ms2lda["model_parameters"]
        train_parameters = self.ms2lda["train_parameters"]
        fingerprint_parameters = self.ms2lda["fingerprint_parameters"]
        dataset_parameters = self.ms2lda["dataset_parameters"]
        n_motifs = self.ms2lda["n_motifs"]
        n_iterations = self.ms2lda["n_iterations"]
        motif_parameter = self.ms2lda["motif_parameter"]

        return preprocessing_parameters,convergence_parameters,annotation_parameters,model_parameters,train_parameters,fingerprint_parameters,dataset_parameters,n_iterations,n_motifs,motif_parameter

settings_path: str = "/lustre/BIF/nobackup/hendr218/mycode/src/myworkflow/settings_pcdb.json"
# # # workflow_dict: dict[str,str] = {"A": settings_path,"B":settings_path}
# # # workflow_holder: list[WorkflowRunner] = [WorkflowRunner(dataset,settings) for dataset,settings in workflow_dict.items()]
# # # for workflow in workflow_holder:
# # #     workflow.run_all()
test_workflow = WorkflowRunner(settings_path)
test_workflow.run_all()
test_workflow = WorkflowRunner(settings_path)
test_workflow.run_all()
# loc = f"{test_workflow.output_folder}/{test_workflow.name}.pickle"
# with open(loc,"wb") as file:
#     pickle.dump(test_workflow,file)

with open("/lustre/BIF/nobackup/hendr218/Data/with_pcdb_copy/first_test/first_test.pickle","rb") as file:
    wf = pickle.load(file)
# # print('ah')


def create_counts_df(df: pd.DataFrame,col_to_plot: str,col_to_filter:str=None,normalize:bool=True) -> pd.DataFrame:
    """Creates a dataframe for"""
    print(col_to_filter)
    df = df[df[col_to_plot]!= "N/A"] 
    if col_to_filter:
        df = df[col_to_plot][(df[col_to_filter] == True) & (df["is_blank"] == False)].value_counts(normalize=normalize)
    else: 
        df = df[col_to_plot][df["is_blank"] == False].value_counts(normalize=normalize)
    df=df.to_dict()
    df = pd.DataFrame.from_dict(df,orient="index")
    if col_to_filter is None:
        return df.rename(columns={0:"All"})
    return df.rename(columns={0:col_to_filter})

def join_counts(*dataframes: Iterable[pd.DataFrame]):
    if len(dataframes) == 1:
        return dataframes
    dfs = [dataframe for dataframe in dataframes]
    df_joined = dfs.pop(0)
    for i,df in enumerate(dfs):
        df_joined = df_joined.join(df,how="outer") #lsuffix=f"_{i+1}",rsuffix=f"_{i+2}"
    return df_joined.fillna(0)

def df_to_barplot(df:pd.DataFrame,columns:list[str],title: str,ytitle:str,fn: str) -> None:
    print(f"df to barplot: {columns}")
    ax = df[columns].plot(kind='bar', title=title, figsize=(15, 10), legend=True, fontsize=12)
    ax.set_ylabel(ytitle, fontsize=12)
    plt.tight_layout()
    plt.savefig(fn)


def df_to_counts_plot(df:pd.DataFrame,filter_cols: list[str],vis_col: str,fn:str,title:str,relative_counts: bool = True,include_nofilter:bool=False) -> None:
    count_dfs = [create_counts_df(df,vis_col,filter_col,normalize=relative_counts) for filter_col in filter_cols]
    cols = filter_cols
    print(f"df to counts plot:{cols} type: {type(cols)}")
    if include_nofilter:
        print(f" including nofilter")
        count_dfs.append(create_counts_df(df,vis_col,normalize=relative_counts))
        cols.append("All")
    joined_df = join_counts(*count_dfs)
    ytitle = "Feature count"
    if relative_counts:
        joined_df = joined_df.apply(lambda x: x*100)
        ytitle = "occurence (%)"
    print(f"df to counts plot:{cols}")
    df_to_barplot(joined_df,columns=cols,title=title,ytitle=ytitle,fn=fn)

def plot_classes(df,loc,relative_counts=True,filter_type="type",class_type="all"):
    filter_dict = {
        "type":["PE","PE_PET","PE_PA"],
        "state":["Removed_by_dec","Introduced_by_dec","Kept_by_dec"]
    } 
    class_dict={
        "canopus":["canopus:CF_subclass", "canopus:CF_class", "canopus:CF_superclass"],
        "classyfire":["CF:subclass", "CF:class", "CF:superclass"],
        "all":["subclass", "class", "superclass"]
    }
    vis_cols = class_dict[class_type]
    filter_cols = filter_dict[filter_type]

    counttype = "relative" if relative_counts else "absolute"
    for vis_col in vis_cols:
        print(F"vis col: {vis_col}")
        fn = f"{loc}/{class_type}_classes_by_{filter_type}_{vis_col}_{counttype}_counts.png"
        print(f"plot_classes pre all: {filter_cols}")
        df_to_counts_plot(df,filter_cols,vis_col,fn=fn,title=vis_col,relative_counts=relative_counts)
        fn = f"{fn[:-4]}_all.png"
        print(f"plot_classes post all: {filter_cols}")
        df_to_counts_plot(df,filter_cols,vis_col,fn=fn,title=vis_col,relative_counts=relative_counts,include_nofilter=True)
    
loc = "/lustre/BIF/nobackup/hendr218/Data/plot_tests"
df = wf.network_df
plot_classes(df,loc,relative_counts=True,filter_type="type",class_type="all")
# for filter_type in ["type","state"]:
#     for class_type in ["canopus","classyfire","all"]:
#         plot_classes(df,loc,relative_counts=True,filter_type=filter_type,class_type=class_type)
#         plot_classes(df,loc,relative_counts=False,filter_type=filter_type,class_type=class_type)
#         plt.close()

print('hi')