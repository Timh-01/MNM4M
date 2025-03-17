from main import WorkflowRunner,WorkflowSettings,validate_dictkeys,Settingserror
import pytest
import types 
import json


@pytest.fixture
def my_settings():
    return {
"paths":{
    "base_output_folder":"/lustre/BIF/nobackup/hendr218/Data/testrun_28-01",
    "mzmine_location":"/lustre/BIF/nobackup/hendr218/Programs/MZMine/bin/mzmine",
    "mzmine_userfile_location":"/lustre/BIF/nobackup/hendr218/Programs/MZMine/tim.mzuser",
    "mzmine_base_batchfile":"",
    "toxtree_module_path":"path",
    "toxtree_path":"path"
    },
"run_tools":{
    "sirius": "True",
    "toxtree": "True",
    "classyfire": "True",
    "ms2lda": "True",
    "mzmine": "True"
    },
"mzmine": {

},
"toxtree": {

},
"classyfire":{

},
"sirius":{

},
"ms2lda":{
    "preprocessing_parameters": {
        "min_mz": 0,
        "max_mz": 1000,
        "max_frags": 1000,
        "min_frags": 3,
        "min_intensity": 0.01,
        "max_intensity": 1
        },
    "convergence_parameters": {
        "step_size": 50,
        "window_size": 10,
        "threshold": 0.001,
        "type": "perplexity_history"
        },
    "annotation_parameters": {
        "criterium": "best",
        "cosine_similarity": 0.70,
        "n_mols_retrieved": 5 
        },
    "model_parameters": {
        "rm_top": 0, 
        "min_cf": 0, 
        "min_df": 3, 
        "alpha": 0.6, 
        "eta": 0.01,
        "seed": 42
        },
    "train_parameters": {
            "parallel": 3,
            "workers": 1 
        },
    "fingerprint_parameters": {
        "fp_type": "rdkit",
        "threshold": 0.8
        },
    "dataset_parameters": {
        "acquisition_type": "DDA",
        "significant_digits": 2,
        "charge": 1,
        "name": "pos_all",
        "output_folder": "output_folder"
        },
    "n_motifs": 200,
    "n_iterations": 2000,
    "motif_parameter": 20
    }
}

@pytest.fixture
def my_settings_no_tools():
    return {"paths":{
    "base_output_folder":"/lustre/BIF/nobackup/hendr218/Data/testrun_28-01",
    "mzmine_location":"/lustre/BIF/nobackup/hendr218/Programs/MZMine/bin/mzmine",
    "mzmine_userfile_location":"/lustre/BIF/nobackup/hendr218/Programs/MZMine/tim.mzuser"
    },
"run_tools":{
    "sirius": "False",
    "cramer": "False",
    "classyfire": "False",
    "ms2lda": "False",
    "mzmine": "False"
    },
"other":{
    "placeholder":"valuess"
    }
}

@pytest.fixture
def my_workflowrunner(my_settings):
    return WorkflowRunner(dataset="A",settings_json=my_settings)
    #return WorkflowRunner(dataset="A",settings_json="/lustre/BIF/nobackup/hendr218/mycode/src/myworkflow/Tests/settings.json")

@pytest.fixture
def my_workflowrunner_missingsettings(my_settings):
    missing_settings = my_settings.pop("run_tools")
    return WorkflowRunner(dataset="A",settings_json=missing_settings)
\
@pytest.fixture
def my_workflowrunner_toomanysettings(my_settings):
    toomany_settings = my_settings
    toomany_settings["another"] = {"keyy":"vvalue"}
    return WorkflowRunner(dataset="A",settings_json=toomany_settings)
    
@pytest.fixture
def my_settingsclass(my_settings):
    return WorkflowSettings(settings_json=my_settings) 
#"/lustre/BIF/nobackup/hendr218/mycode/src/myworkflow/Tests/settings.json")
# def test_my_workflowrunner_from_json():
#     # json_settings = "/lustre/BIF/nobackup/hendr218/mycode/src/myworkflow/Tests/settings.json"
#     # assert 5*2 == 10 #WorkflowRunner(dataset="A",settings_json=json_settings).settings_json == json_settings

def test_validate_dictkeys():
    test_dict = {"base_output_folder":"/lustre/BIF/nobackup/hendr218/Data/testrun_28-01",
    "mzmine_location":"/lustre/BIF/nobackup/hendr218/Programs/MZMine/bin/mzmine",
    "mzmine_userfile_location":"/lustre/BIF/nobackup/hendr218/Programs/MZMine/tim.mzuser"}
    test_list = ["base_output_folder","mzmine_location","mzmine_userfile_location"]
    wrong_list = ["base_output_folder","mzmine_location", "notpresent","mzmine_userfile_location"]
    assert validate_dictkeys(test_dict,test_list) 
    assert not validate_dictkeys(test_dict,wrong_list)

def test_mysettingsclass_validate_settingserror(my_settingsclass,my_settings):
    assert my_settingsclass.validate_settings(my_settings)
    with pytest.raises(Settingserror):
        my_settings["other"]={"wrong":"placeholder"}
        my_settingsclass.validate_settings(my_settings)

def test_workflowsettings_extractsettings(my_settingsclass,my_settings):
     paths, tools, other = my_settingsclass.extract_settings()
     assert (my_settings["paths"], my_settings["run_tools"], my_settings["other"]) == (paths,tools,other)

def test_settingsfilepath():
    with pytest.raises(FileNotFoundError):
        WorkflowRunner(dataset="A",settings_json="somepath/invalid_andstuff/error")

def test_get_setting(my_workflowrunner):
    assert my_workflowrunner.settings.get_other("placeholder") == "valuess"
    assert ("nonexistentpath" not in my_workflowrunner.settings.paths) 



@pytest.fixture
def my_workflowrunner_only_cramer(my_settings_no_tools):
    only_require_cramer = my_settings_no_tools
    only_require_cramer["run_tools"]["cramer"] = "True"
    return WorkflowRunner("A", only_require_cramer)

@pytest.fixture
def my_workflowrunner_with_optional(my_settings_no_tools):
    only_require_cramer = my_settings_no_tools
    only_require_cramer["run_tools"]["sirius"] = "True"
    only_require_cramer["other"]["sirius_instrument"] = "test_instrument"
    return WorkflowRunner("A", only_require_cramer)

def test_required_tool_settings(my_workflowrunner_only_cramer):
    assert my_workflowrunner_only_cramer.settings.tools["cramer"] == "True"
    assert my_workflowrunner_only_cramer.settings.tools["sirius"] == "False"
    with pytest.raises(Settingserror):
        my_workflowrunner_only_cramer.run_all()

def test_get_runnable(my_workflowrunner_only_cramer):
    assert len(list(my_workflowrunner_only_cramer.get_runnable())) == 1

def test_optional_tool_settings(my_workflowrunner_with_optional):
    assert my_workflowrunner_with_optional.settings.tools["mzmine"] == "False"
    assert my_workflowrunner_with_optional.settings.other.get("sirius_instrument") == "test_instrument"

def test_run_all(my_workflowrunner):
    assert my_workflowrunner.run_all()


    # only_require_cramer = my_settings_no_tools
    # only_require_cramer["run_tools"]["cramer"] = "True"
    # assert WorkflowRunner(dataset="A",settings_json=only_require_cramer))
    # print(only_require_cramer)
    # with pytest.raises(Settingserror):
    #     WorkflowRunner(dataset="A",settings_json=only_require_cramer).run_all()

# def test_requied_tool_settings(myworkflowrunner):
#     assert 


# def test_optional_tool_settings(myworkflowrunner):
#     assert 
# def test_run_all_selection(my_workflowrunner):
#     assert my_workflowrunner.run_all(), types.GeneratorType




# def test_workflowsettings(my_settingsclass,my_settings):
#      assert my_settingsclass.settings == my_settings

# # def test_mysettingsclass(my_settingsclass):
#     print(my_settingsclass.settings())

# def test_toomanysettings(my_workflowrunner_toomanysettings):
#     with pytest.raises(Settingserror):
#         my_workflowrunner_toomanysettings.set_settings()

# def test_toofew(my_workflowrunner_toofewsettings):
#     with pytest.raises(Settingserror):
#         my_workflowrunner_toofewsettings.set_settings()

# def test_settingsvalidation(my_settingsclass):
#     assert my_settingsclass.validate_settings()

