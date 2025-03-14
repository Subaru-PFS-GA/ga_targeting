import os
import sys
from shutil import copyfile
from nbparameterise import extract_parameters, parameter_values, replace_definitions
import nbformat
from nbconvert.preprocessors import ClearOutputPreprocessor, ClearMetadataPreprocessor
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError
from nbconvert import HTMLExporter

from ..setup_logger import logger

class NotebookRunner():
    def __init__(self, orig=None):

        if not isinstance(orig, NotebookRunner):
            self.__workdir = None
            self.__parameters = {}
            self.__kernel = None
        else:
            self.__workdir = orig.__workdir
            self.__parameters = orig.__parameters
            self.__kernel = orig.__kernel

        self.__nb = None
        self.__notebook_path = None

    #region Properties

    def __get_workdir(self):
        return self.__workdir
    
    def __set_workdir(self, workdir):
        self.__workdir = workdir

    workdir = property(__get_workdir, __set_workdir)

    def __get_parameters(self):
        return self.__parameters
    
    def __set_parameters(self, parameters):
        self.__parameters = parameters

    parameters = property(__get_parameters, __set_parameters)

    def __get_kernel(self):
        return self.__kernel
    
    def __set_kernel(self, kernel):
        self.__kernel = kernel

    kernel = property(__get_kernel, __set_kernel)

    #endregion

    def open_ipynb(self, fn):
        with open(fn, 'r') as f:
            self.__nb = nbformat.read(f, as_version=4)
            self.__notebook_path = fn

    def save_ipynb(self, fn):
        with open(fn, 'w') as f:
            nbformat.write(self.__nb, f)

    def save_html(self, fn):
        html_exporter = HTMLExporter()
        #html_exporter.template_file = 'basic'
        (body, resources) = html_exporter.from_notebook_node(self.__nb)
        with open(fn, 'w') as f:
            f.write(body)

    def __clear_all_output(self):
        cpp = ClearOutputPreprocessor()
        self.__nb, resources = cpp.preprocess(self.__nb, None)

    def __execute(self):
        cwd = os.getcwd()
        if self.__workdir is not None:
            os.chdir(self.__workdir)

        resources = {
            'metadata': {
                'path': self.__workdir
            }
        }
        epp = ExecutePreprocessor(timeout=None)
        try:
            self.__nb, resources = epp.preprocess(self.__nb, resources)
        except Exception as ex:
            logger.error(f'An error has occurred while executing notebook `{self.__notebook_path}`')
            # Error is not propagated to allow saving notebook
        finally:
            os.chdir(cwd)

    def run(self):        
        # Set kernel
        if self.__kernel is not None:
            ks = self.__nb.metadata.get('kernelspec', {})
            ks['name'] = self.__kernel

        # Set parameters
        orig_params = extract_parameters(self.__nb)
        new_params = parameter_values(orig_params, **self.__parameters)
        self.__nb = replace_definitions(self.__nb, new_params, execute=False)

        self.__clear_all_output()
        self.__execute()
