#include "Python.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifdef __linux
#include <dirent.h>
#elif _WIN32
#include "dirent.h"
#endif

/* BOINC */
#include <boinc_api.h>
#include <boinc_zip.h>

/* Numpy staticly linked */
#include "numpy_builtin.h"
#include "frozen.h"

#include "linker_generator.h"

// CMAKE vars: APP_NAME, APP_VERSION


//#include "config.h"

#define APP_NAME "linker_generator"
#define APP_MAJOR 0
#define APP_MINOR 16

 // Use linux functions
#define getExecPath getExecPath_l
#define createFailIfExists createFailIfExists_l
#define DeleteDirectory DeleteDirectory_l

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define APP_VERSION TOSTRING(APP_MAJOR) "." TOSTRING(APP_MINOR)

#ifdef __linux
#include "linux_functions.h"
#define DIR_SLASH "/"
#elif _WIN32
#include "windows_functions.h"
#pragma comment(lib, "crypt32.lib") 
#define DIR_SLASH "\\"
#endif

#define RESOURCES_ZIP "linkergen_res.zip"
#define STDERR_FILE "stderr.txt"
#define STDOUT_FILE STDERR_FILE
#define OUTPUT_FILE "resultfile.txt"

#define DEBUG


int fileExists(char* name) {
  struct stat buf;
  return (stat(name, &buf) == 0) ? TRUE : FALSE;
}

const char *used_filenames[] = {
								"instructions.txt",
								"atomfile.pdb",
								"\0"
								};


// TODO change to stderr
void print_error(char *place, char *msg, char *param){
	if(param != NULL)
		fprintf(stderr, "ERROR[%-12.12s]: %s %s\n", place, msg, param);
	else
		fprintf(stderr, "ERROR[%-12.12s]: %s\n", place, msg);
}


char *strnchr(char *str, size_t len, int character) {
    char *end = str + len;
    char c = (char)character;
    do {
        if (*str == c) {
            return str;
        }
    } while (++str <= end);
    return NULL;
}


int unzip_resources(char *basedir) {
	char unzip_file[MAX_PATH] = {0};
	char unzip_lock_started[MAX_PATH] = {0};
	char unzip_lock_finished[MAX_PATH] = {0};
	char working_dir[MAX_PATH] = {0};
    int rc;
    FILE *lockfile;
    //unzFile uf;

    /* Check if other process is already extracting */
    strncpy(unzip_lock_started, basedir, MAX_PATH);
	strncat(unzip_lock_started, DIR_SLASH "unzip_started", MAX_PATH - (strlen(basedir) + 1));
    
    strncpy(unzip_lock_finished, basedir, MAX_PATH);
	strncat(unzip_lock_finished, DIR_SLASH "unzip_finished", MAX_PATH - (strlen(basedir) + 1));

	mkdir(basedir,  S_IRWXU | S_IRUSR | S_IWUSR | S_IXUSR);

    if(!fileExists(unzip_lock_finished)) {
    	if(createFailIfExists(unzip_lock_started) < 0){	// Lockfile exists --> other process extracting
    		/* Wait max 5 minutes for process to finish */
    		fprintf(stderr, "Other process allready extracting...\n");
    		boinc_temporary_exit_wrapper(20, "Other process allready extracting\n", FALSE);
    		return -1; 
    	}
    	else {								// Lockfile does not exist --> extract data
    		boinc_resolve_filename(RESOURCES_ZIP, unzip_file, MAX_PATH);
        	if(fileExists(unzip_file)) {
				rc = boinc_zip(UNZIP_IT, unzip_file, basedir);
				#ifdef DEBUG
				fprintf(stderr, "Unzip routine returned %i\n", rc);
				#endif
					
				lockfile = fopen(unzip_lock_finished, "w");
				fclose(lockfile);
				return rc;
        	}
			else {
				fprintf(stderr, "%s does not exist. Aborting...\n", unzip_file);
				remove(unzip_lock_started);
				return -1;
	        }
    	}
    }
    else
    	return 0;
}

void handle_pyerror(const char *errormsg){
	PyObject *ptype, *pvalue, *ptraceback;
	PyErr_Fetch(&ptype, &pvalue, &ptraceback);
	print_error("python", errormsg, PyString_AsString(ptype));
	if (pvalue != NULL)
		print_error("python", "pvalue:\n", PyString_AsString(pvalue));
	if (ptraceback != NULL)
		print_error("python", "Traceback:\n", PyString_AsString(ptraceback));
	Py_XDECREF(ptype);
	Py_XDECREF(pvalue);
	Py_XDECREF(ptraceback);
	Py_Finalize();
}

static PyObject *logFunction(PyObject *self, PyObject *args){
	char *logstring;
	if (!PyArg_ParseTuple(args, "s", &logstring)){
		handle_pyerror("Parsing args tuble in log function");
		return 0;
	}
	int chars_written = fprintf(stderr, "%s\n", logstring);
	return Py_BuildValue("i", chars_written);

}

static PyObject *exitFunction(PyObject *self, PyObject *args){
	char *logstring;
	if (!PyArg_ParseTuple(args, "s", &logstring)){
		handle_pyerror("Parsing args tuble in log function");
		return 0;
	}
	print_error("py_exit", "loader.exit called: ", logstring);
}

static PyMethodDef loader_methods[] = {
	{"log", logFunction, METH_VARARGS, "" },
	{"exit", exitFunction, METH_VARARGS, "" },
	{NULL, NULL, 0, NULL}
};

#ifdef DEBUG
char* getPythonTraceback()
{
	// Python equivilant:
	// import traceback, sys
	// return "".join(traceback.format_exception(sys.exc_type,
	//    sys.exc_value, sys.exc_traceback))

	PyObject *type, *value, *traceback;
	PyObject *tracebackModule;
	char *chrRetval;

	PyErr_Fetch(&type, &value, &traceback);

	tracebackModule = PyImport_ImportModule("traceback");
	if (tracebackModule != NULL)
	{
		PyObject *tbList, *emptyString, *strRetval;

		tbList = PyObject_CallMethod(
			tracebackModule,
			"format_exception",
			"OOO",
			type,
			value == NULL ? Py_None : value,
			traceback == NULL ? Py_None : traceback);

		emptyString = PyString_FromString("");
		strRetval = PyObject_CallMethod(emptyString, "join",
			"O", tbList);

		chrRetval = strdup(PyString_AsString(strRetval));

		Py_DECREF(tbList);
		Py_DECREF(emptyString);
		Py_DECREF(strRetval);
		Py_DECREF(tracebackModule);
	}
	else
	{
		chrRetval = strdup("Unable to import traceback module.");
	}

	Py_DECREF(type);
	Py_XDECREF(value);
	Py_XDECREF(traceback);

	return chrRetval;
}
#endif

int call_python(char *ProgramName, char *python_path) {
	PyObject *module, *loader, *calc;
	char **argv;

	// Reserve space for resolved filenames
	char *resolved_files[4];
	char file1[MAX_PATH] = { 0 };
	char file2[MAX_PATH] = { 0 };
	char file3[MAX_PATH] = { 0 };
	resolved_files[0] = file1;
	resolved_files[1] = file2;
	resolved_files[2] = file3;

	int rc;
	double ram = 0;

	argv = &ProgramName;
	// Prevent import of site.py
	Py_NoSiteFlag = 1;

	add_numpy_builtin();
  	PyImport_FrozenModules = _PyImport_FrozenModules;
  	Py_FrozenFlag = 1;
  	Py_NoSiteFlag = 1;

	Py_SetProgramName(ProgramName);

	Py_SetPythonHome(python_path);
	Py_InitializeEx(0);

	loader = Py_InitModule("loader", loader_methods);

	if(loader == NULL){
#ifdef DEBUG
		char *traceback = getPythonTraceback();
		print_error("python", "importing loader module", traceback);
		free(traceback);
		Py_Finalize();
#else
		handle_pyerror("importing loader module:\n");
#endif
		return -1;
	}

	// Setup Python sys path to include modeller Module
	PySys_SetArgv(1, argv);

	module = PyImport_ImportModule("linker_generator");
	if(module == NULL){
#ifdef DEBUG
		char *traceback = getPythonTraceback();
		print_error("python", "importing linker_generator module", traceback);
		free(traceback);
		Py_Finalize();
#else
		handle_pyerror("importing Module\n");
#endif
		//handle_pyerror("importing Module\n");
		//PyErr_Print();
		return -2;
	}

	calc = PyObject_GetAttrString(module, "calc");
	if(calc == NULL){
#ifdef DEBUG
		char *traceback = getPythonTraceback();
		print_error("python", "getting calc attribute", traceback);
		free(traceback);
		Py_Finalize();
#else
		handle_pyerror("getting calc attribute\n");
#endif
		return -3;
	}
	int i;
	// Don't pass first file
	for (i = 0; i < 2; i++){
		rc = boinc_resolve_filename(used_filenames[i], resolved_files[i], MAX_PATH);
		if (rc){
			print_error("resolve files", "unable to resolve filename: ", used_filenames[i]);
			Py_Finalize();
			return -4;
		}
	}

	rc = boinc_resolve_filename(OUTPUT_FILE, resolved_files[2], MAX_PATH);
	if (rc){
		print_error("resolve files", "unable to resolve filename: ", OUTPUT_FILE);
		Py_Finalize();
		return -5;
	}
	//ram = boinc_get_ram(0.45);
	ram = 1; // Set RAM to 1Gb
	fprintf(stderr, "Calling Python Method with Arguments: %s, %s, %s, %f\n", resolved_files[0], resolved_files[1], resolved_files[2], ram);

	PyObject_CallFunction(calc, "sssd", resolved_files[0], resolved_files[1], resolved_files[2], ram);

	if (PyErr_Occurred()!= NULL){
#ifdef DEBUG
		char *traceback = getPythonTraceback();
		print_error("python", "in modeller execution", traceback);
		free(traceback);
		Py_Finalize();
#else
		handle_pyerror("in modeller execution\n");
#endif
		return -6;
	}

	Py_Finalize();
	return 0;
}


int delete_old_versions(char *basedir){
	struct dirent *entry;
	printf("Searching in dir %s for old versions\n", basedir);
	char delete_path[MAX_PATH] = { 0 };
	DIR *dir = { 0 };
	int app_major = 0, app_minor = 0;

	dir = opendir(basedir);
	if (dir == NULL){
		print_error("delete_old_files", "unable to open basedir", NULL);
		return -1;
	}
	while (entry = readdir(dir)){
		if (entry->d_type == DT_DIR){
			app_major = app_minor = 0;
			int count = sscanf(entry->d_name, APP_NAME"-%i.%i", &app_major, &app_minor);
			if (count == 2)
			{	
				if (app_major <= APP_MAJOR){
					if (app_minor < APP_MINOR || app_major < APP_MAJOR){
						strncpy(delete_path, basedir, MAX_PATH);
						strncat(delete_path, DIR_SLASH, MAX_PATH - strlen(delete_path));
						strncat(delete_path, entry->d_name, MAX_PATH - strlen(delete_path));
						print_error("delete_old_files", "Deleting: ", delete_path);
						DeleteDirectory(delete_path, 1);
					}
				}

			}
		}
	}
	closedir(dir);
	return 0;
}

int main(int argc, char **argv)
{
	char modeller_path[MAX_PATH] = {0};
	char exec_path[MAX_PATH] = { 0 };
	int rc;

	boinc_init();

	// Get Paths and append version path

	getExecPath(modeller_path, MAX_PATH, 1);
	delete_old_versions(modeller_path);
	strncat(modeller_path, DIR_SLASH APP_NAME "-"APP_VERSION, MAX_PATH-strlen(modeller_path));

	getExecPath(exec_path, MAX_PATH, 0);
	strncat(exec_path, DIR_SLASH APP_NAME "-"APP_VERSION, MAX_PATH - strlen(exec_path));


	rc = unzip_resources(modeller_path); // Extract modeller resource files if not existent
	if(rc < 0){
		printf("Temporary exiting...\n");
		boinc_temporary_exit_wrapper(20, "Something went wrong while extracting!\n", FALSE);
	}

	//if(check_file_signings((const char **)used_filenames))
	//	{print_error("main", "verifying files failed!", NULL); boinc_finish(-1);}

	/*	Set Modeller environment variables correctly 						**
	**	and reset PYTHONPATH to avoid loading of locally installed files 	*/


	#ifdef __linux
	//setenv("LIBS_LIB9v14", libs_path, 1);
	setenv("PYTHONPATH", modeller_path, 1);
	//setenv("NUITKA_IMPORT_PATH", pylib_path, 1);
	#elif _WIN32
	//_putenv_s("LIBS_LIB9v14", libs_path);
	_putenv_s("PYTHONPATH", modeller_path);
	//_putenv_s("NUITKA_IMPORT_PATH", pylib_path);
	#endif
	
	#ifdef DEBUG
//		printf("Setting Modeller install dir: \"%s\"\nSetting path to libs file: \"%s\"\n", modeller_path, libs_path);
	#endif

	rc = call_python(exec_path, modeller_path);
	if (rc < 0){
		print_error("main", "somthing went wrong in python, exiting with errorcode!", NULL);
		fprintf(stderr, "Python Errorcode: %i\n", rc);
		boinc_finish(-1);
	}
	boinc_finish(0);
	return -1;
}
