import os
import subprocess
import sys
from functions.new_constants import ProcessingException

def do_cmd(*args):
    
    if len(args) == 1:
        print(args)
        array = args[0].split()
    else:
        array = args

    str_cmd = ""
    for element in array:
        element = str(element)
        if " " in element:
            element = f'"{element}"'
        str_cmd += element + " "
    str_cmd = str_cmd.rstrip()
    print(str_cmd)
    print(f"COMMAND --> {str_cmd}")
    subprocess.run(array)

def log_message(level, *messages):
    if int(os.getenv("VERBOSE", "-1")) >= level:
        print(*messages)

def show_error(*messages):
    print("\033[38;5;9mERROR:\033[0m", *messages)

def show_warning(*messages):
    print("\033[38;5;184mWARNING:\033[0m", *messages)

def show_note(*args):
    print("\033[0;36;10mNOTE:\033[0m", *args)

def show_info(*messages):
    print("\033[38;5;75mINFO:\033[0m", *messages)

def show_title(*messages):
    print("\033[38;5;141m", *messages, "\033[0m")

def allowed_to_regex(array):
    return "|".join(array)

def parse_option_single_value(output_variable, args, allowed_values=None):
    if allowed_values is not None:
        if args[1] not in allowed_values:
            show_error(f"Invalid value for {args[0]}: {args[1]}")
            raise ProcessingException(f"Invalid value for {args[0]}: {args[1]}")
    output_variable = args[1]
    return output_variable, args[2:]

def parse_option_multiple_values(output_variable, args, allowed_values=None, all=None):
    if allowed_values is not None:
        for value in args[1:]:
            if value not in allowed_values:
                show_error(f"Invalid value for {args[0]}: {value}")
                raise ProcessingException(f"Invalid value for {args[0]}: {value}")
    if all is not None and "all" in args[1:]:
        output_variable = all
    else:
        output_variable = args[1:]
    return output_variable, args[len(output_variable)+1:]

def assert_required(option, value, error_message=None):
    if value is None:
        if error_message is None:
            error_message = f"{option} is required"
        show_error(error_message)
        raise ProcessingException(error_message)

def assert_same_size(option1, list1, option2, list2):
    if len(list1) != len(list2):
        show_error(f"{option1} and {option2} must have the same number of elements")
        raise ProcessingException(f"{option1} and {option2} must have the same number of elements")

def assert_exists(path, error_message=None):
    if not os.path.exists(path):
        if error_message is None:
            error_message = f"{path} does not exist"
        show_error(error_message)
        raise ProcessingException(error_message)

def assert_columns_in_csv(csv, required_columns):
    with open(csv, 'r') as f:
        header = f.readline().strip().split(',')
    for column in required_columns:
        if column not in header:
            show_error(f"{csv} is missing column {column}")
            raise ProcessingException(f"{csv} is missing column {column}")

def submit_job(scheduler, *args):
    if scheduler == "qsub":
        do_cmd("qsub", *args)
    elif scheduler == "sbatch":
        do_cmd("sbatch", *args)
    else:
        show_error(f"Unknown scheduler: {scheduler}")
        raise ProcessingException(f"Unknown scheduler: {scheduler}")

def print_version():
    print("Version 1.0")