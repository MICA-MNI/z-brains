import os
import subprocess
from .constants import ProcessingException
import fnmatch
import os
import shutil
import sys
import tempfile
from contextlib import contextmanager
import xml.etree.ElementTree as ET
import re


def remove_doctype_from_xml(file_path):
    # Read the XML file as a string
    with open(file_path, "r") as file:
        xml_string = file.read()

    # Remove the !DOCTYPE line
    xml_string = re.sub(r"<!DOCTYPE[^>]*>", "", xml_string)

    # Parse the modified string back into an XML tree
    tree = ET.ElementTree(ET.fromstring(xml_string))

    # Write the changes back to the file
    tree.write(file_path)


def edit_root_attributes_in_xml(file_path, new_attributes):
    # Parse the XML file
    tree = ET.parse(file_path)
    root = tree.getroot()

    # Update the root's attributes
    root.attrib.update(new_attributes)

    # Write the changes back to the file
    tree.write(file_path)


def replace_field_in_xml(file_path, field_path, new_value):
    # Parse the XML file
    tree = ET.parse(file_path)
    root = tree.getroot()

    # Find the field
    field = root.find(field_path)

    # If the field doesn't exist, raise an error
    if field is None:
        raise ValueError(f"No field found at path '{field_path}'")

    # Replace the field's text with the new value
    field.text = new_value

    # Write the changes back to the file
    tree.write(file_path)


def add_field_to_xml(file_path, parent_field_path, field_name, field_value):
    # Parse the XML file
    tree = ET.parse(file_path)
    root = tree.getroot()

    # Find the parent element
    parent_element = root.find(parent_field_path)

    # If the parent element doesn't exist, raise an error
    if parent_element is None:
        raise ValueError(f"No element found at path '{parent_field_path}'")

    # Create a new element
    new_element = ET.Element(field_name)
    new_element.text = field_value

    # Add the new element to the parent element
    parent_element.append(new_element)

    # Write the changes back to the file
    tree.write(file_path)


def delete_temp_folders(directory):
    for root, dirs, files in os.walk(directory, topdown=False):
        for name in dirs:
            if fnmatch.fnmatch(name, "z_brains_temp.????????"):
                print(f"Deleting temp folder {name}")
                shutil.rmtree(os.path.join(root, name))


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
    return output_variable, args[len(output_variable) + 1 :]


def assert_required(option, value, error_message=None):
    if value is None:
        if error_message is None:
            error_message = f"{option} is required"
        show_error(error_message)
        raise ProcessingException(error_message)


def assert_same_size(option1, list1, option2, list2):
    if len(list1) != len(list2):
        show_error(f"{option1} and {option2} must have the same number of elements")
        raise ProcessingException(
            f"{option1} and {option2} must have the same number of elements"
        )


def assert_exists(path, error_message=None):
    if not os.path.exists(path):
        if error_message is None:
            error_message = f"{path} does not exist"
        show_error(error_message)
        raise ProcessingException(error_message)


def assert_columns_in_csv(csv, required_columns):
    with open(csv, "r") as f:
        header = f.readline().strip().split(",")
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


@contextmanager
def tempdir(SUBJECT_OUTPUT_DIR, prefix):
    path = tempfile.mkdtemp(dir=SUBJECT_OUTPUT_DIR, prefix=prefix)
    try:
        yield path
    finally:
        print(f"Cleaning up temp dir {path}")
        try:
            shutil.rmtree(path)
        except IOError:
            sys.stderr.write(f"Failed to clean up temp dir {path}")


def print_version():
    print("Version 1.0")
