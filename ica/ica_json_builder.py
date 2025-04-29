#/usr/bin/env python3

import json
from typing import Any
import yaml # type: ignore
import sys

def interpret_string(in_type: str, default: str | None = None) -> dict[str, str]:
    """
    Interpret a string type and return a dictionary with the appropriate ICA attributes.
    Arguments:
        in_type: The type of the input.
        default: The default value for the input.
    Returns:
        A dictionary with the ICA attributes for the input type.
    """
    out: dict[str, Any] = {}
    if default:
        out["value"] = default
    if in_type.endswith("?"):
        in_type = in_type.replace("?", "")
    else:
        out["minValues"] = 1
    if "[]" in in_type:
        out["maxValues"] = 1000
        in_type = in_type.replace("[]", "")
    converter: dict[str, str] = {
        "string": "textbox",
        "float": "number",
        "int": "integer",
        "File": "data",
        "Directory": "data"
    }
    if in_type == "null":
        out["minValues"] = 0
    elif in_type == "boolean":
        if "value" in out: del out["value"]
        out.update(interpret_bool(default))
    elif in_type in converter:
        out["type"] = converter[in_type]
    else:
        raise ValueError(f"Unaccounted type: {in_type}")
    if in_type == "File":
        out["dataFilter"] = {"dataType": "file"}
    elif in_type == "Directory":
        out["dataFilter"] = {"dataType": "directory"}
    return out


def interpret_enum(enum_type: dict[str, str], default: str | None = None) -> dict[str, str]:
    """
    Interpret an enum type and return a dictionary with the appropriate ICA attributes.
    Arguments:
        enum_type: A dictionary with the enum type.
        default: The default value for the enum type.
    Returns:
        A dictionary with the ICA attributes for the enum type.
    """
    out: dict[str, Any] = {"type": "select",
           "choices": [{"value": i, "text": i, "selected": i == default} for i in enum_type["symbols"]]}
    return out

def interpret_bool(default: str | None = None) -> dict[str, str]:
    """
    Interpret a boolean type and return a dictionary with the appropriate ICA attributes.
    Arguments:
        default: The default value for the boolean type.
    Returns:
        A dictionary with the ICA attributes for the boolean type.
    """
    out = {"type": "select",
           "choices": [{"value": i, "text": str(i), "selected": i == default} for i in [True, False]]}
    return out

def interpret_list(list_type: list[str], default: str | None = None) -> dict[str, str]:
    """
    Interpret a list type and return a dictionary with the appropriate ICA attributes.
    Arguments:
        list_type: A list with the types of the input.
        default: The default value for the list type.
    Returns:
        A dictionary with the ICA attributes for the list type.
    """
    out = {"minValues": 1}
    for t in list_type:
        if type(t) == str:
            out.update(interpret_string(t))
        elif type(t) == dict:
            out.update(interpret_enum(t, default))
        else:
            raise Exception(f"I don't know what to do with {t}")
    return out


def main():
    with open(sys.argv[-1], "r") as file:
        wf = yaml.safe_load(file)
    inputs = [dict(v, **{"id": k}) for k, v in wf["inputs"].items()]
    payload = {"fields":[]}
    for input in inputs:
        ica_attr = {"id": input["id"]}
        doc = None if "doc" not in input else input["doc"]
        if doc: ica_attr["helpText"] = doc
        default = None if "default" not in input else input["default"]
        if type(input["type"]) == dict:
            ica_attr.update(interpret_enum(input["type"], default))
        elif type(input["type"]) == str:
            ica_attr.update(interpret_string(input["type"], default))
        elif type(input["type"]) == list:
            ica_attr.update(interpret_list(input["type"], default))
        else:
            raise Exception(f"Cannot process input {input["type"]}")
        payload["fields"].append(ica_attr)
    print(json.dumps(payload, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
