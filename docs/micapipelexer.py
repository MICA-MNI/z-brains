# -*- coding: utf-8 -*-
"""
    pygments.styles.micapipe
    ~~~~~~~~~~~~~~~~~~~~~

    A modern style based on the VIM pyte theme.


"""

from pygments.style import Style
from pygments.token import Keyword, Name, Comment, String, Error, \
    Number, Operator, Whitespace, Generic


class micapipeLexerStyle(Style):
    """
    A modern style based on the VIM pyte theme.
    """

    background_color = "#d0d0fe"
    default_style = ""

    styles = {
        Whitespace:                "#bbbbbb",
        Comment:                   "italic #888888",
        Comment.Preproc:           "noitalic #5a2ca0",
        Comment.Special:           "noitalic bg:#e3d7f4",

        Keyword:                   "bold #5a2ca0",
        Keyword.Pseudo:            "nobold",
        Keyword.Type:              "nobold #902000",

        Operator:                  "#000000",
        Operator.Word:             "bold #5a2ca0",

        Name.Builtin:              "#5a2ca0",
        Name.Function:             "#06287e",
        Name.Class:                "bold #9966cc",
        Name.Namespace:            "bold #9966cc",
        Name.Exception:            "#5a2ca0",
        Name.Variable:             "#bb60d5",
        Name.Constant:             "#60add5",
        Name.Label:                "bold #002070",
        Name.Entity:               "bold #d55537",
        Name.Attribute:            "#4140a0",
        Name.Tag:                  "bold #062873",
        Name.Decorator:            "bold #555555",

        String:                    "#4140a0",
        String.Doc:                "italic",
        String.Interpol:           "italic #70a0d0",
        String.Escape:             "bold #4140a0",
        String.Regex:              "#235388",
        String.Symbol:             "#517918",
        String.Other:              "#EF4405",
        Number:                    "#149c58",

        Generic.Heading:           "bold #000080",
        Generic.Subheading:        "bold #800080",
        Generic.Deleted:           "#A00000",
        Generic.Inserted:          "#00A000",
        Generic.Error:             "#FF0000",
        Generic.Emph:              "italic",
        Generic.Strong:            "bold",
        Generic.Prompt:            "bold #EF4405",
        Generic.Output:            "#888",
        Generic.Traceback:         "#04D",

        Error:                     "border:#FF0000"
    }
