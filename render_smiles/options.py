"""
Configuration options for molecular rendering.

This module defines the Options dataclass that controls various aspects
of molecular diagram generation using RDKit.
"""

from dataclasses import dataclass


@dataclass
class Options:
    """
    Configuration options for molecular rendering.

    This class contains all configurable parameters for generating molecular diagrams,
    including dimensions, colors, fonts, and various display options.
    """

    width: int = 500
    height: int = 500
    backgroundColour: tuple = (1, 1, 1)
    dummiesAreAttachments: bool = False
    dummyIsotopeLabels: bool = False
    noAtomLabels: bool = False
    addBondIndices: bool = False
    addAtomIndices: bool = False
    unspecifiedStereoIsUnknown: bool = False
    useMolBlockWedging: bool = False
    baseFontSize: float = 0.6
    fontFile: str | None = None
    additionalAtomLabelPadding: float = 0.1
    bondLineWidth: float = 2
    multipleBondOffset: float = 0.15
    rotate: int = 30
    remove_non_molecular_paths: bool = True
