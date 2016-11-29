from bauhaus import Workflow
from .datasetOps import *
from .mapping import genChunkedMapping
from bauhaus.experiment import (InputType, ResequencingConditionTable)

__all__ = ["NoHQMappingWorkflow"]

class NoHQMappingWorkflow(Workflow):
    @staticmethod
    def name():
        return "NoHQMapping"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        outputDict = {}
        for condition in ct.conditions:
            with pflow.context("condition", condition):
                reference = ct.reference(condition)
                if ct.inputType == InputType.SubreadSet:
                    inputs = ct.inputs(condition)
                    noHQAlnSets = genChunkedMapping(pflow, inputs, reference, extraBlasrArgs="--ignoreHQRegions")
                    outputDict[condition] = noHQAlnSets
                else:
                    raise NotImplementedError, "Support not yet implemented for this input type"
        return outputDict
