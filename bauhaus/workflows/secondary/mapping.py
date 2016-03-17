__all__ = [ "genMappingWorkflow", "genChunkedMappingWorkflow" ]

import os.path as op

# ----- Splitting, merging, consolidating -----

# N.B.: The name for subreads sets generated by dataset split is
#       currently {movieName}.subreadset.{chunk#}.xml, which Martin
#       and I think is wrong -- it should be
#       {movieName}.chunk.{subreadset#}.xml; Martin will fix this
#       soon.  When it is fixed, we don't need separate methods for
#       subreads and alignments.

def subreadsChunkName(datumFilename):
    """
    returns movieName or movieName.chunk#, depending on whehter chunked or not
    """
    base = op.basename(datumFilename)
    fields = base.split(".")
    assert len(fields) in (3, 4) and fields[-1] == "xml"
    movieName = fields[0]
    assert fields[-2].startswith("chunk")
    return "%s.%s" % (movieName, fields[-2])

def alignmentsChunkName(datumFilename):
    base = op.basename(datumFilename)
    fields = base.split(".")
    assert len(fields) in (3, 4) and fields[-1] == "xml"
    movieName = fields[0]
    assert fields[1].startswith("chunk")
    return "%s.%s" % (movieName, fields[1])

def movieName(filename):
    base = op.basename(filename)
    fields = base.split(".")
    assert len(fields) in (3, 4) and fields[-1] == "xml"
    movieName = fields[0]
    return movieName



# -------- Generators ------------

# Conventions (in the absence of types!):
#
#   - Every generator should returns the list of context-resolved
#     outputs from its final build statement
#

# ------- Split/merge -------------


def genSubreadsSetSplit(pflow, subreadsSet, splitFactor):
    # Split by ZMWs.  Returns: [dset]
    assert splitFactor >= 1
    pflow.genRuleOnce(
            "splitByZmw",
            "$grid dataset split --zmws --chunks %d --outdir $outdir $in" % (splitFactor,))
    movie = movieName(subreadsSet)
    splitOutputs =  [ "{condition}/subreads_chunks/%s.chunk%d.subreadset.xml" % (movie, i)
                      for i in xrange(splitFactor) ]
    buildStmt = pflow.genBuildStatement(splitOutputs,
                                        "splitByZmw",
                                        [subreadsSet],
                                        variables={"outdir": "{condition}/subreads_chunks"})
    return buildStmt.outputs


def genAlignmentSetMerge(pflow, alignmentSets):
    # MERGE entails the lightweight operation of merging the
    # dataset XML files--BAM files not merged
    pflow.genRuleOnce("mergeAlignmentSets",
                      "$grid dataset merge $out $in")
    outputs = [ "{condition}/mapping/{movieName}_preconsolidate.alignmentset.xml" ]
    buildStmt = pflow.genBuildStatement(outputs,
                                        "mergeAlignmentSets",
                                        alignmentSets)
    return buildStmt.outputs

def genAlignmentSetConsolidate(pflow, alignmentSets):
    # CONSOLIDATE entails actually merging BAM files
    # We first need to run a MERGE.
    mergedAlignmentSet = genAlignmentSetMerge(pflow, alignmentSets)[0]
    pflow.genRuleOnce("consolidateAlignmentSets",
                      "$grid dataset consolidate $in $out")
    consolidateAlignmentSet = "{condition}/mapping/{movieName}.alignmentset.xml"
    consolidatedBam = "{condition}/mapping/{movieName}.aligned_subreads.bam"
    outputs = [ consolidatedBam, consolidateAlignmentSet ]
    buildStmt = pflow.genBuildStatement(outputs,
                                        "consolidateAlignmentSets",
                                        [mergedAlignmentSet])
    return buildStmt.outputs[1]

# ----------- Mapping --------------------


# Assumption that should be asserted: each input subreads set comes from a separate movie.

def genMapping(pflow, subreadsSets, reference):
    """
    Map the subreads set, without chunking.  This is painfully slow on
    Sequel-scale data.
    """
    mapRule = pflow.genRuleOnce(
        "map",
        "$gridSMP $ncpus pbalign --nproc $ncpus $in $reference $out")
    # TODO: these need to be "merged" within each condition
    for subreadsSet in subreadsSets:
        with pflow.context("movieName", movieName(subreadsSet)):
            buildVariables = dict(reference=reference, ncpus=8)
            pflow.genBuildStatement(["{condition}/mapping/{movieName}.alignmentset.xml"],
                                    "map",
                                    [subreadsSet],
                                    buildVariables)

def genChunkedMapping(pflow, subreadsSets, reference, splitFactor=8):
    """
    Break the subreads set into chunks, map the chunks, then
    consolidate the mapped chunks
    """
    mapRule = pflow.genRuleOnce(
        "map",
        "$gridSMP $ncpus pbalign --nproc $ncpus $in $reference $out")
    # TODO: need to "merge" within each condition
    for subreadsSet in subreadsSets:
        with pflow.context("movieName", movieName(subreadsSet)):
            alignmentSetChunks = []
            subreadsSetChunks = genSubreadsSetSplit(pflow, subreadsSet, splitFactor)
            for (i, subreadsSetChunk) in enumerate(subreadsSetChunks):
                with pflow.context("chunkNum", i):
                    buildVariables = dict(reference=reference, ncpus=8)
                    buildStmt = pflow.genBuildStatement(
                        ["{condition}/mapping_chunks/{movieName}.chunk{chunkNum}.alignmentset.xml"],
                        "map",
                        [subreadsSetChunk],
                        buildVariables)
                    alignmentSetChunks.extend(buildStmt.outputs)
            genAlignmentSetConsolidate(pflow, alignmentSetChunks)



# ---------- Workflows -------------

def genMappingWorkflow(pflow, ct):
    for condition in ct.conditions:
        with pflow.context("condition", condition):
            subreadsSets = ct.inputs(condition)
            reference = ct.reference(condition)
            genMapping(pflow, subreadsSets, reference)


def genChunkedMappingWorkflow(pflow, ct):
    for condition in ct.conditions:
        with pflow.context("condition", condition):
            subreadsSets = ct.inputs(condition)
            reference = ct.reference(condition)
            genChunkedMapping(pflow, subreadsSets, reference, splitFactor=8)




# -------------------- Demo -------------------

# lambda short insert ecoli runs intended for CCS analysis
if __name__ == '__main__':
    from bauhaus.pflow import PFlow
    inputDataByCondition = { "Replicate1" : ["/pbi/collections/315/3150128/r54008_20160308_001811/1_A01/m54008_160308_002050.subreadset.xml"],
                             "Replicate2" : ["/pbi/collections/315/3150128/r54008_20160308_001811/2_B01/m54008_160308_053311.subreadset.xml"] }

    reference = "/mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta"

    pflow = PFlow()
    for (condition, inputSubreadSets) in inputDataByCondition.iteritems():
        with pflow.context("condition", condition):
            #genMapping(pflow, inputSubreadSets, reference)
            genChunkedMapping(pflow, inputSubreadSets, reference, 8)
    pflow.write("build.ninja")
