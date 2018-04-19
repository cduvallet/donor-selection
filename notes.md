Some notes on how we can start thinking about selecting donors for FMT trials.

If there are no existing studies on the indication of interest, then there are two choices:

1. Go with logistics (timeline, amount of material, etc)   
2. Find a "similar" disease and try to extrapolate    

# Defining healthy donors

But if we have data available, we can start thinking about strategically choosing donors.

I think we can approach this from an entire community perspective, or from an OTU-level perspective.

## Healthy donors

If we think about the entire community, we can do something like a build a ML classifier to distinguish healthy vs. disease (or whatever the differential diagnosis of interest is), then use that classifier to "score" each of our potential donors, and pick the one with the highest score.

You might also be able to do something similarity-based, like the donor which is most similar to an average "healthy" person in your target demographic (i.e. children, Rwandans, vegans, etc... Not necessarily just vanilla HMP data, but the "native" "healthy" population of interest.)

## Healthy OTUs

I think you could also build a scoring function that's built on individual OTUs. Either by including "healthy" OTUs, and/or excluding "disease" OTUs.

Some ideas on defining healthy **OTUs**, from the Hsiao et al cholera dataset: (1) correlation with distance to healthy community [i.e. corr(abundance of OTU, mean distance of that community to other healthy communities)], (2) indicator species analysis/differential abundance, and (3) prominent members of adult microbiome with known beneficial functions (i.e. SCFA).

More generally:
1. bugs that, when they're present, indicate a community that's "closer" to healthy   
2. bugs that are important: differential testing, indicator analysis, classifiers   
3. bugs that are known to be health-associated: functions, my non-specifically healthy bugs

Also, clearly, any "core" healthy bugs (i.e. bactera that are health-associated across many different studies).

Then, once you have your list of "healthy" OTUs, you can score each donor based on some combination of each OTU's "importance" and its "probability of engraftment" (probably abundance).

You could also do similar things to define "disease" OTUs (though this might be harder to tease apart? I dunno) and incorporate those into your scoring function as negatives.

## Other random thoughts

Something about pulling from ML/AI world: train one model to optimize for something (i.e. abundance of "good" bugs) and then train another model to reward/punish that first model's answer based on something else (i.e. similarity to "good"/distance from "bad" communities)?